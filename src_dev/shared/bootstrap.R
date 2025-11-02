#' Shared bootstrap variance helpers
#' @description S3 generic + methods to estimate the variance of an estimator
#'   via resampling (IID) or replicate weights (survey). Designed to be reused
#'   across NMAR engines.
#' @details
#'   - For `data.frame` inputs, performs i.i.d. bootstrap by resampling rows and
#'     rerunning `estimator_func`.
#'   - For `survey.design` inputs, converts to a bootstrap replicate-weight
#'     design (`svrep::as_bootstrap_design`), computes replicate estimates by
#'     rebuilding the original design with each replicate weight vector, and
#'     then computes variance with `survey::svrVar` using the replicate scales.
#'   `estimator_func` is typically an engine method (e.g., `el()`), and is called
#'   with the same arguments used for the point estimate, except that the `data`
#'   argument is replaced by the resampled data or replicate design.
#'
#' @section Progress Reporting:
#'   If the \code{progressr} package is installed, progress reporting is available.
#'   Enable it by setting handlers before calling the bootstrap:
#'
#'   \code{library(progressr)}
#'
#'   \code{handlers(global = TRUE)}
#'
#'   \code{handlers("txtprogressbar")  # or "progress", "cli", etc.}
#'
#'   To disable progress in simulations or batch jobs, use \code{handlers("void")}.
#'   If progressr is not installed or no handlers are set, bootstrap runs silently
#'   (default behavior). Progress reporting works with all future backends (sequential,
#'   multisession, cluster, etc.) and does not affect reproducibility.
#'
#' @section Reproducibility:
#'   For reproducible bootstrap results, always set a seed before calling
#'   the estimation function:
#'
#'   \preformatted{
#'   set.seed(123)  # Set seed for reproducibility
#'   result <- nmar(Y ~ X, data = df,
#'                  engine = el_engine(variance_method = "bootstrap",
#'                                     bootstrap_reps = 500))
#'   }
#'
#'   The \code{future} package (via \code{future.seed = TRUE}) ensures each
#'   bootstrap replicate uses an independent L'Ecuyer-CMRG random number stream
#'   derived from this seed, guaranteeing reproducibility across all future
#'   backends (sequential, multisession, cluster, etc.).
#'
#' @param data a `data.frame` or a `survey.design`.
#' @param estimator_func function that returns an S3 result object; the primary
#'   estimate is extracted via `$y_hat` and convergence via `$converged`.
#' @param point_estimate numeric; point estimate used for some survey variance formulas.
#' @param ... passed through to `estimator_func`.
#' @keywords internal
bootstrap_variance <- function(data, estimator_func, point_estimate, ...) {
# Check for replicate designs first (they don't inherit from survey.design)
  if (inherits(data, "svyrep.design")) {
    stop(
      "Cannot bootstrap a replicate design (svyrep.design).\n  ",
      "Bootstrap variance requires a standard survey.design object.\n  ",
      "The input design already contains replicate weights.\n  ",
      "If you need bootstrap variance, start from the original survey.design.",
      call. = FALSE
    )
  }
  if (inherits(data, "survey.design")) {
    return(bootstrap_variance.survey.design(data, estimator_func, point_estimate, ...))
  }
  if (is.data.frame(data)) {
    return(bootstrap_variance.data.frame(data, estimator_func, point_estimate, ...))
  }
  stop("Unsupported data type for bootstrap_variance().", call. = FALSE)
}

#' Default method dispatch (internal safety net)
#' @keywords internal
bootstrap_variance.default <- function(data, estimator_func, point_estimate, ...) {
  stop("Unsupported data type for bootstrap_variance().", call. = FALSE)
}

#' Bootstrap for i.i.d. data.frames
#' @inheritParams bootstrap_variance
#' @param bootstrap_reps integer; number of resamples.
#' @return a list with `se`, `variance`, and the vector of `replicates`.
#' @keywords internal
bootstrap_variance.data.frame <- function(data, estimator_func, point_estimate, bootstrap_reps = 500, ...) {
# Validate bootstrap_reps
  validator$assert_positive_integer(bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  if (bootstrap_reps < 2) {
    stop("`bootstrap_reps` must be at least 2 for variance estimation.", call. = FALSE)
  }

  n_obs <- nrow(data)

# Validate data is non-empty
  if (n_obs == 0) {
    stop("Cannot bootstrap from empty data (0 rows).", call. = FALSE)
  }

  dot_args <- list(...)
  est_fun <- estimator_func
  resample_guard <- NULL
  if (!is.null(dot_args$resample_guard)) {
# Some estimators (exptilt) require the bootstrap replicate to contain at
# least one respondent. Allow callers to supply a simple guard to reject
# unsuitable resamples
    resample_guard <- dot_args$resample_guard
    dot_args$resample_guard <- NULL
  }

  replicate_fn <- function(i) {
    attempts <- 0
    repeat {
      resample_indices <- sample.int(n_obs, n_obs, replace = TRUE)
      attempts <- attempts + 1

      if (!is.null(resample_guard)) {
        guard_ok <- tryCatch(
          resample_guard(resample_indices, data = data),
          error = function(e) FALSE
        )
        if (isTRUE(guard_ok)) break

# Give up after max attempts - return NA immediately
        if (attempts >= 20) {
          return(NA_real_)
        }
      } else {
        break
      }
    }

# Continue with valid resample
    bootstrap_data <- data[resample_indices, , drop = FALSE]
    call_args <- c(list(data = bootstrap_data), dot_args)
    fn_formals <- tryCatch(names(formals(est_fun)), error = function(e) character())
    if ("on_failure" %in% fn_formals && is.null(call_args$on_failure)) {
      call_args$on_failure <- "return"
    }
    fit <- tryCatch(suppressWarnings({ do.call(est_fun, call_args) }), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    if (!is.null(fit$converged) && !fit$converged) return(NA_real_)
    val <- tryCatch(as.numeric(fit$y_hat), error = function(e) NA_real_)
    if (length(val) != 1L || !is.finite(val)) return(NA_real_)
    val
  }

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required for bootstrap variance.", call. = FALSE)
  }

# Use progressr if available for progress reporting (optional)
  use_progress <- requireNamespace("progressr", quietly = TRUE)

  if (use_progress) {
# Wrap in progressr for optional progress reporting
    lst <- progressr::with_progress({
      p <- progressr::progressor(steps = bootstrap_reps)
      future.apply::future_lapply(
        seq_len(bootstrap_reps),
        function(i) {
          res <- replicate_fn(i)
          p() # Signal progress
          res
        },
        future.seed = TRUE
      )
    })
  } else {
# No progressr available - use standard future_lapply
    lst <- future.apply::future_lapply(
      seq_len(bootstrap_reps),
      replicate_fn,
      future.seed = TRUE
    )
  }

  estimates <- vapply(lst, identity, numeric(1))

  failed_reps <- sum(is.na(estimates))
  if (failed_reps > 0 && (failed_reps / bootstrap_reps > 0.1)) {
    warning(sprintf(
      "%d out of %d bootstrap replicates failed to converge. The variance estimate may be unreliable.",
      failed_reps, bootstrap_reps
    ), call. = FALSE)
  }
  if (failed_reps == bootstrap_reps) {
    stop("All bootstrap replicates failed.", call. = FALSE)
  }

  good_n <- sum(is.finite(estimates), na.rm = TRUE)
  if (!is.finite(good_n) || good_n < 2L) {
    stop(sprintf(
      "Too few successful bootstrap replicates (%d/%d). Increase bootstrap_reps or check estimator stability.",
      good_n, bootstrap_reps
    ), call. = FALSE)
  }
  boot_variance <- stats::var(estimates, na.rm = TRUE)
  var_num <- as.numeric(boot_variance)
  if (!is.finite(var_num)) {
    var_num <- NA_real_
  }
  list(
    se = if (is.finite(var_num)) sqrt(pmax(var_num, 0)) else NA_real_,
    variance = var_num,
    replicates = estimates
  )
}

#' Bootstrap for survey designs via replicate weights
#' @inheritParams bootstrap_variance
#' @param bootstrap_reps integer; number of bootstrap replicates.
#' @details This path constructs a replicate-weight design using
#'   [svrep::as_bootstrap_design()] and rebuilds the original sampling design for
#'   each replicate weight vector. The supplied design must have been created
#'   directly with [survey::svydesign()].
#'
#'   \strong{NA Policy:} Survey bootstrap uses a strict NA policy - if any replicate
#'   fails to produce a finite estimate, the entire bootstrap fails with an error.
#'   This ensures the replicate design structure is maintained for design-calibrated
#'   variance via \code{survey::svrVar()}. In contrast, IID bootstrap allows up to
#'   10\% failures before warning, as it uses uncalibrated \code{stats::var()}.
#'
#' @section Limitations:
#'   \strong{Design Reconstruction:} Survey bootstrap currently supports only designs
#'   created directly with \code{survey::svydesign()}. Post-hoc adjustments applied
#'   via \code{survey::calibrate()}, \code{survey::postStratify()}, or
#'   \code{survey::rake()} cannot be reconstructed across bootstrap replicates and
#'   will cause the function to error.
#'
#'   Calibrated or post-stratified designs are not supported by this bootstrap
#'   path. Start from the original `survey::svydesign()` object prior to
#'   calibration/post-stratification.
#'
#'   \strong{Supported Design Features:} The following \code{svydesign()} parameters
#'   are preserved during reconstruction:
#'   \itemize{
#'     \item \code{ids} (or \code{id}): Sampling unit identifiers
#'     \item \code{strata}: Stratification variables
#'     \item \code{fpc}: Finite population correction
#'     \item \code{nest}: Nested vs non-nested strata
#'   }
#'
#'   The following are NOT preserved (they conflict with replicate weights):
#'   \itemize{
#'     \item \code{probs}: Sampling probabilities (incompatible with direct weights)
#'     \item \code{pps}: PPS sampling specification (incompatible with direct weights)
#'   }
#'
#'   \strong{Rationale:} When reconstructing designs for each replicate, we replace
#'   the original weights with bootstrap replicate weights. Specifying both
#'   \code{weights} and \code{probs}/\code{pps} simultaneously is undefined behavior
#'   in \code{survey::svydesign()}. The structural design parameters (ids, strata,
#'   fpc, nest) define the sampling topology and are preserved; the weights define
#'   the analysis weights and are replaced.
#'
#' @param survey_na_policy Character string specifying how to handle replicates that
#'   fail to produce estimates. Options:
#'   \describe{
#'     \item{\code{"strict"}}{(default) Any failed replicate causes an error.
#'       This ensures the full replicate design structure is maintained and
#'       is required for proper calibration-based variance.}
#'     \item{\code{"omit"}}{Failed replicates are omitted. The corresponding
#'       \code{rscales} are also omitted to maintain correct variance scaling.
#'       Use with caution: if failures are non-random, variance may be biased.}
#'   }
#' @return a list with `se`, `variance`, and the vector of `replicates`.
#' @keywords internal
bootstrap_variance.survey.design <- function(data, estimator_func, point_estimate, bootstrap_reps = 500,
                                             survey_na_policy = c("strict", "omit"), ...) {
  if (!requireNamespace("svrep", quietly = TRUE)) {
    stop("Package 'svrep' is required for bootstrap variance with survey objects. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for bootstrap variance with survey objects. Please install it.", call. = FALSE)
  }

# Validate survey_na_policy argument
  survey_na_policy <- match.arg(survey_na_policy)

# Guard against replicate designs (svyrep.design)
# Bootstrapping a replicate design is ill-defined (creates second-order replicates)
  if (inherits(data, "svyrep.design")) {
    stop(
      "Cannot bootstrap a replicate design (svyrep.design).\n  ",
      "Bootstrap variance requires a standard survey.design object.\n  ",
      "The input design already contains replicate weights.\n  ",
      "If you need bootstrap variance, start from the original survey.design.",
      call. = FALSE
    )
  }

# Validate bootstrap_reps
  validator$assert_positive_integer(bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  if (bootstrap_reps < 2) {
    stop("`bootstrap_reps` must be at least 2 for variance estimation.", call. = FALSE)
  }

# Validate survey is non-empty
  if (nrow(data$variables) == 0) {
    stop("Cannot bootstrap from empty survey design (0 observations).", call. = FALSE)
  }

  dot_args <- list(...)
  est_fun <- estimator_func
  bootstrap_settings <- list()
  if (!is.null(dot_args$bootstrap_settings)) {
    if (!is.list(dot_args$bootstrap_settings)) {
      stop("`bootstrap_settings` must be a list of arguments understood by svrep::as_bootstrap_design().", call. = FALSE)
    }
    bootstrap_settings <- utils::modifyList(bootstrap_settings, dot_args$bootstrap_settings)
    dot_args$bootstrap_settings <- NULL
  }
  if (!is.null(dot_args$bootstrap_options)) {
    if (!is.list(dot_args$bootstrap_options)) {
      stop("`bootstrap_options` must be a list.", call. = FALSE)
    }
    bootstrap_settings <- utils::modifyList(bootstrap_settings, dot_args$bootstrap_options)
    dot_args$bootstrap_options <- NULL
  }
  if (!is.null(dot_args$bootstrap_type)) {
    bootstrap_settings$type <- dot_args$bootstrap_type
    dot_args$bootstrap_type <- NULL
  }
  if (!is.null(dot_args$bootstrap_mse)) {
    bootstrap_settings$mse <- dot_args$bootstrap_mse
    dot_args$bootstrap_mse <- NULL
  }

  rep_design <- do.call(svrep::as_bootstrap_design, c(list(design = data, replicates = bootstrap_reps), bootstrap_settings))

  estimator_args <- dot_args
  template_call <- nmar_extract_svydesign_call(data)

# Detect calibrated or post-stratified designs
# These adjustments cannot be reconstructed across bootstrap replicates
# because they are stored as internal state, not in the original svydesign() call
  if (!is.null(data$postStrata) || !is.null(data$calibration)) {
    adjustments <- c(
      if (!is.null(data$postStrata)) "post-stratification" else NULL,
      if (!is.null(data$calibration)) "calibration" else NULL
    )
    stop(
      "Bootstrap variance is not currently supported for calibrated or post-stratified designs.\n",
      "  Detected: ", paste(adjustments, collapse = " and "), "\n",
      "  These adjustments cannot be reconstructed across bootstrap replicates.\n",
      "  The replicate designs would be uncalibrated, leading to incorrect variance.",
      call. = FALSE
    )
  }

# Extract replicate weights matrix from the replicate design using survey API.
# This returns a matrix of replicate weights (one column per replicate).
# weights() is a generic from stats; S3 dispatch calls weights.svyrep.design() from survey
  repw <- weights(rep_design, type = "analysis", rep = TRUE)

# Check replicate count (may differ from request in stratified designs)
# The variance formula is statistically valid for the actual count produced
  J_actual <- ncol(repw)
  if (J_actual != bootstrap_reps) {
    warning(sprintf(
      paste0(
        "Replicate design produced %d replicates (requested %d).\n  ",
        "This can occur with stratified designs having small strata.\n  ",
        "Proceeding with %d replicates; variance is statistically valid\n  ",
        "but precision may differ from expectations based on %d replicates."
      ),
      J_actual, bootstrap_reps, J_actual, bootstrap_reps
    ), call. = FALSE, immediate. = TRUE)
  }

# Save variance scaling factors before freeing rep_design
# These are needed later for survey::svrVar()
  rep_scale <- rep_design$scale
  rep_rscales <- rep_design$rscales
  rep_mse <- rep_design$mse

# Extract data frame only (not full survey design) to reduce serialization
  data_vars <- data$variables

# Free design object (keep repw matrix and data_vars)
  rm(rep_design)

# Define replicate evaluation function that takes replicate index
# This allows us to iterate by index and export repw once per worker
  replicate_eval <- function(j) {
    pw <- repw[, j]
    data_subset <- data_vars
    data_subset$..replicate_weights.. <- pw
    temp_design <- nmar_reconstruct_design(template_call, data_subset)
    call_args <- c(list(data = temp_design), estimator_args)
    fn_formals <- tryCatch(names(formals(est_fun)), error = function(e) character())
    if ("on_failure" %in% fn_formals && is.null(call_args$on_failure)) {
      call_args$on_failure <- "return"
    }
    fit <- tryCatch(
      suppressWarnings({ do.call(est_fun, call_args) }),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    if (!is.null(fit$converged) && !fit$converged) return(NA_real_)
    val <- tryCatch(as.numeric(fit$y_hat), error = function(e) NA_real_)
    if (length(val) != 1L || !is.finite(val)) return(NA_real_)
    val
  }

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required for bootstrap variance.", call. = FALSE)
  }

# Use progressr if available for progress reporting (optional)
  use_progress <- requireNamespace("progressr", quietly = TRUE)

# Replicate indices (iterate by index, not by pre-split weights)
  J <- seq_len(ncol(repw))

  if (use_progress) {
# Wrap in progressr for optional progress reporting
    lst <- progressr::with_progress({
      p <- progressr::progressor(steps = length(J))
      future.apply::future_lapply(
        J,
        function(j) {
          res <- replicate_eval(j)
          p() # Signal progress
          res
        },
        future.seed = TRUE,
        future.globals = list(
          replicate_eval = replicate_eval,
          repw = repw,
          data_vars = data_vars,
          nmar_reconstruct_design = nmar_reconstruct_design,
          template_call = template_call,
          estimator_args = estimator_args,
          est_fun = est_fun,
          p = p
        ),
        future.packages = "survey"
      )
    })
  } else {
# No progressr available - use standard future_lapply
    lst <- future.apply::future_lapply(
      J,
      replicate_eval,
      future.seed = TRUE,
      future.globals = list(
        replicate_eval = replicate_eval,
        repw = repw,
        data_vars = data_vars,
        nmar_reconstruct_design = nmar_reconstruct_design,
        template_call = template_call,
        estimator_args = estimator_args,
        est_fun = est_fun
      ),
      future.packages = "survey"
    )
  }

  replicate_estimates <- vapply(lst, identity, numeric(1))

# Handle NA replicates according to policy
  if (anyNA(replicate_estimates)) {
    failed_idx <- which(!is.finite(replicate_estimates))
    n_failed <- length(failed_idx)

# Show pattern of failures for debugging
    if (n_failed <= 10) {
      pattern_msg <- sprintf("Replicate indices: %s", paste(failed_idx, collapse = ", "))
    } else {
      pattern_msg <- sprintf(
        "First 10 failed replicates: %s\n  (%d more failures)",
        paste(failed_idx[1:10], collapse = ", "),
        n_failed - 10
      )
    }

    if (survey_na_policy == "strict") {
# Strict policy: any failure is an error
      stop(sprintf(
        paste0(
          "%d/%d survey bootstrap replicates failed to produce finite estimates.\n  %s\n\n  ",
          "Survey bootstrap uses strict NA policy (survey_na_policy='strict').\n\n  ",
          "Troubleshooting:\n  ",
          "  - Check if estimator converges reliably on full design\n  ",
          "  - Try variance_method = 'delta' instead\n  ",
          "  - Increase solver tolerances (control = list(xtol = 1e-6))\n  ",
          "  - Reduce bootstrap_reps if stratification creates small replicates\n  ",
          "  - Or set survey_na_policy = 'omit' to allow some failures"
        ),
        n_failed, bootstrap_reps, pattern_msg
      ), call. = FALSE)
    } else if (survey_na_policy == "omit") {
# Omit policy: subset to successful replicates
      keep_idx <- which(is.finite(replicate_estimates))

      if (length(keep_idx) < 2) {
        stop(sprintf(
          paste0(
            "Too few successful survey bootstrap replicates (%d/%d succeeded).\n  ",
            "At least 2 finite estimates required for variance calculation.\n  ",
            "Consider: reducing bootstrap_reps, improving estimator stability,\n  ",
            "or using variance_method = 'delta' instead."
          ),
          length(keep_idx), length(replicate_estimates)
        ), call. = FALSE)
      }

      warning(sprintf(
        paste0(
          "%d/%d survey bootstrap replicates failed (omitted from variance).\n  ",
          "%s\n  ",
          "Proceeding with %d replicates. Effective sample size reduced.\n  ",
          "Variance estimate may be biased if failures are associated with\n  ",
          "specific design features. Consider investigating failure pattern."
        ),
        n_failed, length(replicate_estimates), pattern_msg, length(keep_idx)
      ), call. = FALSE, immediate. = TRUE)

# Subset replicate_estimates and rep_rscales to match
# This is critical for mathematical correctness
      replicate_estimates <- replicate_estimates[keep_idx]
      rep_rscales <- rep_rscales[keep_idx]
# rep_scale and rep_mse remain unchanged (scalars for bootstrap)
    }
  }

  boot_variance <- survey::svrVar(
    thetas = replicate_estimates,
    scale = rep_scale,
    rscales = rep_rscales,
    mse = rep_mse,
    coef = point_estimate
  )
  var_num <- as.numeric(boot_variance)

  if (!is.finite(var_num)) {
    stop(
      "survey::svrVar() returned non-finite variance (", var_num, "). ",
      "This indicates a problem with the bootstrap design or replicate estimates. ",
      "Ensure your survey design is correctly specified and the estimator ",
      "converges reliably across bootstrap replicates.",
      call. = FALSE
    )
  }

  list(
    se = sqrt(pmax(var_num, 0)),
    variance = var_num,
    replicates = as.vector(replicate_estimates)
  )
}

nmar_extract_svydesign_call <- function(design) {
  design_call <- try(getCall(design), silent = TRUE)
  if (inherits(design_call, "try-error") || is.null(design_call)) {
    stop("Unable to retrieve the original call for the provided survey design. Ensure it was created with survey::svydesign().", call. = FALSE)
  }
  fun <- design_call[[1]]
  is_svydesign <- identical(fun, as.name("svydesign")) || identical(fun, quote(survey::svydesign))
  if (!is_svydesign) {
    stop("bootstrap_variance() currently supports survey designs created directly with survey::svydesign().", call. = FALSE)
  }
  design_call[[1]] <- quote(survey::svydesign)
  design_call
}

nmar_reconstruct_design <- function(template_call, data_subset, weight_var = "..replicate_weights..") {
# Rebuild a fresh survey::svydesign call using the structural pieces from the
# original call, but avoiding eval() in a deep parent frame to prevent stack
# growth under replicate evaluation.
  tc <- template_call
# Extract known args if present
  args <- as.list(tc)[-1]
  get_arg <- function(nm) if (!is.null(args[[nm]])) args[[nm]] else NULL
# Check both 'id' and 'ids' (svydesign accepts both)
# Preserve the original parameter name used in the call
  ids_val <- get_arg("ids")
  id_val <- get_arg("id")
  id_param_name <- if (!is.null(ids_val)) "ids" else if (!is.null(id_val)) "id" else NULL
  id_param_val <- if (!is.null(ids_val)) ids_val else id_val

  strata <- get_arg("strata")
  fpc <- get_arg("fpc")
  nest <- get_arg("nest")
# DO NOT extract probs/pps - they conflict with replicate weights
# probs= and pps= are used to derive initial weights in svydesign()
# We are replacing weights with replicate weights, so probs/pps must be excluded

# Assemble call, including only non-NULL components
  call_list <- list(quote(survey::svydesign))
  if (!is.null(id_param_name)) {
    call_list <- c(call_list, setNames(list(id_param_val), id_param_name))
  }
  if (!is.null(strata)) call_list <- c(call_list, list(strata = strata))
  if (!is.null(fpc)) call_list <- c(call_list, list(fpc = fpc))
  if (!is.null(nest)) call_list <- c(call_list, list(nest = nest))
# DO NOT add probs or pps - they conflict with weights parameter
  call_list <- c(call_list,
    list(data = quote(data_subset), weights = as.formula(paste0("~", weight_var))))
  new_call <- as.call(call_list)
  eval(new_call, envir = list(data_subset = data_subset))
}
