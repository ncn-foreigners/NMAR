#' Shared bootstrap variance helpers
#' @description Internal helpers to estimate the variance of a scalar estimator
#'   via bootstrap resampling (IID data) or bootstrap replicate weights
#'   (survey designs). Designed to be reused across NMAR engines.
#' @details
#'   \itemize{
#'     \item For \code{data.frame} inputs, performs IID bootstrap by resampling
#'       rows and rerunning \code{estimator_func} on each resample, then
#'       computing the empirical variance of the replicate estimates.
#'     \item For \code{survey.design} inputs, converts the design to a bootstrap
#'       replicate-weight design with \code{svrep::as_bootstrap_design()},
#'       evaluates \code{estimator_func} on each replicate weight vector (by
#'       injecting the replicate analysis weights into a copy of the input design), and
#'       passes the resulting replicate estimates and replicate scaling factors
#'       to \code{survey::svrVar()}.
#'   }
#'
#'   \code{estimator_func} is typically an engine-level estimator (for example
#'   the EL engine) and is called with the same arguments used for the point
#'   estimate, except that the \code{data} argument is replaced by the resampled
#'   data (IID) or a replicate-weighted \code{survey.design} (survey). Arguments
#'   reserved for the bootstrap implementation are stripped from \code{...}
#'   before forwarding.
#'
#' @section Bootstrap-specific options:
#'   \describe{
#'     \item{\code{resample_guard}}{IID bootstrap only. A function
#'       \code{function(indices, data)} that returns \code{TRUE} to accept a
#'       resample and \code{FALSE} to reject it.}
#'     \item{\code{bootstrap_settings}}{Survey bootstrap only. A list of
#'       arguments forwarded to \code{svrep::as_bootstrap_design()}.}
#'     \item{\code{bootstrap_options}}{Alias for \code{bootstrap_settings}.}
#'     \item{\code{bootstrap_type}}{Shortcut for the \code{type} argument to
#'       \code{svrep::as_bootstrap_design()}.}
#'     \item{\code{bootstrap_mse}}{Shortcut for the \code{mse} argument to
#'       \code{svrep::as_bootstrap_design()}.}
#'   }
#'
#' @section Progress Reporting:
#'   If the optional \code{progressr} package is installed, bootstrap calls
#'   signal progress via a \code{progressr::progressor} inside
#'   \code{progressr::with_progress()}. Users control whether progress is shown
#'   (and how) by registering handlers with \code{progressr::handlers()}. When
#'   \code{progressr} is not installed or no handlers are active, bootstrap runs
#'   silently. Progress reporting is compatible with all future backends.
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
#'   If the optional \code{future.apply} package is installed, bootstrap uses
#'   \code{future.apply::future_lapply(future.seed = TRUE)} which provides
#'   backend-independent, parallel-safe random number streams under the
#'   \code{future} framework. If \code{future.apply} is not installed, bootstrap
#'   falls back to sequential evaluation via \code{base::lapply()}, which is
#'   still reproducible under \code{set.seed()} but may not match the
#'   \code{future.seed} stream.
#'
#' @param data A \code{data.frame} or a \code{survey.design}.
#' @param estimator_func Function returning an object with a numeric scalar
#'   component \code{y_hat} and an optional logical component \code{converged}.
#' @param point_estimate Numeric scalar; used for survey bootstrap variance
#'   (passed to \code{survey::svrVar()} as \code{coef}).
#' @param ... Additional arguments. Some are consumed by \code{bootstrap_variance()}
#'   itself (for example \code{resample_guard} for IID bootstrap or
#'   \code{bootstrap_settings}/\code{bootstrap_options}/\code{bootstrap_type}/\code{bootstrap_mse}
#'   for survey bootstrap); remaining arguments are forwarded to \code{estimator_func}.
#' @keywords internal
bootstrap_variance <- function(data, estimator_func, point_estimate, ...) {
# Check for replicate designs first (they do not inherit from survey.design).
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

# Internal helper: detect whether future.apply is available.
nmar_has_future_apply <- function() {
  requireNamespace("future.apply", quietly = TRUE)
}

# Internal helper: warn once per R session that we are falling back to sequential
# evaluation because future.apply is not installed.
nmar_warn_no_future_apply_once <- function() {
  opt <- "NMAR.bootstrap.warned_no_future_apply"
  if (isTRUE(getOption(opt, FALSE))) return(invisible(FALSE))
  warning(
    "Package 'future.apply' is not installed. Running bootstrap sequentially via base::lapply().\n  ",
    "Install 'future.apply' to enable future-based parallel execution and future-seeded RNG (future.seed = TRUE).",
    call. = FALSE,
    immediate. = TRUE
  )
  options(setNames(list(TRUE), opt))
  invisible(TRUE)
}

# Internal helper: warn when survey bootstrap assumptions may be violated.
nmar_warn_survey_bootstrap_assumptions <- function(design) {
  if (!inherits(design, "survey.design")) return(invisible(FALSE))

  allprob <- design$allprob
  has_multistage_probs <- is.data.frame(allprob) && ncol(allprob) > 1L
  has_pps <- isTRUE(design$pps)
  has_fpc_arg <- FALSE
  has_fpc_popsize <- FALSE
  has_probs_arg <- FALSE
  has_pps_arg <- FALSE

  fpc <- design$fpc
  has_fpc_popsize <- is.list(fpc) && !is.null(fpc$popsize)

  dc <- try(getCall(design), silent = TRUE)
  if (!inherits(dc, "try-error") && !is.null(dc)) {
    args <- as.list(dc)[-1]
    has_fpc_arg <- !is.null(args$fpc) || !is.null(args$fpctype)
    has_probs_arg <- !is.null(args$probs)
    has_pps_arg <- !is.null(args$pps)
  }

  if (has_multistage_probs || has_pps || has_fpc_arg || has_fpc_popsize || has_probs_arg || has_pps_arg) {
    warning(
      "Survey bootstrap injects replicate analysis weights into the design.\n  ",
      "This is valid when the estimator depends on weights (and optionally strata/cluster)\n  ",
      "but does not recompute stage-specific probabilities or FPC. This design appears\n  ",
      "to include PPS/multistage probabilities or FPC; if the estimator uses those\n  ",
      "fields directly, bootstrap variance may be incorrect.",
      call. = FALSE
    )
    return(invisible(TRUE))
  }
  invisible(FALSE)
}
# Internal helper: apply over X using future.apply (if installed) or base::lapply
# (sequential fallback). Progress is reported via progressr if installed.
nmar_bootstrap_apply <- function(X, FUN, use_progress, future_globals = NULL, future_packages = NULL) {
  use_future <- nmar_has_future_apply()
  if (!use_future) nmar_warn_no_future_apply_once()

  if (isTRUE(use_progress) && requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p <- progressr::progressor(steps = length(X))
      wrapper <- function(x) {
        res <- FUN(x)
        p()
        res
      }
      if (use_future) {
        fg <- future_globals
        if (is.null(fg)) {
          fg <- TRUE
        } else if (isTRUE(fg)) {
          fg <- TRUE
        } else if (is.list(fg)) {
          fg <- c(fg, list(FUN = FUN, p = p))
        }
        future.apply::future_lapply(
          X,
          wrapper,
          future.seed = TRUE,
          future.globals = fg,
          future.packages = future_packages
        )
      } else {
        lapply(X, wrapper)
      }
    })
  } else {
    if (use_future) {
      if (is.null(future_globals)) future_globals <- TRUE
      future.apply::future_lapply(
        X,
        FUN,
        future.seed = TRUE,
        future.globals = future_globals,
        future.packages = future_packages
      )
    } else {
      lapply(X, FUN)
    }
  }
}

#' Default method dispatch (internal safety net)
#' @keywords internal
bootstrap_variance.default <- function(data, estimator_func, point_estimate, ...) {
  stop("Unsupported data type for bootstrap_variance().", call. = FALSE)
}

#' Bootstrap for IID data frames
#' @inheritParams bootstrap_variance
#' @param data A \code{data.frame}.
#' @param point_estimate Unused for IID bootstrap; included for signature
#'   consistency.
#' @param bootstrap_reps integer; number of resamples.
#' @return A list with components \code{se}, \code{variance}, and \code{replicates}.
#' @keywords internal
bootstrap_variance.data.frame <- function(data, estimator_func, point_estimate, bootstrap_reps = 500, ...) {
  validator_assert_positive_integer(bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  if (bootstrap_reps < 2) {
    stop("`bootstrap_reps` must be at least 2 for variance estimation.", call. = FALSE)
  }

  n_obs <- nrow(data)

  if (n_obs == 0) {
    stop("Cannot bootstrap from empty data (0 rows).", call. = FALSE)
  }

  dot_args <- list(...)
  est_fun <- estimator_func
  resample_guard <- NULL
  if (!is.null(dot_args$bootstrap_cores)) dot_args$bootstrap_cores <- NULL
  if (!is.null(dot_args$bootstrap_workers)) dot_args$bootstrap_workers <- NULL
  if (!is.null(dot_args$bootstrap_settings)) dot_args$bootstrap_settings <- NULL
  if (!is.null(dot_args$bootstrap_options)) dot_args$bootstrap_options <- NULL
  if (!is.null(dot_args$bootstrap_type)) dot_args$bootstrap_type <- NULL
  if (!is.null(dot_args$bootstrap_mse)) dot_args$bootstrap_mse <- NULL
  if (!is.null(dot_args$survey_na_policy)) dot_args$survey_na_policy <- NULL
  if (!is.null(dot_args$resample_guard)) {
# Some NMAR estimators require each resample to contain at least one respondent.
# Allow callers to supply a guard that rejects unsuitable resamples.
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

# Cap attempts to avoid infinite loops when guards are too strict.
        if (attempts >= 20) {
          return(NA_real_)
        }
      } else {
        break
      }
    }

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

  use_progress <- requireNamespace("progressr", quietly = TRUE)
  lst <- nmar_bootstrap_apply(
    X = seq_len(bootstrap_reps),
    FUN = replicate_fn,
    use_progress = use_progress
  )

  estimates <- vapply(lst, identity, numeric(1))

  failed_reps <- sum(is.na(estimates))
  if (failed_reps > 0 && failed_reps < bootstrap_reps && (failed_reps / bootstrap_reps > 0.1)) {
    warning(sprintf(
      "%d out of %d bootstrap replicates failed to produce finite estimates. The variance estimate may be unreliable.",
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
#' @param data A \code{survey.design}.
#' @param bootstrap_reps integer; number of bootstrap replicates.
#' @details This path constructs a replicate-weight design using
#'   \code{svrep::as_bootstrap_design()} and evaluates the estimator on each set of
#'   bootstrap replicate analysis weights.
#'
#'   Replicate evaluation starts from a shallow template copy of the input survey
#'   design (including its ids/strata/fpc structure) and injects each replicate's
#'   analysis weights by
#'   updating the design's probability slots (\code{prob}/\code{allprob}) so that
#'   \code{weights(design)} returns the desired replicate weights (with
#'   zero weights represented as \code{prob = Inf}). This avoids replaying or
#'   reconstructing a \code{svydesign()} call and therefore supports designs
#'   created via \code{subset()} and \code{update()}.
#'
#'   \strong{NA policy:} By default, survey bootstrap uses a strict NA policy:
#'   if any replicate fails to produce a finite estimate, the entire bootstrap
#'   fails with an error. Setting \code{survey_na_policy = "omit"} drops failed
#'   replicates (and their corresponding \code{rscales}) and proceeds with the
#'   remaining replicates.
#'
#' @section Limitations:
#'   \strong{Calibrated/post-stratified designs:} Post-hoc adjustments applied
#'   via \code{survey::calibrate()}, \code{survey::postStratify()}, or
#'   \code{survey::rake()} are not supported here and will cause the function to
#'   error. These adjustments are not recomputed when replicate weights are
#'   injected, so the replicate designs would not reflect the intended
#'   calibrated/post-stratified analysis.
#'
#' @param survey_na_policy Character string specifying how to handle replicates
#'   that fail to produce estimates. Options:
#'   \describe{
#'     \item{\code{"strict"}}{(default) Any failed replicate causes an error.
#'       This is a conservative default that makes instability explicit.}
#'     \item{\code{"omit"}}{Failed replicates are omitted. The corresponding
#'       \code{rscales} are also omitted to maintain correct variance scaling.
#'       Use with caution: if failures are non-random, variance may be biased.}
#'   }
#' @return A list with components \code{se}, \code{variance}, and \code{replicates}.
#' @keywords internal
bootstrap_variance.survey.design <- function(data, estimator_func, point_estimate, bootstrap_reps = 500,
                                             survey_na_policy = c("strict", "omit"), ...) {
  if (!requireNamespace("svrep", quietly = TRUE)) {
    stop("Package 'svrep' is required for bootstrap variance with survey objects. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for bootstrap variance with survey objects. Please install it.", call. = FALSE)
  }

  survey_na_policy <- match.arg(survey_na_policy)

# Guard against replicate designs (svyrep.design). Bootstrapping a replicate
# design would create a second-order bootstrap and is not supported here.
  if (inherits(data, "svyrep.design")) {
    stop(
      "Cannot bootstrap a replicate design (svyrep.design).\n  ",
      "Bootstrap variance requires a standard survey.design object.\n  ",
      "The input design already contains replicate weights.\n  ",
      "If you need bootstrap variance, start from the original survey.design.",
      call. = FALSE
    )
  }

  validator_assert_positive_integer(bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  if (bootstrap_reps < 2) {
    stop("`bootstrap_reps` must be at least 2 for variance estimation.", call. = FALSE)
  }

  if (nrow(data$variables) == 0) {
    stop("Cannot bootstrap from empty survey design (0 observations).", call. = FALSE)
  }

  dot_args <- list(...)
  est_fun <- estimator_func
  if (!is.null(dot_args$bootstrap_cores)) dot_args$bootstrap_cores <- NULL
  if (!is.null(dot_args$bootstrap_workers)) dot_args$bootstrap_workers <- NULL
  if (!is.null(dot_args$resample_guard)) dot_args$resample_guard <- NULL
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

  if (!is.null(bootstrap_settings$design) || !is.null(bootstrap_settings$replicates)) {
    stop(
      "`bootstrap_settings` must not include `design` or `replicates`.\n  ",
      "Use `data = <survey.design>` and `bootstrap_reps = <int>` instead.",
      call. = FALSE
    )
  }

  rep_design <- do.call(svrep::as_bootstrap_design, c(list(design = data, replicates = bootstrap_reps), bootstrap_settings))

  estimator_args <- dot_args

# Detect calibrated or post-stratified designs. These adjustments are tied to
# the original analysis weights and are not recomputed when we inject replicate
# weights, so bootstrap variance would be incorrect.
  if (!is.null(data$postStrata) || !is.null(data$calibration)) {
    adjustments <- c(
      if (!is.null(data$postStrata)) "post-stratification" else NULL,
      if (!is.null(data$calibration)) "calibration" else NULL
    )
    stop(
      "Bootstrap variance is not currently supported for calibrated or post-stratified designs.\n",
      "  Detected: ", paste(adjustments, collapse = " and "), "\n",
      "  These adjustments are not recomputed when replicate weights are injected.\n",
      "  The replicate designs would be unadjusted, leading to incorrect variance.",
      call. = FALSE
    )
  }

  nmar_warn_survey_bootstrap_assumptions(data)

# Extract replicate analysis weights matrix (one column per replicate).
  repw <- weights(rep_design, type = "analysis")

# Check replicate count (may differ from the requested number in stratified
# designs). The variance formula is valid for the actual count produced.
  J_actual <- if (is.null(dim(repw))) 1L else ncol(repw)
  if (!is.finite(J_actual) || J_actual < 2L) {
    stop(
      sprintf(
        "Bootstrap replicate design produced %d replicate(s). Variance estimation requires at least 2 replicates.\n  ",
        J_actual
      ),
      "This can happen with very small strata or restrictive design settings. ",
      "Consider increasing per-stratum sample size or adjusting bootstrap settings.",
      call. = FALSE
    )
  }
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

# Save variance scaling factors before freeing rep_design; these are needed
# when calling survey::svrVar().
  rep_scale <- rep_design$scale
  rep_rscales <- rep_design$rscales
  rep_mse <- rep_design$mse
  if (isTRUE(rep_mse) && (!is.numeric(point_estimate) || length(point_estimate) != 1L || !is.finite(point_estimate))) {
    stop(
      "`point_estimate` must be a finite numeric scalar when bootstrap replicate design uses mse = TRUE.",
      call. = FALSE
    )
  }

# Extract data frame only (not full survey design) to reduce serialization.
  data_vars <- data$variables

# Keep a shallow design template for per-replicate weight injection.
# Drop variables to avoid duplicating the (potentially large) data.frame.
  design_template <- data
  design_template$variables <- NULL

  rm(rep_design)

# Replicates are indexed to avoid materializing a list of weight vectors and to
# keep the set of exported globals small and predictable under future.apply.
  replicate_eval <- function(j) {
    pw <- repw[, j]
    temp_design <- nmar_inject_design_weights(
      template_design = design_template,
      variables = data_vars,
      analysis_weights = pw
    )
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

  use_progress <- requireNamespace("progressr", quietly = TRUE)

  J <- seq_len(J_actual)

  lst <- nmar_bootstrap_apply(
    X = J,
    FUN = replicate_eval,
    use_progress = use_progress,
    future_globals = list(
      replicate_eval = replicate_eval,
      repw = repw,
      data_vars = data_vars,
      nmar_inject_design_weights = nmar_inject_design_weights,
      design_template = design_template,
      estimator_args = estimator_args,
      est_fun = est_fun
    ),
    future_packages = "survey"
  )

  replicate_estimates <- vapply(lst, identity, numeric(1))

  if (anyNA(replicate_estimates)) {
    total_reps <- length(replicate_estimates)
    failed_idx <- which(!is.finite(replicate_estimates))
    n_failed <- length(failed_idx)

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
      stop(sprintf(
        paste0(
          "%d/%d survey bootstrap replicates failed to produce finite estimates.\n  %s\n\n  ",
          "Survey bootstrap uses strict NA policy (survey_na_policy='strict').\n\n  ",
          "Troubleshooting:\n  ",
          "  - Check if estimator converges reliably on full design\n  ",
          "  - If your engine supports a non-bootstrap variance method, consider switching to it\n  ",
          "  - Increase solver tolerances (control = list(xtol = 1e-6))\n  ",
          "  - Reduce bootstrap_reps if stratification creates small replicates\n  ",
          "  - Or set survey_na_policy = 'omit' to allow some failures"
        ),
        n_failed, total_reps, pattern_msg
      ), call. = FALSE)
    } else if (survey_na_policy == "omit") {
      keep_idx <- which(is.finite(replicate_estimates))

      if (length(keep_idx) < 2) {
        stop(sprintf(
          paste0(
            "Too few successful survey bootstrap replicates (%d/%d succeeded).\n  ",
            "At least 2 finite estimates required for variance calculation.\n  ",
            "Consider: reducing bootstrap_reps, improving estimator stability,\n  ",
            "or using a non-bootstrap variance method (if available in your engine)."
          ),
          length(keep_idx), total_reps
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
        n_failed, total_reps, pattern_msg, length(keep_idx)
      ), call. = FALSE, immediate. = TRUE)

# When omitting replicates, subset rscales to the same indices so that
# survey::svrVar() applies the correct replicate-specific scaling.
      replicate_estimates <- replicate_estimates[keep_idx]
      rep_rscales <- rep_rscales[keep_idx]
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

# Internal helper: inject analysis weights into a survey.design by updating
# prob/allprob so that weights(design) returns the desired weights (with
# zero weights represented as prob = Inf).
nmar_inject_design_weights <- function(template_design, variables, analysis_weights) {
  if (!inherits(template_design, "survey.design")) {
    stop("Internal error: template_design must be a survey.design.", call. = FALSE)
  }
  if (!is.data.frame(variables)) {
    stop("Internal error: variables must be a data.frame.", call. = FALSE)
  }
  n <- nrow(variables)
  if (length(analysis_weights) != n) {
    stop("Internal error: replicate weights must align with design variables.", call. = FALSE)
  }
  analysis_weights <- as.numeric(analysis_weights)
  if (any(!is.finite(analysis_weights)) || any(analysis_weights < 0)) {
    stop("Internal error: replicate weights must be finite and nonnegative.", call. = FALSE)
  }

  prob <- rep(Inf, n)
  pos <- analysis_weights > 0
  prob[pos] <- 1 / analysis_weights[pos]

  tmp <- template_design
  tmp$variables <- variables
  tmp$prob <- prob
# Keep allprob consistent with prob (weights = 1/prob).
  tmp$allprob <- data.frame(prob = prob)
  tmp
}
