#' Bootstrap variance estimation module
#'
#' @description Estimates the variance of a scalar estimator via bootstrap
#' resampling for IID data or bootstrap replicate weights for survey designs.
#'
#' @details
#' \itemize{
#' \item For \code{data.frame} inputs, performs IID bootstrap by resampling
#' rows and rerunning \code{estimator_func} on each resample, then
#' computing the empirical variance of the replicate estimates.
#' \item For \code{survey.design} inputs, converts the design to a bootstrap
#' replicate-weight design with \code{svrep::as_bootstrap_design()},
#' evaluates \code{estimator_func} on each replicate weight vector by
#' injecting the replicate analysis weights into a copy of the input design,
#' and passes the resulting replicate estimates and replicate scaling factors
#' to \code{survey::svrVar()}.
#' }
#'
#' @section Bootstrap-specific options:
#' \describe{
#' \item{\code{resample_guard}}{IID bootstrap only. A function
#' \code{function(indices, data)} that returns \code{TRUE} to accept a resample
#' and \code{FALSE} to reject it.}
#' \item{\code{bootstrap_settings}}{A list of arguments forwarded to \code{svrep::as_bootstrap_design()}.}
#' \item{\code{bootstrap_options}}{Alias for \code{bootstrap_settings}.}
#' \item{\code{bootstrap_type}}{The \code{type} argument for \code{svrep::as_bootstrap_design()}.}
#' \item{\code{bootstrap_mse}}{The \code{mse} argument for \code{svrep::as_bootstrap_design()}.}
#' }
#'
#' @section Progress Reporting:
#' If the optional \code{progressr} package is installed, bootstrap calls
#' indicate progress via a \code{progressr::progressor} inside
#' \code{progressr::with_progress()}. Users control if and how progress is shown
#' by registering handlers with \code{progressr::handlers()}. When
#' \code{progressr} is not installed or no handlers are active, bootstrap runs
#' silently.
#'
#' @section Parallelization:
#' By default, bootstrap replicate evaluation runs sequentially via
#' \code{base::lapply()} for both IID resampling and survey replicate-weight bootstrap.
#' If the optional \code{future.apply} package is installed, bootstrap can use
#' \code{future.apply::future_lapply(future.seed = TRUE)} when the user has set
#' a parallel \code{future::plan()}.
#' The backend is controlled by the package option \code{nmar.bootstrap_apply}:
#' \describe{
#' \item{\code{"auto"}}{(default) Use \code{base::lapply()} unless the current
#' future plan has more than one worker, in which case use
#' \code{future.apply::future_lapply()} if available.}
#' \item{\code{"base"}}{Always use \code{base::lapply()}, even
#' if \code{future.apply} is installed.}
#' \item{\code{"future"}}{Always use \code{future.apply::future_lapply()}.}
#' }
#' When \code{future.apply} is used, random-number streams are parallel-safe and
#' backend-independent under the \code{future} framework. When \code{base::lapply()}
#' is used, results are reproducible under \code{set.seed()} but will
#' likely not match the \code{future.seed} streams.
#'
#' @param data A \code{data.frame} or a \code{survey.design}.
#' @param estimator_func Function returning an object with a numeric scalar
#' component \code{y_hat} and an optional logical component \code{converged}.
#' @param point_estimate Numeric scalar; used for survey bootstrap variance
#' (passed to \code{survey::svrVar()} as \code{coef}).
#' @param ... Additional arguments. Some are consumed by \code{bootstrap_variance()}
#' itself (for example \code{resample_guard} for IID bootstrap or
#' \code{bootstrap_settings}/\code{bootstrap_options}/\code{bootstrap_type}/\code{bootstrap_mse}
#' or survey bootstrap). Remaining arguments are forwarded to \code{estimator_func}.
#'
#' @keywords internal
bootstrap_variance <- function(data, estimator_func, point_estimate, ...) {
  UseMethod("bootstrap_variance")
}

#' Default dispatch
#' @keywords internal
bootstrap_variance.default <- function(data, estimator_func, point_estimate, ...) {
  stop("Unsupported data type for bootstrap_variance().", call. = FALSE)
}

#' Replicate-weight designs not supported
#' @keywords internal
bootstrap_variance.svyrep.design <- function(data, estimator_func, point_estimate, ...) {
  stop(
    "Cannot bootstrap a replicate design (svyrep.design).\n  ",
    "Bootstrap variance requires a standard survey.design object.\n  ",
    "The input design already contains replicate weights.\n  ",
    "If you need bootstrap variance, start from the original survey.design.",
    call. = FALSE
  )
}

#' Bootstrap for IID data frames
#' @inheritParams bootstrap_variance
#' @param data A \code{data.frame}.
#' @param point_estimate Unused for IID bootstrap, included for signature
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
  if (!is.null(dot_args$bootstrap_settings)) dot_args$bootstrap_settings <- NULL
  if (!is.null(dot_args$bootstrap_options)) dot_args$bootstrap_options <- NULL
  if (!is.null(dot_args$bootstrap_type)) dot_args$bootstrap_type <- NULL
  if (!is.null(dot_args$bootstrap_mse)) dot_args$bootstrap_mse <- NULL
  if (!is.null(dot_args$survey_na_policy)) dot_args$survey_na_policy <- NULL
  if (!is.null(dot_args$resample_guard)) {
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

#' Bootstrap for survey designs
#' @inheritParams bootstrap_variance
#' @param data A \code{survey.design}.
#' @param bootstrap_reps integer; number of bootstrap replicates.
#' @details This path constructs a replicate-weight design using
#' \code{svrep::as_bootstrap_design()} and evaluates the estimator on each set of
#' bootstrap replicate analysis weights.
#' Replicate evaluation starts from a shallow template copy of the input survey
#' design (including its ids/strata/fpc structure) and injects each replicate's
#' analysis weights by updating the design's probability slots (\code{prob}/\code{allprob}) so that
#' \code{weights(design)} returns the desired replicate weights.
#' This avoids replaying or reconstructing a \code{svydesign()} call and therefore
#' supports designs created via \code{subset()} and \code{update()}.
#' \strong{NA policy:} By default, survey bootstrap uses a strict NA policy:
#' if any replicate fails to produce a finite estimate, the entire bootstrap
#' fails with an error. Setting \code{survey_na_policy = "omit"} drops failed
#' replicates and proceeds with the remaining replicates.
#'
#' @section Limitations:
#' \strong{Calibrated/post-stratified designs:} Post-hoc adjustments applied
#' via \code{survey::calibrate()}, \code{survey::postStratify()}, or
#' \code{survey::rake()} are not supported here and will cause the function to
#' error. These adjustments are not recomputed when replicate weights are
#' injected, so the replicate designs would not reflect the intended
#' calibrated/post-stratified analysis.
#'
#' @param survey_na_policy Character string specifying how to handle replicates
#' that fail to produce estimates. Options:
#' \describe{
#' \item{\code{"strict"}}{(default) Any failed replicate causes an error.
#' This is a conservative default that makes instability explicit.}
#' \item{\code{"omit"}}{Failed replicates are omitted. The corresponding
#' \code{rscales} are also omitted to maintain correct variance scaling.
#' Use with caution: if failures are non-random, variance may be biased.}
#' }
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

  validator_assert_positive_integer(bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  if (bootstrap_reps < 2) {
    stop("`bootstrap_reps` must be at least 2 for variance estimation.", call. = FALSE)
  }

  if (nrow(data$variables) == 0) {
    stop("Cannot bootstrap from empty survey design (0 observations).", call. = FALSE)
  }

  dot_args <- list(...)
  est_fun <- estimator_func
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

  repw <- weights(rep_design, type = "analysis")

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

# These will be needed when calling svyVar
  rep_scale <- rep_design$scale
  rep_rscales <- rep_design$rscales
  rep_mse <- rep_design$mse
  if (isTRUE(rep_mse) && (!is.numeric(point_estimate) || length(point_estimate) != 1L || !is.finite(point_estimate))) {
    stop(
      "`point_estimate` must be a finite numeric scalar when bootstrap replicate design uses mse = TRUE.",
      call. = FALSE
    )
  }

# Data frame only instead of full survey design to reduce serialization
  data_vars <- data$variables

# Shallow design template for per-replicate weight injection avoids duplicating
# potentially large data frames
  design_template <- data
  design_template$variables <- NULL

  rm(rep_design)

# Indexing replicates avoids making list of weight vectors
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

# Ensures svyVar applies the correct scaling when omitting replicates
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
  tmp$allprob <- data.frame(prob = prob)

  tmp
}

nmar_has_future_apply <- function() {
  requireNamespace("future.apply", quietly = TRUE)
}

nmar_future_workers <- function() {
  if (!requireNamespace("future", quietly = TRUE)) return(1L)
  w <- tryCatch(future::nbrOfWorkers(), error = function(e) 1L)
  w <- suppressWarnings(as.integer(w))
  if (!is.finite(w) || w < 1L) w <- 1L
  w
}

nmar_bootstrap_apply_backend <- function() {
  opt <- getOption("nmar.bootstrap_apply", "auto")
  if (is.null(opt)) opt <- "auto"
  if (!is.character(opt) || length(opt) != 1L || is.na(opt) || !nzchar(opt)) {
    stop("Option `nmar.bootstrap_apply` must be one of: 'auto', 'base', 'future'.", call. = FALSE)
  }
  opt <- tolower(opt)
  if (!opt %in% c("auto", "base", "future")) {
    stop("Option `nmar.bootstrap_apply` must be one of: 'auto', 'base', 'future'.", call. = FALSE)
  }
  opt
}

nmar_warn_no_future_apply_once <- function(context = NULL) {
  opt <- "NMAR.bootstrap.warned_no_future_apply"
  if (isTRUE(getOption(opt, FALSE))) return(invisible(FALSE))
  msg <- paste0(
    "Package 'future.apply' is not installed. Running bootstrap sequentially via base::lapply().\n  ",
    "Install 'future.apply' to enable future-based parallel execution and future-seeded RNG (future.seed = TRUE)."
  )
  if (is.character(context) && length(context) == 1L && nzchar(context)) {
    msg <- paste0(msg, "\n  Context: ", context)
  }
  warning(msg, call. = FALSE, immediate. = TRUE)
  options(setNames(list(TRUE), opt))
  invisible(TRUE)
}

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

nmar_bootstrap_apply <- function(X, FUN, use_progress, future_globals = NULL, future_packages = NULL) {
  backend <- nmar_bootstrap_apply_backend()

  use_future <- FALSE
  if (backend == "future") {
    if (!nmar_has_future_apply()) {
      stop(
        "Option `nmar.bootstrap_apply = 'future'` requires the suggested package 'future.apply'.",
        call. = FALSE
      )
    }
    use_future <- TRUE
  } else if (backend == "auto") {
    if (nmar_future_workers() > 1L) {
      if (nmar_has_future_apply()) {
        use_future <- TRUE
      } else {
        nmar_warn_no_future_apply_once(context = "future plan has >1 worker")
      }
    }
  } else if (backend == "base") {
    use_future <- FALSE
  }

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
