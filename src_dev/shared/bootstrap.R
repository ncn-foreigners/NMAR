#' Shared bootstrap variance helpers
#' @description S3 generic + methods to estimate the variance of an estimator
#'   via resampling (IID) or replicate weights (survey). Designed to be reused
#'   across NMAR engines.
#' @details
#'   - For `data.frame` inputs, performs i.i.d. bootstrap by resampling rows and
#'     rerunning `estimator_func`.
#'   - For `survey.design` inputs, converts to a bootstrap replicate‑weight
#'     design (`svrep::as_bootstrap_design`) and uses `survey::withReplicates` to
#'     compute replicate estimates; variance is computed with `survey::svrVar`.
#'   `estimator_func` is typically an engine method (e.g., `el()`), and is called
#'   with the same arguments used for the point estimate, except that the `data`
#'   argument is replaced by the resampled data or replicate design.
#' @param data a `data.frame` or a `survey.design`.
#' @param estimator_func function that returns an S3 result object; the primary
#'   estimate is extracted via `estimate()` and convergence via `$converged`.
#' @param point_estimate numeric; point estimate used for some survey variance formulas.
#' @param ... passed through to `estimator_func`.
#' @keywords internal
bootstrap_variance <- function(data, estimator_func, point_estimate, ...) {
  UseMethod("bootstrap_variance")
}

#' Bootstrap for i.i.d. data.frames
#' @inheritParams bootstrap_variance
#' @param bootstrap_reps integer; number of resamples.
#' @return a list with `se`, `variance`, and the vector of `replicates`.
#' @exportS3Method
bootstrap_variance.data.frame <- function(data, estimator_func, point_estimate, bootstrap_reps = 500, ...) {
  n_obs <- nrow(data)
  estimates <- numeric(bootstrap_reps)
  dot_args <- list(...)

  for (i in seq_len(bootstrap_reps)) {
    resample_indices <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
    bootstrap_data <- data[resample_indices, , drop = FALSE]

    call_args <- c(list(data = bootstrap_data, on_failure = "return"), dot_args)

    fit <- suppressWarnings({
      do.call(estimator_func, call_args)
    })

    estimates[i] <- if (!is.null(fit$converged) && !fit$converged) NA else as.numeric(estimate(fit))
  }

  failed_reps <- sum(is.na(estimates))
  if (failed_reps > 0 && (failed_reps / bootstrap_reps > 0.1)) {
    warning(sprintf(
      "%d out of %d bootstrap replicates failed to converge. The variance estimate may be unreliable.",
      failed_reps, bootstrap_reps
    ), call. = FALSE)
  }

  boot_variance <- stats::var(estimates, na.rm = TRUE)
  list(
    se = sqrt(boot_variance),
    variance = boot_variance,
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
#' @return a list with `se`, `variance`, and the vector of `replicates`.
#' @exportS3Method
bootstrap_variance.survey.design <- function(data, estimator_func, point_estimate, bootstrap_reps = 500, ...) {
  if (!requireNamespace("svrep", quietly = TRUE)) {
    stop("Package 'svrep' is required for bootstrap variance with survey objects. Please install it.", call. = FALSE)
  }

  dot_args <- list(...)
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

  rep_design <- do.call(
    svrep::as_bootstrap_design,
    c(list(design = data, replicates = bootstrap_reps), bootstrap_settings)
  )

  estimator_args <- dot_args
  template_call <- nmar_extract_svydesign_call(data)

  replicate_results <- survey::withReplicates(
    design = rep_design,
    theta = function(pw, data_subset) {
      data_subset$..replicate_weights.. <- pw
      temp_design <- nmar_reconstruct_design(template_call, data_subset)
      call_args <- c(list(data = temp_design), estimator_args)
      fn_formals <- tryCatch(names(formals(estimator_func)), error = function(e) character())
      if ("on_failure" %in% fn_formals && is.null(call_args$on_failure)) {
        call_args$on_failure <- "return"
      }
      fit <- suppressWarnings({
        do.call(estimator_func, call_args)
      })
      if (!is.null(fit$converged) && !fit$converged) return(NA_real_)
      as.numeric(estimate(fit))
    },
    return.replicates = TRUE
  )

  replicate_estimates <- replicate_results$replicates
  if (anyNA(replicate_estimates)) {
    failed_reps <- sum(is.na(replicate_estimates))
    stop(sprintf(
      "%d out of %d bootstrap replicates failed to produce an estimate (NA). Consider adjusting solver settings or bootstrap configuration.",
      failed_reps, bootstrap_reps
    ), call. = FALSE)
  }

  boot_variance <- survey::svrVar(
    thetas = replicate_estimates,
    scale = rep_design$scale,
    rscales = rep_design$rscales,
    mse = rep_design$mse,
    coef = point_estimate
  )
  list(
    se = as.numeric(sqrt(boot_variance)),
    variance = as.numeric(boot_variance),
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
  call_copy <- template_call
  call_copy$data <- quote(data_subset)
  call_copy$weights <- as.formula(paste0("~", weight_var))
  eval(call_copy, envir = list(data_subset = data_subset), enclos = parent.frame())
}
