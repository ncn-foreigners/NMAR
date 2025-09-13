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
#' @param estimator_func function that returns a fit with `y_hat` and `converged`.
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

    estimates[i] <- if (!is.null(fit$converged) && !fit$converged) NA else fit$y_hat
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
#' @details The replicate‑weight path tries to reconstruct a temporary design for
#'   each replicate weight vector; it first attempts to replay the original
#'   `svydesign` call with updated weights and then falls back to simple designs
#'   if needed. Replicates that fail to converge are recorded as `NA`; a warning
#'   is issued if many fail.
#' @return a list with `se`, `variance`, and the vector of `replicates`.
#' @exportS3Method
bootstrap_variance.survey.design <- function(data, estimator_func, point_estimate, bootstrap_reps = 500, ...) {
  if (!requireNamespace("svrep", quietly = TRUE)) {
    stop("Package 'svrep' is required for bootstrap variance with survey objects. Please install it.", call. = FALSE)
  }

  rep_design <- svrep::as_bootstrap_design(design = data, replicates = bootstrap_reps)
  other_args <- list(...)

  replicate_results <- survey::withReplicates(
    design = rep_design,
    theta = function(pw, data_subset) {
      data_subset$..replicate_weights.. <- pw
      temp_design <- NULL
      temp_call <- try(getCall(data), silent = TRUE)
      if (!inherits(temp_call, "try-error") && !is.null(temp_call)) {
        temp_call$data <- quote(data_subset)
        temp_call$weights <- quote(~..replicate_weights..)
        temp_call$probs <- NULL
        temp_call$fpc <- NULL
        temp_design <- try(eval(temp_call), silent = TRUE)
        if (inherits(temp_design, "try-error")) temp_design <- NULL
      }
      if (is.null(temp_design)) {
        temp_design <- try(
          survey::svydesign(
            ids = data$call$ids,
            strata = data$call$strata,
            fpc = data$call$fpc,
            data = data_subset,
            weights = ~..replicate_weights..,
            nest = TRUE
          ),
          silent = TRUE
        )
        if (inherits(temp_design, "try-error")) temp_design <- NULL
      }
      if (is.null(temp_design)) {
        temp_design <- try(
          survey::svydesign(
            ids = ~1,
            data = data_subset,
            weights = ~..replicate_weights..
          ),
          silent = TRUE
        )
        if (inherits(temp_design, "try-error")) temp_design <- NULL
      }
      if (is.null(temp_design)) stop("Failed to reconstruct survey design for bootstrap replicates.")
      call_args <- c(list(data = temp_design, on_failure = "return"), other_args)
      fit <- suppressWarnings({
        do.call(estimator_func, call_args)
      })
      if (!is.null(fit$converged) && !fit$converged) {
        return(NA)
      }
      fit$y_hat
    },
    return.replicates = TRUE
  )

  replicate_estimates <- replicate_results$replicates
  failed_reps <- sum(is.na(replicate_estimates))
  if (failed_reps > 0 && (failed_reps / bootstrap_reps > 0.1)) {
    warning(sprintf(
      "%d out of %d bootstrap replicates failed to converge. The variance estimate may be unreliable.",
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
