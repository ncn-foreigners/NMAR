#' Induced-logit variance estimators
#'
#' @keywords internal
#' @noRd
il_compute_variance <- function(data,
                                formula,
                                point_estimate,
                                variance_method,
                                bootstrap_reps,
                                standardize,
                                control,
                                survey_design_policy,
                                resample_guard = NULL) {
  if (identical(variance_method, "none")) {
    return(list(
      se = NA_real_,
      message = "Variance skipped (variance_method='none')",
      diagnostics = list()
    ))
  }

  boot_runner <- function(data, formula, ...) {
    induced_logit(
      data = data,
      formula = formula,
      variance_method = "none",
      bootstrap_reps = 2,
      standardize = standardize,
      control = control,
      on_failure = "return",
      survey_design_policy = survey_design_policy,
      keep_fits = FALSE
    )
  }

  boot_args <- list(
    data = data,
    estimator_func = boot_runner,
    point_estimate = point_estimate,
    bootstrap_reps = bootstrap_reps,
    formula = formula
  )
  if (!is.null(resample_guard)) {
    boot_args$resample_guard <- resample_guard
  }

  boot_try <- tryCatch(
    list(result = do.call(bootstrap_variance, boot_args), message = "Calculation successful"),
    error = function(e) list(result = NULL, message = paste("Bootstrap failed:", e$message))
  )

  if (is.null(boot_try$result)) {
    return(list(se = NA_real_, message = boot_try$message, diagnostics = list()))
  }

  reps <- boot_try$result$replicates
  list(
    se = as.numeric(boot_try$result$se),
    message = boot_try$message,
    diagnostics = list(
      bootstrap_reps = as.integer(bootstrap_reps),
      bootstrap_failed_reps = as.integer(sum(!is.finite(reps))),
      bootstrap_success_reps = as.integer(sum(is.finite(reps)))
    )
  )
}
