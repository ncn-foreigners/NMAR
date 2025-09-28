#' Run method for EL engine
#' @keywords internal
#' @exportS3Method run_engine nmar_engine_el
run_engine.nmar_engine_el <- function(engine, spec) {
  response_predictors <- spec$response_predictors
  if (length(response_predictors) == 0) response_predictors <- NULL

  args <- list(
    data = spec$original_data,
    formula = spec$formula,
    response_predictors = response_predictors,
    auxiliary_means = engine$auxiliary_means,
    standardize = engine$standardize,
    trim_cap = engine$trim_cap,
    control = engine$control,
    on_failure = engine$on_failure,
    variance_method = engine$variance_method,
    variance_jacobian = engine$variance_jacobian,
    solver_jacobian = engine$solver_jacobian,
    variance_pseudoinverse = engine$variance_pseudoinverse,
    bootstrap_reps = engine$bootstrap_reps,
    suppress_warnings = engine$suppress_warnings,
    family = engine$family
  )

# Dispatch to EL implementation (data.frame or survey.design)
  res <- do.call(el, args)

# Ensure class includes the NMAR parent for downstream compatibility
  if (!inherits(res, "nmar_result_el")) {
    stop("EL engine did not return an 'nmar_result_el' object.")
  }
  res
}
