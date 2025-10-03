#' Run method for EL engine
#' @keywords internal
#' @exportS3Method run_engine nmar_engine_el
run_engine.nmar_engine_el <- function(engine, task) {
  # Reuse the shared design preparation so EL mirrors the ET workflow for
  # survey designs, scaling, and auxiliary moment injection
  design_info <- prepare_nmar_design(
    task,
    standardize = engine$standardize,
    auxiliary_means = engine$auxiliary_means,
    include_response = TRUE,
    include_auxiliary = TRUE
  )

  response_predictors <- design_info$response_predictors
  if (length(response_predictors) == 0) response_predictors <- NULL

  args <- list(
    data = design_info$survey_design %||% design_info$data,
    formula = task$formula,
    response_predictors = response_predictors,
    auxiliary_means = design_info$auxiliary_means,
    standardize = design_info$standardize,
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
