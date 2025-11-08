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

# Reconstruct a formula carrying response-only predictors to the right of `|`
  f_use <- nmar_rebuild_partitioned_formula(
    base_formula = design_info$engine_formula,
    response_rhs_lang = design_info$response_rhs_lang,
    aux_rhs_lang = design_info$aux_rhs_lang,
    env = task$environment
  )

  args <- list(
    data = design_info$survey_design %||% design_info$data,
    formula = f_use,
    user_formula = design_info$user_formula,
    auxiliary_means = design_info$auxiliary_means,
    standardize = design_info$standardize,
    design_matrices = design_info$design_matrices,
    outcome_label = design_info$outcome_label,
    n_total = engine$n_total,
    start = engine$start,
    trim_cap = engine$trim_cap,
    control = engine$control,
    on_failure = engine$on_failure,
    variance_method = engine$variance_method,
    bootstrap_reps = engine$bootstrap_reps,
    family = engine$family,
    trace_level = task$trace_level
  )

# Dispatch to EL implementation (data.frame or survey.design)
  res <- do.call(el, args)

# Ensure class includes the NMAR parent for downstream compatibility
  if (!inherits(res, "nmar_result_el")) {
    stop("EL engine did not return an 'nmar_result_el' object.")
  }
  res
}
