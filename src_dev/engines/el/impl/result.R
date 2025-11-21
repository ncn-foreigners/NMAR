#' Build EL result object (success or failure)
#' @keywords internal
el_build_result <- function(core_results, inputs, call, formula, engine_name = "empirical_likelihood") {
  diag_list <- core_results$diagnostics %||% list()

  if (!core_results$converged) {
    msg <- diag_list$message %||% NA_character_
    is_svy <- inherits(inputs$analysis_data, "survey.design")
    result <- new_nmar_result(
      estimate = NA_real_,
      estimate_name = inputs$outcome_expr,
      se = NA_real_,
      converged = FALSE,
      model = list(coefficients = NULL, vcov = NULL),
      weights_info = list(values = numeric(0), trimmed_fraction = NA_real_),
      sample = list(
        n_total = inputs$N_pop,
        n_respondents = sum(inputs$respondent_mask),
        is_survey = is_svy,
        design = if (isTRUE(is_svy)) inputs$analysis_data else NULL
      ),
      inference = list(variance_method = inputs$variance_method %||% NA_character_, df = NA_real_, message = msg),
      diagnostics = diag_list,
      meta = list(engine_name = engine_name, call = call, formula = formula),
      extra = list(nmar_scaling_recipe = core_results$nmar_scaling_recipe),
      class = "nmar_result_el"
    )
    return(validate_nmar_result(result, "nmar_result_el"))
  }

  result <- new_nmar_result_el(
    y_hat = core_results$y_hat,
    se = core_results$se,
    weights = core_results$weights,
    coefficients = core_results$coefficients,
    vcov = core_results$vcov,
    converged = TRUE,
    diagnostics = diag_list,
    inputs = inputs,
    nmar_scaling_recipe = core_results$nmar_scaling_recipe,
    fitted_values = core_results$fitted_values,
    call = call,
    formula = formula
  )
  validate_nmar_result(result, "nmar_result_el")
}
