#' Build EL result object (success or failure)
#' @keywords internal
el_build_result <- function(core_results, context, call, formula, engine_name = "empirical_likelihood") {
  meta <- list(
    outcome_var = context$outcome_var,
    n_total = context$N_pop,
    nobs_resp = length(context$respondent_indices),
    is_survey = context$is_survey,
    design = context$design,
    variance_method = context$variance_method
  )
  diag_list <- core_results$diagnostics %||% list()

  if (!core_results$converged) {
    msg <- diag_list$message %||% NA_character_
    result <- new_nmar_result(
      estimate = NA_real_,
      estimate_name = meta$outcome_var,
      se = NA_real_,
      converged = FALSE,
      model = list(coefficients = NULL, vcov = NULL),
      weights_info = list(values = numeric(0), trimmed_fraction = NA_real_),
      sample = list(
        n_total = meta$n_total,
        n_respondents = meta$nobs_resp,
        is_survey = meta$is_survey,
        design = meta$design
      ),
      inference = list(variance_method = meta$variance_method, df = NA_real_, message = msg),
      diagnostics = diag_list,
      meta = list(engine_name = engine_name, call = call, formula = formula),
      extra = list(nmar_scaling_recipe = core_results$nmar_scaling_recipe),
      class = "nmar_result_el"
    )
    return(validate_nmar_result(result, "nmar_result_el"))
  }

  meta$call <- call
  meta$formula <- formula

  new_nmar_result_el(
    y_hat = core_results$y_hat,
    se = core_results$se,
    weights = core_results$weights,
    coefficients = core_results$coefficients,
    vcov = core_results$vcov,
    converged = TRUE,
    diagnostics = diag_list,
    data_info = meta,
    nmar_scaling_recipe = core_results$nmar_scaling_recipe,
    fitted_values = core_results$fitted_values,
    call = call,
    formula = formula
  ) |>
    validate_nmar_result("nmar_result_el")
}
