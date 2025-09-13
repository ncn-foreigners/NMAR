#' @exportS3Method NULL
run_engine.nmar_engine_exptilt_nonparam <- function(engine, formula, data,response_predictors) {

  # common_covariates_from_formula <- all.vars(formula$covariates_outcome)
  # instrumental_covariates_from_formula <- all.vars(formula$covariates_missingness)
  outcome_cols <- all.vars(formula[[2]])
  common_covariates_from_formula <- all.vars(formula[[3]])
  instrumental_covariates_from_formula <- response_predictors


  if (length(common_covariates_from_formula) == 0) {
    stop("`covariates_outcome` in formula must specify at least one common covariate for the nonparametric EM method.")
  }


  if (length(instrumental_covariates_from_formula) == 0) {
    instrumental_covariates_used <- common_covariates_from_formula
  } else {
    instrumental_covariates_used <- instrumental_covariates_from_formula
  }


  model_results <- run_em_nmar_nonparametric(
    data = data,
    # outcome_cols = engine$outcome_cols,
    outcome_cols = outcome_cols,
    refusal_col = engine$refusal_col,
    common_covariates = common_covariates_from_formula,
    instrumental_covariates = instrumental_covariates_used,
    max_iter = engine$max_iter,
    tol = engine$tol_value
  )


  return(structure(list(
    O_values = model_results$final_O_values,
    estimated_m_ij = model_results$final_m_ij,
    iterations = model_results$iterations,
    processed_data = model_results$final_data
  ), class = "nmar_result"))
}
