#' @exportS3Method NULL
run_engine.nmar_engine_exptilt_nonparam <- function(engine, task) {
# Nonparametric ET reuses the shared design prep to match the EL/ET data
# workflow (common scaling, auxiliary handling, and survey support)
  design_info <- prepare_nmar_design(
    task,
    standardize = FALSE,
    auxiliary_means = NULL,
    include_response = TRUE,
    include_auxiliary = TRUE
  )

# Use the normalized outcome list from the prepared design for consistency
  outcome_cols <- design_info$outcome
  common_covariates_from_formula <- design_info$auxiliary_vars
  instrumental_covariates_from_formula <- design_info$response_predictors

  if (length(common_covariates_from_formula) == 0) {
    stop("`covariates_outcome` in formula must specify at least one common covariate for the nonparametric EM method.")
  }

  if (length(instrumental_covariates_from_formula) == 0) {
    instrumental_covariates_used <- common_covariates_from_formula
  } else {
    instrumental_covariates_used <- instrumental_covariates_from_formula
  }

  model_results <- run_em_nmar_nonparametric(
    data = design_info$data,
    outcome_cols = outcome_cols,
    refusal_col = engine$refusal_col,
    common_covariates = common_covariates_from_formula,
    instrumental_covariates = instrumental_covariates_used,
    max_iter = engine$max_iter,
    tol = engine$tol_value
  )

  outcome_label <- paste(outcome_cols, collapse = " + ")
  n_total <- nrow(design_info$data)
  outcome_var <- outcome_cols[[1]]
  respondents <- if (length(outcome_cols)) {
    sum(rowSums(!is.na(design_info$data[outcome_cols])) > 0)
  } else {
    NA_integer_
  }

  sample_info <- list(
    n_total = n_total,
    n_respondents = respondents,
    is_survey = design_info$is_survey,
    design = design_info$survey_design
  )

  diagnostics <- list(
    iterations = model_results$iterations
  )

  meta <- list(
    engine_name = "exponential_tilting_nonparam",
    call = match.call(),
    formula = nmar_rebuild_partitioned_formula(
      base_formula = design_info$engine_formula,
      response_rhs_lang = design_info$response_rhs_lang,
      aux_rhs_lang = design_info$aux_rhs_lang,
      env = task$environment
    )
  )

  result <- new_nmar_result(
    estimate = NA_real_,
    estimate_name = outcome_label,
    se = NA_real_,
    converged = FALSE,
    model = list(coefficients = NULL, vcov = NULL),
    weights_info = list(values = numeric(0), trimmed_fraction = NA_real_),
    sample = sample_info,
    inference = list(
      variance_method = NA_character_,
      df = NA_real_,
      message = "Nonparametric ET diagnostic result (estimation not yet integrated into nmar_result schema)."
    ),
    diagnostics = diagnostics,
    meta = meta,
    extra = list(
      O_values = model_results$final_O_values,
      estimated_m_ij = model_results$final_m_ij,
      iterations = model_results$iterations,
      processed_data = model_results$final_data
    ),
    class = "nmar_result_exptilt_nonparam"
  )

  validate_nmar_result(result, "nmar_result_exptilt_nonparam")
}
