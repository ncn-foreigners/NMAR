#' Construct EL Result Object
#' @keywords internal
new_nmar_result_el <- function(y_hat, se, weights, coefficients, vcov,
                               converged, diagnostics, data_info,
                               nmar_scaling_recipe, fitted_values, call) {
  diagnostics <- diagnostics %||% list()
  if (is.null(data_info$method)) data_info$method <- "Empirical Likelihood (EL)"
  outcome_name <- data_info$outcome_var %||% NA_character_
  trim_fraction <- diagnostics$trimmed_fraction %||% NA_real_

  sample <- list(
# Preserve the engine-supplied analysis-scale population total (N_pop).
# Use data_info$n_total when provided (set by dataframe/survey methods), not
# data_info$nobs. This ensures weights(object, scale = "population") sums
# to N_pop and survey-scale reporting is correct.
    n_total = data_info$n_total %||% data_info$nobs %||% NA_integer_,
    n_respondents = data_info$nobs_resp %||% NA_integer_,
    is_survey = isTRUE(data_info$is_survey),
    design = if (isTRUE(data_info$is_survey)) data_info$design else NULL
  )

  df_val <- NA_real_
  if (isTRUE(sample$is_survey) && !is.null(sample$design) && requireNamespace("survey", quietly = TRUE)) {
    df_val <- tryCatch(survey::degf(sample$design), error = function(e) NA_real_)
  }

  if (!is.list(coefficients)) {
    response_coeffs <- NULL
    nuisance_terms <- list()
  } else {
    response_coeffs <- coefficients$response_model %||% NULL
    nuisance_terms <- coefficients$nuisance %||% list()
  }

  inference <- list(
    variance_method = data_info$variance_method %||% NA_character_,
    df = df_val,
    message = diagnostics$vcov_message %||% NA_character_
  )

  meta <- list(
    engine_name = "empirical_likelihood",
    call = call,
    formula = data_info$formula
  )

  result <- new_nmar_result(
    estimate = y_hat,
    estimate_name = outcome_name,
    se = se,
    converged = converged,
    model = list(
      coefficients = response_coeffs,
      vcov = vcov,
      nuisance = nuisance_terms
    ),
    weights_info = list(values = weights, trimmed_fraction = trim_fraction),
    sample = sample,
    inference = inference,
    diagnostics = diagnostics,
    meta = meta,
    extra = list(
      nuisance = nuisance_terms,
      nmar_scaling_recipe = nmar_scaling_recipe,
      fitted_values = fitted_values
    ),
    class = "nmar_result_el"
  )

  result
}

#' Prepare inputs for EL estimation
#' @details Validates the two-sided outcome formula and constructs three
#'   internal formulas: outcome (~ outcome_var), response (for the missingness
#'   model), and auxiliary (RHS only, no intercept). Response-only predictors
#'   may include variables not on the outcome RHS; such variables enter only the
#'   response model (no auxiliary moment constraint). Only variables on the
#'   outcome RHS are treated as auxiliaries and, when provided, must match the
#'   names in `auxiliary_means`. See Qin, Leung and Shao (2002) for the EL
#'   formulation.
#' @keywords internal
#' Build internal EL formulas and data
#' @keywords internal
## Validation moved to src_dev/engines/el/impl/validation.R
