#' Induced-logit result constructors
#'
#' @keywords internal
#' @noRd
il_build_extra_payload <- function(fit, spec = NULL, keep_fits = FALSE) {
  out <- list(
    fitted_values = fit$models$fitted_values %||% numeric(0),
    scaling = fit$diagnostics$scaling %||% NULL
  )
  if (isTRUE(keep_fits)) {
    out$raw <- list(
      spec = spec,
      mu_fit = fit$models$mu_fit,
      induced_glm = fit$models$induced_glm
    )
  }
  out
}

#' @keywords internal
#' @noRd
il_build_success_result_iid <- function(spec,
                                        fit,
                                        variance_out,
                                        formula,
                                        call,
                                        variance_method,
                                        survey_design_policy,
                                        keep_fits) {
  new_nmar_result(
    estimate = as.numeric(fit$point$tau_hat),
    estimate_name = spec$outcome,
    se = variance_out$se,
    converged = TRUE,
    model = list(coefficients = fit$models$induced_glm_coef %||% NULL, vcov = fit$models$induced_glm_vcov %||% NULL),
    weights_info = list(values = NULL, trimmed_fraction = NA_real_),
    sample = list(
      n_total = nrow(spec$model_frame),
      n_respondents = fit$sample$n_respondents %||% sum(spec$respondent_mask),
      is_survey = FALSE,
      design = NULL
    ),
    inference = list(variance_method = variance_method, df = NA_real_, message = variance_out$message),
    diagnostics = il_build_engine_diagnostics(
      fit = fit,
      variance_diagnostics = variance_out$diagnostics,
      survey_design_policy = survey_design_policy
    ),
    meta = list(engine_name = "induced_logistic", call = call, formula = formula),
    extra = il_build_extra_payload(fit = fit, spec = spec, keep_fits = keep_fits),
    class = "nmar_result_induced_logit"
  )
}

#' @keywords internal
#' @noRd
il_build_success_result_survey <- function(spec,
                                           fit,
                                           design,
                                           backend,
                                           survey_assumptions,
                                           variance_out,
                                           formula,
                                           call,
                                           variance_method,
                                           survey_design_policy,
                                           keep_fits) {
  df_val <- il_survey_df(design)
  sum_w <- backend$sum_weights
  sum_w_resp <- sum(backend$weights_full[spec$respondent_mask])

  new_nmar_result(
    estimate = as.numeric(fit$point$tau_hat),
    estimate_name = spec$outcome,
    se = variance_out$se,
    converged = TRUE,
    model = list(coefficients = fit$models$induced_glm_coef %||% NULL, vcov = fit$models$induced_glm_vcov %||% NULL),
    weights_info = list(values = NULL, trimmed_fraction = NA_real_),
    sample = list(
      n_total = sum_w,
      n_respondents = fit$sample$n_respondents %||% sum(spec$respondent_mask),
      is_survey = TRUE,
      design = design
    ),
    inference = list(variance_method = variance_method, df = df_val, message = variance_out$message),
    diagnostics = il_build_engine_diagnostics(
      fit = fit,
      variance_diagnostics = variance_out$diagnostics,
      survey_design_policy = survey_design_policy,
      survey_assumptions = survey_assumptions,
      extras = list(
        sum_weights = as.numeric(sum_w),
        sum_respondent_weights = as.numeric(sum_w_resp)
      )
    ),
    meta = list(engine_name = "induced_logistic", call = call, formula = formula),
    extra = il_build_extra_payload(fit = fit, spec = spec, keep_fits = keep_fits),
    class = "nmar_result_induced_logit"
  )
}

#' @keywords internal
#' @noRd
il_build_failure_result <- function(message,
                                    formula,
                                    call,
                                    variance_method,
                                    survey_design_policy,
                                    n_total,
                                    is_survey,
                                    data_for_respondents,
                                    design = NULL,
                                    df = NA_real_) {
  outcome_name <- il_try_outcome_name(formula)
  n_resp <- il_try_n_respondents(data_for_respondents, outcome_name)

  new_nmar_result(
    estimate = NA_real_,
    estimate_name = outcome_name,
    se = NA_real_,
    converged = FALSE,
    model = list(coefficients = NULL, vcov = NULL),
    weights_info = list(values = NULL, trimmed_fraction = NA_real_),
    sample = list(
      n_total = n_total,
      n_respondents = n_resp,
      is_survey = is_survey,
      design = design
    ),
    inference = list(variance_method = variance_method, df = df, message = message),
    diagnostics = il_failure_diagnostics(message = message, survey_design_policy = survey_design_policy),
    meta = list(engine_name = "induced_logistic", call = call, formula = formula),
    extra = list(),
    class = "nmar_result_induced_logit"
  )
}
