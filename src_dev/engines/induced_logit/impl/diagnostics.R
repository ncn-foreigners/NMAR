#' Induced-logit diagnostics builders
#'
#' @keywords internal
#' @noRd
il_compact_warnings <- function(warnings) {
  list(
    mu = unique(warnings$mu %||% character()),
    induced_glm = unique(warnings$induced_glm %||% character())
  )
}

#' @keywords internal
#' @noRd
il_build_engine_diagnostics <- function(fit,
                                        variance_diagnostics = list(),
                                        survey_design_policy = "strict",
                                        survey_assumptions = NULL,
                                        extras = list()) {
  pt <- fit$point %||% list()
  rm <- fit$diagnostics$response_model %||% list()
  warn <- fit$warnings %||% list()
  out <- c(
    list(
      eta_hat = as.numeric(pt$eta_hat),
      gamma_hat_paper = as.numeric(pt$gamma_hat_paper),
      mu_bar = as.numeric(pt$mu_bar),
      m2_over_m1 = as.numeric(pt$m2_over_m1),
      log_m1_hat = as.numeric(pt$log_m1_hat),
      alpha_hat_paper = as.numeric(pt$alpha_hat_paper),
      alpha0_hat_paper = as.numeric(pt$alpha0_hat_paper),
      response_model_rank = as.numeric(rm$rank %||% NA_real_),
      response_model_ncol = as.numeric(rm$ncol %||% NA_real_),
      response_model_condition_number = as.numeric(rm$condition_number %||% NA_real_),
      variance = variance_diagnostics %||% list(),
      survey_design_policy = survey_design_policy,
      warnings = il_compact_warnings(warn)
    ),
    extras
  )
  if (!is.null(survey_assumptions)) {
    out$survey_assumptions <- survey_assumptions
  }
  out
}

#' @keywords internal
#' @noRd
il_failure_diagnostics <- function(message, survey_design_policy = "strict") {
  list(
    message = message,
    survey_design_policy = survey_design_policy
  )
}
