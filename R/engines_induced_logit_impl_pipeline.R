#' Induced-logit fit pipeline
#'
#' @keywords internal
#' @noRd
il_fit_from_backend <- function(spec, backend, standardize, control = list()) {
  mu_res <- backend$fit_mu(spec, standardize = standardize)
  mu_hat <- mu_res$mu_hat

  glm_control <- il_glm_control(control)
  resp_res <- backend$fit_resp(spec, mu_hat = mu_hat, standardize = standardize, glm_control = glm_control)
  beta_glm <- resp_res$beta_glm
  glm_intercept <- beta_glm[["(Intercept)"]]
  if (is.null(glm_intercept) || !is.finite(glm_intercept)) {
    stop("Internal error: induced-logit model is missing a finite intercept coefficient.", call. = FALSE)
  }

  tau_res <- il_compute_estimand(
    spec = spec,
    y_vec = backend$vars[[spec$outcome]],
    mu_hat = mu_hat,
    gamma_hat_paper = resp_res$gamma_hat_paper,
    glm_intercept = glm_intercept,
    weights_full = backend$weights_full
  )

  scaling_out <- NULL
  if (isTRUE(standardize)) {
    scaling_out <- list(mu = mu_res$mu_recipe %||% NULL, x1 = resp_res$x1_recipe %||% NULL)
  }

  list(
    point = c(
      list(gamma_hat_paper = as.numeric(resp_res$gamma_hat_paper)),
      tau_res
    ),
    models = list(
      mu_fit = mu_res$mu_fit,
      induced_glm = resp_res$induced_glm,
      induced_glm_coef = beta_glm,
      induced_glm_vcov = resp_res$vcov,
      fitted_values = resp_res$fitted
    ),
    warnings = list(
      mu = mu_res$warnings %||% character(),
      induced_glm = resp_res$warnings %||% character()
    ),
    diagnostics = list(
      response_model = resp_res$identifiability %||% list(),
      scaling = scaling_out
    ),
    sample = list(
      n_total = nrow(backend$vars),
      n_respondents = sum(spec$respondent_mask)
    )
  )
}
