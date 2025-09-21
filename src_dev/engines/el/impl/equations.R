#' Empirical likelihood estimating equations
#' @details Returns a function evaluating the stacked EL equations for
#'   \eqn{\beta}, \eqn{W} and (optionally) auxiliary multipliers. The
#'   response‑model score uses the derivative of the Bernoulli log‑likelihood
#'   with respect to the linear predictor, \eqn{\partial\log p(\eta)/\partial\eta
#'   = \mu_\eta(\eta)/p(\eta)}, which is valid for both logit and probit links.
#'   This matches the semiparametric EL system in Qin, Leung and Shao (2002).
#' @keywords internal
el_build_equation_system <- function(family, response_model_matrix, auxiliary_matrix,
                                  respondent_weights, N_pop, n_resp_weighted, mu_x_scaled) {
  force(family)
  force(response_model_matrix)
  force(auxiliary_matrix)
  force(respondent_weights)
  force(N_pop)
  force(n_resp_weighted)
  force(mu_x_scaled)
  K_beta <- ncol(response_model_matrix)
  K_aux <- if (is.null(auxiliary_matrix) || ncol(auxiliary_matrix) == 0) 0 else ncol(auxiliary_matrix)
  X_centered <- if (K_aux > 0) sweep(auxiliary_matrix, 2, mu_x_scaled, "-") else matrix(nrow = nrow(response_model_matrix), ncol = 0)
  function(params) {
    beta_vec <- params[1:K_beta]
    z <- params[K_beta + 1]
    W <- plogis(z)
    W <- min(max(W, 1e-12), 1 - 1e-12)
    lambda_x <- if (K_aux > 0) params[(K_beta + 2):length(params)] else numeric(0)
    W_bounded <- W
    lambda_W <- ((N_pop / n_resp_weighted) - 1) / (1 - W_bounded)
    ETA_CAP <- get_eta_cap()
    eta_raw <- as.vector(response_model_matrix %*% beta_vec)
    eta_i <- pmax(pmin(eta_raw, ETA_CAP), -ETA_CAP)
    w_i <- family$linkinv(eta_i)
    score_i <- family$score_eta(eta_i, delta = 1)
    mu_eta_i <- family$mu.eta(eta_i)
    denominator <- 1 + lambda_W * (w_i - W_bounded)
    if (K_aux > 0) denominator <- denominator + as.vector(X_centered %*% lambda_x)
    inv_denominator <- 1 / pmax(denominator, 1e-8)
    scalar_beta_term <- score_i - lambda_W * mu_eta_i * inv_denominator
    eq_betas <- t(response_model_matrix) %*% (respondent_weights * scalar_beta_term)
    eq_W <- sum(respondent_weights * (w_i - W_bounded) * inv_denominator)
    eq_constraints <- if (K_aux > 0) t(X_centered) %*% (respondent_weights * inv_denominator) else numeric(0)
    c(as.vector(eq_betas), eq_W, as.vector(eq_constraints))
  }
}
