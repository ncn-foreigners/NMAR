#' Empirical likelihood estimating equations
#' @details Returns a function that evaluates the full stacked EL system for
#'   \eqn{\theta = (\beta, z, \lambda_x)} with \eqn{z = \operatorname{logit}(W)};
#'   the auxiliary block is omitted when no constraints are present. The
#'   response-model score with respect to the linear predictor uses the
#'   derivative of the Bernoulli log-likelihood, which is valid for both logit
#'   and probit links. This matches the semiparametric EL system in Qin, Leung
#'   and Shao (2002). Denominator guards are applied to avoid invalid weights.
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193-200.
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
    mu_eta_i <- family$mu.eta(eta_i)
# For logit, d/deta log w equals 1 - w; for probit, use family score implementation
    dlw_i <- family$score_eta(eta_i, delta = 1)
    denominator <- 1 + lambda_W * (w_i - W_bounded)
    if (K_aux > 0) denominator <- denominator + as.vector(X_centered %*% lambda_x)
    inv_denominator <- 1 / pmax(denominator, nmar_get_el_denom_floor())
    scalar_beta_term <- dlw_i - lambda_W * mu_eta_i * inv_denominator
    eq_betas <- shared_weighted_Xty(response_model_matrix, respondent_weights, scalar_beta_term)
    eq_W <- sum(respondent_weights * (w_i - W_bounded) * inv_denominator)
    eq_constraints <- if (K_aux > 0) shared_weighted_Xty(X_centered, respondent_weights, inv_denominator) else numeric(0)
    c(as.vector(eq_betas), eq_W, as.vector(eq_constraints))
  }
}
