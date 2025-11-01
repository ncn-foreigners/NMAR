#' Analytical Jacobian for empirical likelihood
#' @details Builds the block Jacobian \eqn{A=\partial F/\partial \theta} for the
#'   full EL system with \eqn{\theta=(\beta, z, \lambda_x)} and \eqn{z=\operatorname{logit}(W)}.
#'   The response-model score derivative with respect to the linear predictor
#'   is used for both logit and probit links. This is consistent with the EL
#'   formulation in Qin, Leung and Shao (2002). Denominator guards are applied
#'   when forming terms that depend on the empirical likelihood weights.
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193-200.
#' @keywords internal
build_el_jacobian <- function(family, response_model_matrix, auxiliary_matrix,
                              respondent_weights, N_pop, n_resp_weighted, mu_x_scaled) {
  force(family)
  force(response_model_matrix)
  force(auxiliary_matrix)
  force(respondent_weights)
  force(N_pop)
  force(n_resp_weighted)
  force(mu_x_scaled)
  if (is.null(family) || is.null(family$d2mu.deta2)) {
    return(NULL)
  }
  if (is.null(response_model_matrix) || !is.matrix(response_model_matrix)) stop("response_model_matrix must be a matrix.")
  n_resp <- nrow(response_model_matrix)
  K_beta <- ncol(response_model_matrix)
  K_aux <- if (is.null(auxiliary_matrix) || ncol(auxiliary_matrix) == 0) 0 else ncol(auxiliary_matrix)
  if (length(respondent_weights) != n_resp) stop("Length of respondent_weights must equal nrow(response_model_matrix).")
  if (K_aux == 0) {
    auxiliary_matrix_mat <- matrix(nrow = n_resp, ncol = 0)
    mu_x_scaled_vec <- numeric(0)
  } else {
    auxiliary_matrix_mat <- as.matrix(auxiliary_matrix)
    if (nrow(auxiliary_matrix_mat) != n_resp) stop("auxiliary_matrix must match rows.")
    if (is.null(names(mu_x_scaled))) stop("mu_x_scaled must be named.")
    if (!all(colnames(auxiliary_matrix_mat) %in% names(mu_x_scaled))) stop("Names of mu_x_scaled must include all columns of auxiliary_matrix.")
    mu_x_scaled_vec <- as.numeric(mu_x_scaled[colnames(auxiliary_matrix_mat)])
    names(mu_x_scaled_vec) <- colnames(auxiliary_matrix_mat)
  }
  C_const <- (N_pop / n_resp_weighted) - 1
  function(params) {
    if (length(params) != (K_beta + 1 + K_aux)) stop("Parameter vector length mismatch.")
    beta_vec <- as.numeric(params[1:K_beta])
    z <- as.numeric(params[K_beta + 1])
    W <- plogis(z)
    W <- min(max(W, 1e-12), 1 - 1e-12)
    dWb_dTheta <- if (W > 1e-12 && W < 1 - 1e-12) W * (1 - W) else 0
    lambda_x <- if (K_aux > 0) as.numeric(params[(K_beta + 2):length(params)]) else numeric(0)
    W_bounded <- W
    lambda_W <- C_const / (1 - W_bounded)
    d_lambda_W_dWb <- C_const / (1 - W_bounded)^2
    d_lambda_W_dTheta <- d_lambda_W_dWb * dWb_dTheta
    ETA_CAP <- get_eta_cap()
    eta_raw <- as.vector(response_model_matrix %*% beta_vec)
    eta_i <- pmax(pmin(eta_raw, ETA_CAP), -ETA_CAP)
    p_i <- family$linkinv(eta_i)
    m_i <- family$mu.eta(eta_i)
    m2_i <- family$d2mu.deta2(eta_i)
    X_centered <- if (K_aux > 0) sweep(auxiliary_matrix_mat, 2, mu_x_scaled_vec, "-") else matrix(nrow = n_resp, ncol = 0)
# QLS Eq. (5): Di = 1 + lambda_W * (w_i - W) + (Xc %*% lambda_x)
    denominator <- 1 + lambda_W * (p_i - W_bounded)
    if (K_aux > 0) denominator <- denominator + as.vector(X_centered %*% lambda_x)
    denom_floor <- nmar_get_el_denom_floor()
# Active mask for max(Di, floor) kink: derivative is zero when clamped
    active <- as.numeric(denominator > denom_floor)
    denominator <- pmax(denominator, denom_floor)
    inv_denom <- 1 / denominator
    inv_denom_sq <- inv_denom^2
    dden_dTheta <- active * (d_lambda_W_dTheta * (p_i - W_bounded) - lambda_W * dWb_dTheta)
    dden_deta <- active * (lambda_W * m_i)
# Score wrt eta for log-likelihood: d/deta log p(eta) = m_i / p_i
    p_i_clipped <- pmin(pmax(p_i, 1e-12), 1 - 1e-12)
    dlw_i <- m_i / p_i_clipped
# beta block term (QLS Eq. 9)
    beta_eq_term <- dlw_i - lambda_W * m_i * inv_denom
# Derivative wrt eta: d/deta(m/p) = (m2 * p - m^2) / p^2
    d_deta_logw <- (m2_i * p_i_clipped - m_i^2) / (p_i_clipped^2)
    d_betaeq_deta <- d_deta_logw - lambda_W * m2_i * inv_denom + (lambda_W^2) * (m_i^2) * inv_denom_sq
    d_betaeq_dTheta <- -d_lambda_W_dTheta * m_i * inv_denom + lambda_W * m_i * inv_denom_sq * dden_dTheta
    d_betaeq_dlambda_mat <- if (K_aux > 0) lambda_W * m_i * inv_denom_sq * X_centered else matrix(nrow = n_resp, ncol = 0)
    w_eff_11 <- as.numeric(respondent_weights * d_betaeq_deta)
    J11 <- shared_weighted_gram(response_model_matrix, w_eff_11)
    J12 <- shared_weighted_Xty(response_model_matrix, respondent_weights, d_betaeq_dTheta)
    J13 <- if (K_aux > 0) shared_weighted_XtY(response_model_matrix, respondent_weights, as.matrix(d_betaeq_dlambda_mat)) else matrix(nrow = K_beta, ncol = 0)
    term21 <- m_i * inv_denom - (p_i - W_bounded) * inv_denom_sq * (lambda_W * m_i) * active
    J21 <- t(shared_weighted_Xty(response_model_matrix, respondent_weights, term21))
    term22 <- -dWb_dTheta * inv_denom - (p_i - W_bounded) * inv_denom_sq * dden_dTheta
    J22 <- sum(as.numeric(respondent_weights * term22))
    J23 <- if (K_aux > 0) t(shared_weighted_Xty(X_centered, respondent_weights, (-(p_i - W_bounded) * inv_denom_sq * active))) else matrix(nrow = 1, ncol = 0)
    if (K_aux > 0) {
      term31 <- -dden_deta * inv_denom_sq
      J31 <- shared_weighted_XtY(X_centered, as.numeric(respondent_weights * term31), response_model_matrix)
      term32 <- -dden_dTheta * inv_denom_sq
      J32 <- shared_weighted_Xty(X_centered, respondent_weights, term32)
# J33 = - Xc' diag(inv_denom^2) Xc scaled by respondent weights; SPD path (gated by active)
      J33 <- -shared_weighted_gram(X_centered, as.numeric(respondent_weights * (inv_denom_sq * active)))
    } else {
      J31 <- matrix(nrow = 0, ncol = K_beta)
      J32 <- matrix(nrow = 0, ncol = 1)
      J33 <- matrix(nrow = 0, ncol = 0)
    }
    top <- cbind(J11, J12, J13)
    middle <- cbind(J21, matrix(J22, nrow = 1, ncol = 1), J23)
    bottom <- cbind(J31, J32, J33)
    full_mat <- rbind(top, middle, bottom)
    param_names <- c(colnames(response_model_matrix), "(W) (logit)", if (K_aux > 0) paste0("lambda_", colnames(X_centered)) else NULL)
    if (!is.null(param_names) && length(param_names) == ncol(full_mat)) {
      colnames(full_mat) <- rownames(full_mat) <- param_names
    }
    as.matrix(full_mat)
  }
}
