#' EL variance components and sandwich assembly
#' @details Assembles the respondent score matrix, computes the covariance of
#'   total scores (IID crossproduct or survey design-based), and constructs the
#'   delta variance for \eqn{\hat Y} via \eqn{\nabla g\,A^{-1} B A^{-T}\,\nabla g^T}.
#'   Consistent with the asymptotic results in Qin, Leung and Shao (2002).
#' @keywords internal
el_compute_score_contrib <- function(family,
                                     response_model_matrix_scaled,
                                     auxiliary_matrix_scaled,
                                     mu_x_scaled,
                                     eta_i_hat,
                                     w_i_hat,
                                     W_hat,
                                     denominator_hat,
                                     lambda_W_hat) {
  K_aux <- if (is.null(auxiliary_matrix_scaled) || ncol(auxiliary_matrix_scaled) == 0) 0 else ncol(auxiliary_matrix_scaled)

# Score wrt eta for log-likelihood: d/deta log p(eta) = mu.eta / p
  p_hat <- family$linkinv(eta_i_hat)
  p_hat <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
  denom_guard <- pmax(as.numeric(denominator_hat), 1e-8)
  U_beta <- response_model_matrix_scaled * as.vector(
    (family$mu.eta(eta_i_hat) / p_hat) - lambda_W_hat * family$mu.eta(eta_i_hat) / denom_guard
  )
  U_W <- (w_i_hat - W_hat) / denom_guard
  U_aux <- if (K_aux > 0) sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-") * as.vector(1 / denom_guard) else matrix(nrow = nrow(response_model_matrix_scaled), ncol = 0)
  cbind(U_beta, score_W = U_W, U_aux)
}

#' @keywords internal
compute_B_matrix <- function(U_matrix_resp, full_data, compute_score_variance_func) {
  compute_score_variance_func(U_matrix_resp, full_data)
}

#' @keywords internal
el_build_mean_fn <- function(family,
                             response_model_matrix_scaled,
                             auxiliary_matrix_scaled,
                             mu_x_scaled,
                             respondent_weights,
                             N_pop,
                             n_resp_weighted,
                             trim_cap,
                             outcome_vec) {
  K_beta <- ncol(response_model_matrix_scaled)
  K_aux <- if (is.null(auxiliary_matrix_scaled) || ncol(auxiliary_matrix_scaled) == 0) 0 else ncol(auxiliary_matrix_scaled)
  function(p) {
    eta <- response_model_matrix_scaled %*% p[1:K_beta]
    cap <- get_eta_cap()
    eta <- pmax(pmin(eta, cap), -cap)
    w <- family$linkinv(eta)
    W_raw <- plogis(p[K_beta + 1])
    W_b <- max(min(W_raw, 1 - 1e-12), 1e-12)
    lambda_W <- ((N_pop / n_resp_weighted) - 1) / (1 - W_b)
    denom <- 1 + lambda_W * (w - W_b)
    if (K_aux > 0) denom <- denom + as.vector(sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-") %*% p[(K_beta + 2):length(p)])
    denom <- pmax(as.numeric(denom), 1e-8)
    pi_final <- trim_weights(respondent_weights / denom, cap = trim_cap)$weights
    s_pi <- sum(pi_final)
    if (!is.finite(s_pi) || s_pi <= 0) return(NA_real_)
    sum(pi_final * outcome_vec) / s_pi
  }
}

#' @keywords internal
el_compute_delta_variance <- function(A_matrix,
                                      family,
                                      response_model_matrix_scaled,
                                      auxiliary_matrix_scaled,
                                      mu_x_scaled,
                                      eta_i_hat,
                                      w_i_hat,
                                      W_hat,
                                      denominator_hat,
                                      lambda_W_hat,
                                      full_data,
                                      compute_score_variance_func,
                                      respondent_weights,
                                      N_pop,
                                      n_resp_weighted,
                                      trim_cap,
                                      outcome_vec,
                                      estimates,
                                      variance_ridge = FALSE,
                                      variance_pseudoinverse = FALSE) {
  inv_res <- invert_jacobian(A_matrix,
    variance_ridge = variance_ridge,
    variance_pseudoinverse = variance_pseudoinverse
  )
  A_inv <- inv_res$inv
  used_pseudoinverse <- isTRUE(inv_res$used_pinv)
  used_ridge <- isTRUE(inv_res$used_ridge)
  invert_rule <- inv_res$invert_rule
  U_matrix_resp <- el_compute_score_contrib(
    family = family,
    response_model_matrix_scaled = response_model_matrix_scaled,
    auxiliary_matrix_scaled = auxiliary_matrix_scaled,
    mu_x_scaled = mu_x_scaled,
    eta_i_hat = eta_i_hat,
    w_i_hat = w_i_hat,
    W_hat = W_hat,
    denominator_hat = denominator_hat,
    lambda_W_hat = lambda_W_hat
  )
  B_matrix <- compute_B_matrix(U_matrix_resp, full_data, compute_score_variance_func)
  vcov_matrix_sandwich_scaled <- (A_inv %*% B_matrix %*% t(A_inv))
# Prefer analytic gradient when no trimming (smooth ratio)
  grad_g <- NULL
  if (is.infinite(trim_cap)) {
    X_beta <- response_model_matrix_scaled
    K_aux <- if (is.null(auxiliary_matrix_scaled) || ncol(auxiliary_matrix_scaled) == 0) 0 else ncol(auxiliary_matrix_scaled)
    Xc <- if (K_aux > 0) sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-") else matrix(nrow = nrow(X_beta), ncol = 0)
    grad_g <- el_grad_g_analytic(
      family = family,
      X_beta = X_beta,
      Xc = Xc,
      mu_x_scaled = mu_x_scaled,
      respondent_weights = respondent_weights,
      eta_i_hat = eta_i_hat,
      w_i_hat = w_i_hat,
      W_hat = W_hat,
      denominator_hat = denominator_hat,
      lambda_W_hat = lambda_W_hat,
      outcome_vec = outcome_vec,
      n_resp_weighted = n_resp_weighted,
      N_pop = N_pop
    )
  }
  if (is.null(grad_g)) {
    g_fn <- el_build_mean_fn(
      family = family,
      response_model_matrix_scaled = response_model_matrix_scaled,
      auxiliary_matrix_scaled = auxiliary_matrix_scaled,
      mu_x_scaled = mu_x_scaled,
      respondent_weights = respondent_weights,
      N_pop = N_pop,
      n_resp_weighted = n_resp_weighted,
      trim_cap = trim_cap,
      outcome_vec = outcome_vec
    )
    grad_g <- grad_numeric(estimates, g_fn)
  }
  var_y_hat <- as.numeric(t(grad_g) %*% vcov_matrix_sandwich_scaled %*% grad_g)
  list(
    var_y_hat = var_y_hat,
    vcov_matrix_sandwich_scaled = vcov_matrix_sandwich_scaled,
    used_pseudoinverse = used_pseudoinverse,
    used_ridge = used_ridge,
    invert_rule = invert_rule
  )
}
