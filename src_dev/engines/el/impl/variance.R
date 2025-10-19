#' EL variance components and sandwich assembly
#'
#' @details
#' Assembles the respondent score matrix, computes the covariance of total scores
#' (IID or design-based), and constructs the delta variance for \eqn{\widehat{Y}}
#' via \eqn{\nabla g^\top A^{-1} B A^{-\top} \nabla g}. The Jacobian \eqn{A} and
#' the score covariance \eqn{B} are evaluated at the empirical-likelihood (EL)
#' solution of the full stacked system. For IID data, \eqn{B} is computed from
#' centered respondent score totals (asymptotically equivalent at the root); for
#' survey data, \eqn{B} uses the design-based variance of totals. This follows
#' Qin, Leung and Shao (2002). When trimming is active, the mean functional becomes
#' non-smooth; the delta method then uses a numeric gradient and may return \code{NA},
#' in which case bootstrap variance is recommended.
#'
#' @name el_variance
#' @keywords internal
NULL

#' @keywords internal
el_compute_score_contrib <- function(family,
                                     response_model_matrix_scaled,
                                     auxiliary_matrix_scaled,
                                     mu_x_scaled,
                                     eta_i_hat,
                                     w_i_hat,
                                     W_hat,
                                     denominator_hat,
                                     lambda_W_hat,
                                     respondent_weights,
                                     full_data) {
  K_aux <- if (is.null(auxiliary_matrix_scaled) || ncol(auxiliary_matrix_scaled) == 0) 0 else ncol(auxiliary_matrix_scaled)

# Score wrt eta for log-likelihood: d/deta log p(eta) = mu.eta / p
  p_hat <- family$linkinv(eta_i_hat)
  p_hat <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
  denom_guard <- pmax(as.numeric(denominator_hat), 1e-8)
# For IID data.frames, B is built from totals over respondents, so include respondent_weights
# For survey designs, totals are computed by the survey design (svytotal), so do not multiply here
  use_resp_w <- !inherits(full_data, "survey.design")
  w_fac <- if (use_resp_w) as.vector(respondent_weights) else 1
  U_beta <- (response_model_matrix_scaled * as.vector(
    (family$mu.eta(eta_i_hat) / p_hat) - lambda_W_hat * family$mu.eta(eta_i_hat) / denom_guard
  )) * w_fac
  U_W <- ((w_i_hat - W_hat) / denom_guard) * w_fac
  U_aux <- if (K_aux > 0) (sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-") * as.vector(1 / denom_guard)) * w_fac else matrix(nrow = nrow(response_model_matrix_scaled), ncol = 0)
  out <- cbind(U_beta, score_W = U_W, U_aux)
# Ensure stable, aligned column names for B with A (beta, z, lambda order)
  cn_beta <- colnames(response_model_matrix_scaled)
  cn <- c(cn_beta, "(W) (logit)", if (K_aux > 0) paste0("lambda_", colnames(auxiliary_matrix_scaled)) else NULL)
  if (!is.null(cn) && ncol(out) == length(cn)) colnames(out) <- cn
  out
}

#' @keywords internal
el_compute_B_matrix <- function(U_matrix_resp, full_data, compute_score_variance_func) {
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
                                      estimates) {
  U_matrix_resp <- el_compute_score_contrib(
    family = family,
    response_model_matrix_scaled = response_model_matrix_scaled,
    auxiliary_matrix_scaled = auxiliary_matrix_scaled,
    mu_x_scaled = mu_x_scaled,
    eta_i_hat = eta_i_hat,
    w_i_hat = w_i_hat,
    W_hat = W_hat,
    denominator_hat = denominator_hat,
    lambda_W_hat = lambda_W_hat,
    respondent_weights = respondent_weights,
    full_data = full_data
  )
  B_matrix <- el_compute_B_matrix(U_matrix_resp, full_data, compute_score_variance_func)
# Diagnostics: record B min eigenvalue for PSD check
  B_min_eig <- tryCatch({
    ev <- eigen(0.5 * (B_matrix + t(B_matrix)), symmetric = TRUE, only.values = TRUE)$values
    suppressWarnings(min(ev))
  }, error = function(e) NA_real_)
# Assemble sandwich via two linear solves
  vcov_matrix_sandwich_scaled <- tryCatch({
# Form Sigma = A^{-1} B A^{-T} via two solves without explicit inverse
    X <- solve(A_matrix, B_matrix) # X = A^{-1} B
    Sig <- X %*% solve(t(A_matrix)) # Sigma = A^{-1} B A^{-T}
# Symmetrize lightly to reduce roundoff asymmetry
    0.5 * (Sig + t(Sig))
  }, error = function(e) NULL)
# Preserve parameter names on Sigma to enable safe alignment of gradients
  if (!is.null(colnames(A_matrix))) {
    dimnames(vcov_matrix_sandwich_scaled) <- list(colnames(A_matrix), colnames(A_matrix))
  }
# Gradient of mean functional: prefer analytic when no trimming,
# but fall back to numeric if the analytic gradient is suspiciously small
  grad_g <- NULL
  analytic_used <- FALSE
  if (is.infinite(trim_cap)) {
    X_beta <- response_model_matrix_scaled
    K_aux <- if (is.null(auxiliary_matrix_scaled) || ncol(auxiliary_matrix_scaled) == 0) 0 else ncol(auxiliary_matrix_scaled)
    Xc <- if (K_aux > 0) sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-") else matrix(nrow = nrow(X_beta), ncol = 0)
    grad_g_try <- el_grad_g_analytic(
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
    if (is.numeric(grad_g_try) && all(is.finite(grad_g_try)) && sum(abs(grad_g_try)) > 1e-8) {
      grad_g <- grad_g_try
      analytic_used <- TRUE
    }
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
  var_y_hat <- if (is.null(vcov_matrix_sandwich_scaled)) NA_real_ else as.numeric(t(grad_g) %*% vcov_matrix_sandwich_scaled %*% grad_g)
  list(
    var_y_hat = var_y_hat,
    vcov_matrix_sandwich_scaled = vcov_matrix_sandwich_scaled,
    B = B_matrix,
    B_min_eig = B_min_eig
  )
}
