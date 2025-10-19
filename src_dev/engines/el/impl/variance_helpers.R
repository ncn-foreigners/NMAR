#' Variance helper utilities

#' Delta variance computation helpers
#'
#' Assembles the analytic sandwich variance at the EL solution using numerically
#' robust two-solve identities, and returns the mean standard error and the
#' response-model coefficient covariance. Specifically:
#' - Mean: var(g) = x' B x with t(A) x = ∇g, which equals g' (A^{-1} B A^{-T}) g
#'   but avoids explicit inversion and preserves non-negativity when B is PSD.
#' - Coefficients: V_beta = X_beta' B X_beta with t(A) X_beta = E_beta, which is
#'   the beta-block of A^{-1} B A^{-T} assembled without forming the full matrix.
#' B is the covariance of respondent score totals (IID) or the design-based
#' variance of totals (survey). A is the analytic Jacobian of the stacked EL
#' system. In trimming or fragile regimes, the delta method returns NA with a
#' clear warning; bootstrap is recommended in those cases. See Qin, Leung and
#' Shao (2002).
#'
#' @return List with entries: `se_y_hat`, `vcov_unscaled`, `vcov_result`, and
#'   `diag` (a list of diagnostics: A_condition and grad_source).
#' @keywords internal
el_variance_delta <- function(equation_system_func,
                              analytical_jac_func,
                              estimates,
                              family,
                              response_model_matrix_scaled,
                              response_model_matrix_unscaled,
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
                              K_beta,
                              standardize,
                              nmar_scaling_recipe) {
  if (is.finite(trim_cap)) {
    return(list(
      se_y_hat = NA_real_,
      vcov_unscaled = matrix(NA_real_, K_beta, K_beta, dimnames = list(colnames(response_model_matrix_unscaled), colnames(response_model_matrix_unscaled))),
      vcov_result = NULL,
      vcov_message = "Delta variance not applicable with trimming; use bootstrap.",
      diag = list(
        A_condition = NA_real_
      )
    ))
  }

  se_y_hat <- NA_real_
  vcov_unscaled <- matrix(NA_real_, K_beta, K_beta, dimnames = list(colnames(response_model_matrix_unscaled), colnames(response_model_matrix_unscaled)))
  A_condition <- tryCatch(kappa(analytical_jac_func(estimates)), error = function(e) NA_real_)

  vcov_try <- tryCatch({
    A_matrix <- analytical_jac_func(estimates)
    res <- el_compute_delta_variance(
      A_matrix = A_matrix,
      family = family,
      response_model_matrix_scaled = response_model_matrix_scaled,
      auxiliary_matrix_scaled = auxiliary_matrix_scaled,
      mu_x_scaled = mu_x_scaled,
      eta_i_hat = eta_i_hat,
      w_i_hat = w_i_hat,
      W_hat = W_hat,
      denominator_hat = denominator_hat,
      lambda_W_hat = lambda_W_hat,
      full_data = full_data,
      compute_score_variance_func = compute_score_variance_func,
      respondent_weights = respondent_weights,
      N_pop = N_pop,
      n_resp_weighted = n_resp_weighted,
      trim_cap = trim_cap,
      outcome_vec = outcome_vec,
      estimates = estimates
    )
    list(A = A_matrix, result = res, message = "Calculation successful")
  }, error = function(e) list(A = NULL, result = NULL, message = e$message))

  vcov_result <- vcov_try$result
  vcov_message <- vcov_try$message

  grad_source <- NA_character_
  diag_extra <- list(var_y_hat_val = NA_real_, var_anal2 = NA_real_, grad_l1 = NA_real_, sigma_min_eig = NA_real_)
  B_min_eig_diag <- NA_real_
  if (!is.null(vcov_result)) {
    var_y_hat_val <- as.numeric(vcov_result$var_y_hat)
    diag_extra$var_y_hat_val <- var_y_hat_val
# Analytic gradient against the same Sigma
    Sigma <- vcov_result$vcov_matrix_sandwich_scaled
# propagate B minimum eigenvalue when available
    B_min_eig_diag <- tryCatch(vcov_result$B_min_eig, error = function(e) NA_real_)
    var_anal2 <- NA_real_
    if (is.matrix(Sigma) && all(dim(Sigma) == length(estimates))) {
# record Sigma min eigenvalue (symmetric)
      diag_extra$sigma_min_eig <- tryCatch({
        ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
        suppressWarnings(min(ev))
      }, error = function(e) NA_real_)
      if (is.infinite(trim_cap)) {
        K_aux <- if (is.null(auxiliary_matrix_scaled) || ncol(auxiliary_matrix_scaled) == 0) 0 else ncol(auxiliary_matrix_scaled)
        X_beta <- response_model_matrix_scaled
        Xc <- if (K_aux > 0) sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-") else matrix(nrow = nrow(X_beta), ncol = 0)
        g_anal <- tryCatch(
          el_grad_g_analytic(
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
          ), error = function(e) NULL
        )
        if (is.numeric(g_anal) && length(g_anal) == ncol(Sigma) && all(is.finite(g_anal))) {
# Align gradient order to Sigma colnames if names available
          if (!is.null(names(g_anal)) && !is.null(colnames(Sigma))) {
            idx <- match(colnames(Sigma), names(g_anal))
            if (all(is.finite(idx))) g_anal <- g_anal[idx]
          }
          diag_extra$grad_l1 <- sum(abs(g_anal))
# Robust two-solve variance: solve t(A) x = g, then var = x' B x
          var_anal2 <- tryCatch({
            A_mat <- vcov_try$A
            if (is.null(A_mat)) stop("Jacobian missing for two-solve variance")
            x <- solve(t(A_mat), g_anal)
# Use B from vcov_result path (PSD by construction)
            as.numeric(t(x) %*% vcov_result$B %*% x)
          }, error = function(e) NA_real_)
          se_y_hat <- if (is.finite(var_anal2)) sqrt(pmax(var_anal2, 0)) else se_y_hat
          grad_source <- "analytic"
        }
      }
      diag_extra$var_anal2 <- var_anal2
      if (!is.finite(se_y_hat) && is.finite(var_y_hat_val)) se_y_hat <- sqrt(pmax(var_y_hat_val, 0))
    } else {
      if (is.finite(var_y_hat_val)) se_y_hat <- sqrt(pmax(var_y_hat_val, 0))
    }
# Robust vcov(beta): V_beta = X_beta' B X_beta with t(A) X_beta = E_beta
    vcov_beta_scaled <- NULL
    A_mat <- vcov_try$A
    B_mat <- vcov_result$B
    if (is.matrix(A_mat) && is.matrix(B_mat)) {
      p <- ncol(A_mat)
      if (p == nrow(A_mat) && p >= K_beta) {
        E_beta <- matrix(0, nrow = p, ncol = K_beta)
        E_beta[seq_len(K_beta), ] <- diag(K_beta)
        X_beta <- tryCatch(solve(t(A_mat), E_beta), error = function(e) NULL)
        if (is.matrix(X_beta) && nrow(X_beta) == p && ncol(X_beta) == K_beta) {
          Vb <- tryCatch({
            M <- crossprod(X_beta, B_mat %*% X_beta)
            0.5 * (M + t(M))
          }, error = function(e) NULL)
          if (is.matrix(Vb) && all(dim(Vb) == K_beta)) vcov_beta_scaled <- Vb
        }
      }
    }
# Fallback to Σ block if robust path failed
    if (is.null(vcov_beta_scaled) && is.matrix(Sigma)) {
      vcov_beta_scaled <- Sigma[1:K_beta, 1:K_beta, drop = FALSE]
      vcov_beta_scaled <- 0.5 * (vcov_beta_scaled + t(vcov_beta_scaled))
    }
    if (is.null(vcov_beta_scaled)) vcov_beta_scaled <- matrix(NA_real_, K_beta, K_beta)
# Ensure dimnames for beta block and names for scaled coefficients
    beta_names_scaled <- colnames(response_model_matrix_scaled)
    if (!is.null(beta_names_scaled) && all(length(beta_names_scaled) == K_beta)) {
      dimnames(vcov_beta_scaled) <- list(beta_names_scaled, beta_names_scaled)
    }
    coeffs_scaled <- estimates[1:K_beta]
    if (!is.null(beta_names_scaled) && length(coeffs_scaled) == K_beta) names(coeffs_scaled) <- beta_names_scaled
    vcov_unscaled <- if (standardize) unscale_coefficients(coeffs_scaled, vcov_beta_scaled, nmar_scaling_recipe)$vcov else vcov_beta_scaled
  } else {

  }
  if (!is.finite(se_y_hat)) se_y_hat <- NA_real_

  list(
    se_y_hat = se_y_hat,
    vcov_unscaled = vcov_unscaled,
    vcov_result = vcov_result,
    vcov_message = vcov_message,
    diag = list(
      A_condition = A_condition,
      grad_source = grad_source,
      var_y_hat_val = diag_extra$var_y_hat_val,
      var_anal2 = diag_extra$var_anal2,
      grad_l1 = diag_extra$grad_l1,
      sigma_min_eig = diag_extra$sigma_min_eig,
      B_min_eig = B_min_eig_diag
    )
  )
}

#' Analytic gradient of the mean functional g(theta)
#'
#' Computes the analytic gradient of the respondent-weighted mean functional
#'   g = sum_i(pi_i y_i)/sum_i(pi_i) under smooth conditions (no trimming),
#'   for the EL reparameterization (beta, z = logit(W), lambda_x). This mirrors the
#'   guarded denominators used in post-solution weight construction.
#'
#' @keywords internal
el_grad_g_analytic <- function(family,
                               X_beta,
                               Xc,
                               mu_x_scaled,
                               respondent_weights,
                               eta_i_hat,
                               w_i_hat,
                               W_hat,
                               denominator_hat,
                               lambda_W_hat,
                               outcome_vec,
                               n_resp_weighted,
                               N_pop) {
  K_beta <- ncol(X_beta)
  K_aux <- if (is.null(Xc) || ncol(Xc) == 0) 0 else ncol(Xc)
  denom <- pmax(as.numeric(denominator_hat), 1e-8)
  pi <- respondent_weights / denom
  B_sum <- sum(pi)
  if (!is.finite(B_sum) || B_sum <= 0) return(NULL)
  g_val <- sum(pi * outcome_vec) / B_sum
  mu_eta <- family$mu.eta(eta_i_hat)
# d denom / d beta = lambda_W * mu_eta * X
  dden_deta <- lambda_W_hat * mu_eta
  w_resp <- respondent_weights
  factor <- (outcome_vec - g_val) * (-w_resp) / (denom^2)
  grad_beta <- as.numeric(colSums(X_beta * as.numeric(dden_deta * factor))) / B_sum
# Grad wrt z (W): d denom / d z = d lambda_W / d z * (w_i - W) - lambda_W * dW/dz
  C_const <- (N_pop / n_resp_weighted) - 1
  dW_dz <- W_hat * (1 - W_hat)
  dlam_dz <- C_const * W_hat / (1 - W_hat) # = lambda_W * W/(1-W)
  dden_dz <- dlam_dz * (w_i_hat - W_hat) - lambda_W_hat * dW_dz
  grad_z <- sum((outcome_vec - g_val) * (-w_resp * dden_dz / (denom^2))) / B_sum
# Grad wrt lambda_x: d denom / d lambda_x = Xc
  if (K_aux > 0) {
    grad_lambda <- as.numeric(colSums(Xc * as.numeric((outcome_vec - g_val) * (-w_resp) / (denom^2)))) / B_sum
  } else {
    grad_lambda <- numeric(0)
  }
# Name and order gradient components consistently with A/B matrices
  nm_beta <- colnames(X_beta)
  nm_z <- "(W) (logit)"
  nm_lambda <- if (K_aux > 0) paste0("lambda_", colnames(Xc)) else character(0)
  out <- c(grad_beta, grad_z, grad_lambda)
  names(out) <- c(nm_beta, nm_z, nm_lambda)
  out
}
