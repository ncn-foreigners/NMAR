#' Analytical Jacobian for empirical likelihood
#' @details Builds the block Jacobian \eqn{A = \partial F/\partial \theta} for the
#'   EL system with \eqn{\theta = (\beta, z, \lambda_x)} and \eqn{z = \operatorname{logit}(W)}.
#'   Blocks follow Qin, Leung, and Shao (2002, Eqs. 7-10). The derivative with
#'   respect to the linear predictor for the missingness (response) model uses the Bernoulli score form
#'   \eqn{\partial/\partial\eta\, \log w(\eta) = \mu.\eta(\eta)/w(\eta)} with
#'   link-inverse clipping. Denominator guards are applied consistently when
#'   forming terms depending on \eqn{D_i(\theta)}.
#'
#'   Guarding policy (must remain consistent across equations/Jacobian/post):
#'   - Cap eta: eta <- pmax(pmin(eta, get_eta_cap()), -get_eta_cap())
#'   - Compute w <- family$linkinv(eta); clip to [1e-12, 1-1e-12] when used in ratios
#'   - Denominator floor: Di <- pmax(Di_raw, nmar_get_el_denom_floor());
#'     multiply terms that depend on d(1/Di)/d(.) by active = 1(Di_raw > floor)
#'
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193-200.
#'
#' @keywords internal
el_build_jacobian <- function(family, missingness_model_matrix, auxiliary_matrix,
                              respondent_weights, N_pop, n_resp_weighted, mu_x_scaled) {
  force(family)
  force(missingness_model_matrix)
  force(auxiliary_matrix)
  force(respondent_weights)
  force(N_pop)
  force(n_resp_weighted)
  force(mu_x_scaled)
  if (is.null(family) || is.null(family$d2mu.deta2)) {
    return(NULL)
  }
  if (is.null(missingness_model_matrix) || !is.matrix(missingness_model_matrix)) stop("missingness_model_matrix must be a matrix.")
  n_resp <- nrow(missingness_model_matrix)
  K_beta <- ncol(missingness_model_matrix)
  K_aux <- if (is.null(auxiliary_matrix) || ncol(auxiliary_matrix) == 0) 0 else ncol(auxiliary_matrix)
  if (length(respondent_weights) != n_resp) stop("Length of respondent_weights must equal nrow(missingness_model_matrix).")
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
# Precompute centered auxiliary design once (avoids per-call sweep)
  X_centered_base <- if (K_aux > 0) sweep(auxiliary_matrix_mat, 2, mu_x_scaled_vec, "-") else NULL
  C_const <- (N_pop / n_resp_weighted) - 1
  ETA_CAP <- get_eta_cap()
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
    eta_raw <- as.vector(missingness_model_matrix %*% beta_vec)
    eta_state <- el_core_eta_state(family, eta_raw, ETA_CAP)
    eta_i <- eta_state$eta
    w_i <- eta_state$w
    mu_eta_i <- eta_state$mu_eta
    d2mu_eta2_i <- eta_state$d2mu
    X_centered <- if (K_aux > 0) X_centered_base else NULL
# QLS Eq. (5): Di = 1 + lambda_W * (w_i - W) + (Xc %*% lambda_x)
    Xc_lambda <- if (K_aux > 0) as.vector(X_centered %*% lambda_x) else 0
    dpack <- el_denominator(lambda_W, W_bounded, Xc_lambda, w_i, nmar_get_el_denom_floor())
    denominator <- dpack$denom
    inv_denom <- dpack$inv
    inv_denom_sq <- dpack$inv_sq
    active <- dpack$active
    dden_dTheta <- active * (d_lambda_W_dTheta * (w_i - W_bounded) - lambda_W * dWb_dTheta)
    dden_deta <- active * (lambda_W * mu_eta_i)
# Score wrt eta and its derivative (via shared core helper)
    s_eta_i <- eta_state$s_eta
    ds_eta_deta <- eta_state$ds_eta_deta
# beta block term (QLS Eq. 9)
    beta_eq_term <- s_eta_i - lambda_W * mu_eta_i * inv_denom
# Respect denominator floor: terms depending on d(1/Di)/deta vanish when clamped
    d_betaeq_deta <- ds_eta_deta - lambda_W * d2mu_eta2_i * inv_denom + active * (lambda_W^2) * (mu_eta_i^2) * inv_denom_sq
    d_betaeq_dTheta <- -d_lambda_W_dTheta * mu_eta_i * inv_denom + lambda_W * mu_eta_i * inv_denom_sq * dden_dTheta
# Build blocks with preallocation to avoid cbind/rbind overhead
    p_dim <- K_beta + 1 + K_aux
    full_mat <- matrix(0, nrow = p_dim, ncol = p_dim)
# Indices
    idx_beta <- seq_len(K_beta)
    idx_W <- K_beta + 1L
    idx_lambda <- if (K_aux > 0) (K_beta + 2L):(K_beta + 1L + K_aux) else integer(0)

# J11: d beta eqs / d beta
    w_eff_11 <- as.numeric(respondent_weights * d_betaeq_deta)
    full_mat[idx_beta, idx_beta] <- shared_weighted_gram(missingness_model_matrix, w_eff_11)
# J12: d beta eqs / d W
    j12_vec <- shared_weighted_Xty(missingness_model_matrix, respondent_weights, d_betaeq_dTheta)
    full_mat[idx_beta, idx_W] <- as.numeric(j12_vec)
# J13: d beta eqs / d lambda
    if (K_aux > 0) {
# Respect denominator floor via 'active' mask
      d_betaeq_dlambda_mat <- active * lambda_W * mu_eta_i * inv_denom_sq * X_centered
      full_mat[idx_beta, idx_lambda] <- shared_weighted_XtY(missingness_model_matrix, respondent_weights, as.matrix(d_betaeq_dlambda_mat))
    }
# J21: d W eq / d beta
    term21 <- mu_eta_i * inv_denom - (w_i - W_bounded) * inv_denom_sq * (lambda_W * mu_eta_i) * active
    full_mat[idx_W, idx_beta] <- as.numeric(t(shared_weighted_Xty(missingness_model_matrix, respondent_weights, term21)))
# J22: d W eq / d W
    term22 <- -dWb_dTheta * inv_denom - (w_i - W_bounded) * inv_denom_sq * dden_dTheta
    full_mat[idx_W, idx_W] <- as.numeric(crossprod(respondent_weights, term22))
# J23: d W eq / d lambda
    if (K_aux > 0) {
      j23 <- t(shared_weighted_Xty(X_centered, respondent_weights, (-(w_i - W_bounded) * inv_denom_sq * active)))
      full_mat[idx_W, idx_lambda] <- as.numeric(j23)
    }
    if (K_aux > 0) {
# J31: d aux eq / d beta
      term31 <- -dden_deta * inv_denom_sq
      full_mat[idx_lambda, idx_beta] <- shared_weighted_XtY(X_centered, as.numeric(respondent_weights * term31), missingness_model_matrix)
# J32: d aux eq / d W
      term32 <- -dden_dTheta * inv_denom_sq
      full_mat[idx_lambda, idx_W] <- as.numeric(shared_weighted_Xty(X_centered, respondent_weights, term32))
# J33: d aux eq / d lambda (SPD, negative definite of weighted Gram)
      full_mat[idx_lambda, idx_lambda] <- -shared_weighted_gram(X_centered, as.numeric(respondent_weights * (inv_denom_sq * active)))
    }
# Optional names (kept minimal to reduce overhead)
    param_names <- c(colnames(missingness_model_matrix), "(W) (logit)", if (K_aux > 0) paste0("lambda_", colnames(auxiliary_matrix_mat)) else NULL)
    if (!is.null(param_names) && length(param_names) == ncol(full_mat)) {
      colnames(full_mat) <- rownames(full_mat) <- param_names
    }
    full_mat
  }
}

#' Analytical Jacobian for survey EL system (design-weighted QLS analogue)
#'
#' @details Builds the block Jacobian \eqn{A = \partial g/\partial \theta} for the
#'   survey EL system with \eqn{\theta = (\beta, z, \lambda_W, \lambda_x)} and
#'   \eqn{z = \operatorname{logit}(W)}. Blocks follow the design-weighted analogue
#'   of Qin, Leung, and Shao (2002) used in \code{el_build_equation_system_survey()}.
#'   Guarding policy matches the IID Jacobian:
#'   \itemize{
#'     \item cap eta: \code{eta <- pmax(pmin(eta, get_eta_cap()), -get_eta_cap())}
#'     \item compute \code{w <- family$linkinv(eta)} and clip to \code{[1e-12, 1-1e-12]}
#'       when used in ratios
#'     \item denominator floor: \code{Di <- pmax(Di_raw, nmar_get_el_denom_floor())};
#'       multiply terms depending on \code{d(1/Di)/d(.)} by \code{active = 1(Di_raw > floor)}
#'   }
#'
#'   The Jacobian uses the same score and second-derivative machinery as
#'   \code{el_build_jacobian()}; when \code{family$d2mu.deta2} is missing, this
#'   function returns \code{NULL} and the solver falls back to numeric/Broyden
#'   Jacobians.
#'
#' @keywords internal
el_build_jacobian_survey <- function(family, missingness_model_matrix, auxiliary_matrix,
                                     respondent_weights, N_pop, n_resp_weighted, mu_x_scaled) {
  force(family)
  force(missingness_model_matrix)
  force(auxiliary_matrix)
  force(respondent_weights)
  force(N_pop)
  force(n_resp_weighted)
  force(mu_x_scaled)
  if (is.null(family) || is.null(family$d2mu.deta2)) {
    return(NULL)
  }
  if (is.null(missingness_model_matrix) || !is.matrix(missingness_model_matrix)) {
    stop("missingness_model_matrix must be a matrix.")
  }
  n_resp <- nrow(missingness_model_matrix)
  K_beta <- ncol(missingness_model_matrix)
  K_aux <- if (is.null(auxiliary_matrix) || ncol(auxiliary_matrix) == 0) 0 else ncol(auxiliary_matrix)
  if (length(respondent_weights) != n_resp) {
    stop("Length of respondent_weights must equal nrow(missingness_model_matrix).")
  }
  if (K_aux == 0) {
    auxiliary_matrix_mat <- matrix(nrow = n_resp, ncol = 0)
    mu_x_scaled_vec <- numeric(0)
  } else {
    auxiliary_matrix_mat <- as.matrix(auxiliary_matrix)
    if (nrow(auxiliary_matrix_mat) != n_resp) stop("auxiliary_matrix must match rows.")
    if (is.null(names(mu_x_scaled))) stop("mu_x_scaled must be named.")
    if (!all(colnames(auxiliary_matrix_mat) %in% names(mu_x_scaled))) {
      stop("Names of mu_x_scaled must include all columns of auxiliary_matrix.")
    }
    mu_x_scaled_vec <- as.numeric(mu_x_scaled[colnames(auxiliary_matrix_mat)])
    names(mu_x_scaled_vec) <- colnames(auxiliary_matrix_mat)
  }
  X_centered_base <- if (K_aux > 0) sweep(auxiliary_matrix_mat, 2, mu_x_scaled_vec, "-") else NULL
  T0 <- N_pop - n_resp_weighted
  ETA_CAP <- get_eta_cap()
  function(params) {
    expected_len <- K_beta + 1L + 1L + K_aux
    if (length(params) != expected_len) {
      stop("Parameter vector length mismatch in survey EL Jacobian.", call. = FALSE)
    }
    beta_vec <- as.numeric(params[1:K_beta])
    z <- as.numeric(params[K_beta + 1L])
    lambda_W <- as.numeric(params[K_beta + 2L])
    lambda_x <- if (K_aux > 0) as.numeric(params[(K_beta + 3L):expected_len]) else numeric(0)
    W <- plogis(z)
    W <- min(max(W, 1e-12), 1 - 1e-12)
    W_bounded <- W
    dW_dz <- if (W > 1e-12 && W < 1 - 1e-12) W * (1 - W) else 0
    eta_raw <- as.vector(missingness_model_matrix %*% beta_vec)
    eta_state <- el_core_eta_state(family, eta_raw, ETA_CAP)
    eta_i <- eta_state$eta
    w_i <- eta_state$w
    mu_eta_i <- eta_state$mu_eta
    d2mu_eta2_i <- eta_state$d2mu
    X_centered <- if (K_aux > 0) X_centered_base else NULL
    Xc_lambda <- if (K_aux > 0) as.vector(X_centered %*% lambda_x) else 0
    dpack <- el_denominator(lambda_W, W_bounded, Xc_lambda, w_i, nmar_get_el_denom_floor())
    denominator <- dpack$denom
    inv_denom <- dpack$inv
    inv_denom_sq <- dpack$inv_sq
    active <- dpack$active

# Score wrt eta and its derivative (same pattern as IID Jacobian)
    s_eta_i <- eta_state$s_eta
    ds_eta_deta <- eta_state$ds_eta_deta

# Precompute some terms
    w_minus_W <- w_i - W_bounded
# g_beta block: g_beta = sum d_i x_i [s_eta_i - lambda_W * mu_eta_i / D_i]
    beta_eq_term <- s_eta_i - lambda_W * mu_eta_i * inv_denom
# d b_i / d eta_i
    d_betaeq_deta <- ds_eta_deta - lambda_W * d2mu_eta2_i * inv_denom +
      active * (lambda_W^2) * (mu_eta_i^2) * inv_denom_sq
# d b_i / d z
    d_betaeq_dz <- -active * (lambda_W^2) * mu_eta_i * dW_dz * inv_denom_sq
# d b_i / d lambda_W
    d_betaeq_dlambdaW <- -mu_eta_i * inv_denom +
      active * lambda_W * mu_eta_i * w_minus_W * inv_denom_sq
# d b_i / d lambda_x (vector per i)
    if (K_aux > 0) {
      d_betaeq_dlambda_mat <- active * lambda_W * mu_eta_i * inv_denom_sq * X_centered
    }

    p_dim <- K_beta + 1L + 1L + K_aux
    full_mat <- matrix(0, nrow = p_dim, ncol = p_dim)
    idx_beta <- seq_len(K_beta)
    idx_z <- K_beta + 1L
    idx_lambdaW <- K_beta + 2L
    idx_lambda_x <- if (K_aux > 0) (K_beta + 3L):p_dim else integer(0)

# J_beta,beta
    w_eff_11 <- as.numeric(respondent_weights * d_betaeq_deta)
    full_mat[idx_beta, idx_beta] <- shared_weighted_gram(missingness_model_matrix, w_eff_11)
# J_beta,z
    j_beta_z <- shared_weighted_Xty(missingness_model_matrix, respondent_weights, d_betaeq_dz)
    full_mat[idx_beta, idx_z] <- as.numeric(j_beta_z)
# J_beta,lambda_W
    j_beta_lambdaW <- shared_weighted_Xty(missingness_model_matrix, respondent_weights, d_betaeq_dlambdaW)
    full_mat[idx_beta, idx_lambdaW] <- as.numeric(j_beta_lambdaW)
# J_beta,lambda_x
    if (K_aux > 0) {
      full_mat[idx_beta, idx_lambda_x] <- shared_weighted_XtY(
        missingness_model_matrix,
        respondent_weights,
        as.matrix(d_betaeq_dlambda_mat)
      )
    }

# g_W block: g_W = sum d_i (w_i - W)/D_i
# J_W,beta
    term_W_beta <- mu_eta_i * inv_denom - active * w_minus_W * inv_denom_sq * (lambda_W * mu_eta_i)
    full_mat[idx_z, idx_beta] <- as.numeric(
      t(shared_weighted_Xty(missingness_model_matrix, respondent_weights, term_W_beta))
    )
# J_W,z
    term_W_z <- dW_dz * (active * lambda_W * w_minus_W * inv_denom_sq - inv_denom)
    full_mat[idx_z, idx_z] <- as.numeric(crossprod(respondent_weights, term_W_z))
# J_W,lambda_W
    term_W_lambdaW <- -active * (w_minus_W^2) * inv_denom_sq
    full_mat[idx_z, idx_lambdaW] <- as.numeric(crossprod(respondent_weights, term_W_lambdaW))
# J_W,lambda_x
    if (K_aux > 0) {
      j_W_lambda_x <- t(shared_weighted_Xty(
        X_centered,
        respondent_weights,
        -(w_minus_W * inv_denom_sq * active)
      ))
      full_mat[idx_z, idx_lambda_x] <- as.numeric(j_W_lambda_x)
    }

    if (K_aux > 0) {
# g_aux block: g_aux = sum d_i t_i / D_i
# J_aux,beta
      term_aux_beta <- -active * lambda_W * mu_eta_i * inv_denom_sq
      full_mat[idx_lambda_x, idx_beta] <- shared_weighted_XtY(
        X_centered,
        as.numeric(respondent_weights * term_aux_beta),
        missingness_model_matrix
      )
# J_aux,z
      term_aux_z <- active * lambda_W * dW_dz * inv_denom_sq
      full_mat[idx_lambda_x, idx_z] <- as.numeric(
        shared_weighted_Xty(X_centered, respondent_weights, term_aux_z)
      )
# J_aux,lambda_W
      term_aux_lambdaW <- -active * w_minus_W * inv_denom_sq
      full_mat[idx_lambda_x, idx_lambdaW] <- as.numeric(
        shared_weighted_Xty(X_centered, respondent_weights, term_aux_lambdaW)
      )
# J_aux,lambda_x
      full_mat[idx_lambda_x, idx_lambda_x] <- -shared_weighted_gram(
        X_centered,
        as.numeric(respondent_weights * (inv_denom_sq * active))
      )
    }

# g_link block: g_link = T0/(1-W) - lambda_W * sum d_i / D_i
    S <- sum(respondent_weights * inv_denom)
# J_link,beta
    term_link_beta <- lambda_W^2 * mu_eta_i * inv_denom_sq * active
    full_mat[p_dim, idx_beta] <- as.numeric(
      t(shared_weighted_Xty(missingness_model_matrix, respondent_weights, term_link_beta))
    )
# J_link,z: d/dz[T0/(1-W) - lambda_W * S(z)]
# First term: T0 * W / (1-W); second term: -lambda_W * dS/dz with
# dS/dz = sum d_i * lambda_W * dW_dz * active / D_i^2.
    term_link_z_S <- lambda_W^2 * dW_dz * inv_denom_sq * active
    dS_dz <- sum(respondent_weights * term_link_z_S)
    full_mat[p_dim, idx_z] <- T0 * W_bounded / (1 - W_bounded) - dS_dz
# J_link,lambda_W
    term_dS_dlambdaW <- -active * w_minus_W * inv_denom_sq
    dS_dlambdaW <- sum(respondent_weights * term_dS_dlambdaW)
    full_mat[p_dim, idx_lambdaW] <- -S - lambda_W * dS_dlambdaW
# J_link,lambda_x
    if (K_aux > 0) {
      j_link_lambda_x <- lambda_W * shared_weighted_Xty(
        X_centered,
        respondent_weights,
        inv_denom_sq * active
      )
      full_mat[p_dim, idx_lambda_x] <- as.numeric(j_link_lambda_x)
    }

    param_names <- c(
      colnames(missingness_model_matrix),
      "(W) (logit)",
      "(lambda_W)",
      if (K_aux > 0) paste0("lambda_", colnames(auxiliary_matrix_mat)) else NULL
    )
    if (!is.null(param_names) && length(param_names) == ncol(full_mat)) {
      colnames(full_mat) <- rownames(full_mat) <- param_names
    }
    full_mat
  }
}
