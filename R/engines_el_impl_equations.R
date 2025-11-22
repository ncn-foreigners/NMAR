#' Empirical likelihood estimating equations
#' @details Returns a function that evaluates the stacked EL system for
#'   \eqn{\theta = (\beta, z, \lambda_x)} with \eqn{z = \operatorname{logit}(W)}.
#'   Blocks correspond to: (i) missingness (response) model score equations in \eqn{\beta},
#'   (ii) the response-rate equation in \eqn{W}, and (iii) auxiliary moment
#'   constraints in \eqn{\lambda_x}. When no auxiliaries are present the last
#'   block is omitted. The system matches Qin, Leung, and Shao (2002, Eqs. 7-10)
#'   with empirical masses \eqn{m_i = d_i/D_i(\theta)}, \eqn{D_i} as in the paper.
#'   We cap \eqn{\eta}, clip \eqn{p}, and guard \eqn{D_i} away from zero to
#'   ensure numerical stability; these safeguards are applied consistently in
#'   equations, Jacobian, and post-solution weights.
#'
#'   Guarding policy (must remain consistent across equations/Jacobian/post):
#'   - Cap eta: eta <- pmax(pmin(eta, get_eta_cap()), -get_eta_cap())
#'   - Compute w <- family$linkinv(eta); clip to [1e-12, 1-1e-12] when used in ratios
#'   - Denominator floor: Di <- pmax(Di_raw, nmar_get_el_denom_floor()); in the
#'     Jacobian, multiply terms that depend on d(1/Di)/d(.) by
#'     active = 1(Di_raw > floor)
#'
#'   The score with respect to the linear predictor uses the Bernoulli form
#'   \eqn{s_{\eta,i}(\beta) = \partial \log w_i / \partial \eta_i = \mu.\eta(\eta_i)/w_i},
#'   which is valid for both logit and probit links when \eqn{w_i} is clipped.
#'
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193-200.
#'
#' @keywords internal
el_build_equation_system <- function(family, missingness_model_matrix, auxiliary_matrix,
                                     respondent_weights, N_pop, n_resp_weighted, mu_x_scaled) {
  force(family)
  force(missingness_model_matrix)
  force(auxiliary_matrix)
  force(respondent_weights)
  force(N_pop)
  force(n_resp_weighted)
  force(mu_x_scaled)
  K_beta <- ncol(missingness_model_matrix)
  K_aux <- if (is.null(auxiliary_matrix) || ncol(auxiliary_matrix) == 0) 0 else ncol(auxiliary_matrix)
# Hoist centered auxiliaries and constants outside parameter closure
  X_centered <- if (K_aux > 0) sweep(auxiliary_matrix, 2, mu_x_scaled, "-") else matrix(nrow = nrow(missingness_model_matrix), ncol = 0)
  C_const <- (N_pop / n_resp_weighted) - 1
  ETA_CAP <- get_eta_cap()
  function(params) {
    beta_vec <- params[1:K_beta]
    z <- params[K_beta + 1]
    W <- plogis(z)
    W <- min(max(W, 1e-12), 1 - 1e-12)
    lambda_x <- if (K_aux > 0) params[(K_beta + 2):length(params)] else numeric(0)
    W_bounded <- W
# QLS Eq. (10): lambda_W = (N/n - 1) / (1 - W)
    lambda_W <- el_lambda_W(C_const, W_bounded)
    eta_raw <- as.vector(missingness_model_matrix %*% beta_vec)
    eta_i <- pmax(pmin(eta_raw, ETA_CAP), -ETA_CAP)
    w_i <- family$linkinv(eta_i)
    mu_eta_i <- family$mu.eta(eta_i)
# Unified score w.r.t. eta for delta=1: prefer numerically stable family score
    w_i_clipped <- pmin(pmax(w_i, 1e-12), 1 - 1e-12)
    if (!is.null(family$name) && identical(family$name, "logit")) {
# For logit, s_eta = 1 - w
      s_eta_i <- 1 - w_i_clipped
    } else if (!is.null(family$name) && identical(family$name, "probit")) {
# For probit, use Mills ratio in log domain
      log_phi <- stats::dnorm(eta_i, log = TRUE)
      log_Phi <- stats::pnorm(eta_i, log.p = TRUE)
      s_eta_i <- exp(log_phi - log_Phi)
    } else if (!is.null(family$score_eta)) {
      s_eta_i <- family$score_eta(eta_i, 1)
    } else {
      s_eta_i <- mu_eta_i / w_i_clipped
    }
# QLS Eq. (5): Di = 1 + lambda_W * (w_i - W) + (Xc %*% lambda_x)
    Xc_lambda <- if (K_aux > 0) as.vector(X_centered %*% lambda_x) else 0
    dpack <- el_denominator(lambda_W, W_bounded, Xc_lambda, w_i, nmar_get_el_denom_floor())
    inv_denominator <- dpack$inv
# beta block (QLS Eq. 9): s_eta(eta) - lambda_W * mu.eta(eta) / Di
    beta_eq_term <- s_eta_i - lambda_W * mu_eta_i * inv_denominator
    eq_betas <- shared_weighted_Xty(missingness_model_matrix, respondent_weights, beta_eq_term)
# W equation (QLS Eq. 8)
    eq_W <- as.numeric(crossprod(respondent_weights * inv_denominator, (w_i - W_bounded)))
# Auxiliary constraints (QLS Eq. 7)
    if (K_aux > 0) {
      eq_constraints <- shared_weighted_Xty(X_centered, respondent_weights, inv_denominator)
    } else {
# Avoid allocating zero-width matrices repeatedly
      eq_constraints <- numeric(0)
    }
    c(as.vector(eq_betas), eq_W, as.vector(eq_constraints))
  }
}

#' Empirical likelihood equations for survey designs (design-weighted QLS system)
#'
#' @details
#'   Returns a function that evaluates the stacked EL system for complex survey
#'   designs using design weights. Unknowns are
#'   \eqn{\theta = (\beta, z, \lambda_W, \lambda_x)} with \eqn{z = \operatorname{logit}(W)}.
#'   Blocks correspond to:
#'   \itemize{
#'     \item response-model score equations in \eqn{\beta},
#'     \item the response-rate equation in \eqn{W} based on \eqn{\sum d_i (w_i - W)/D_i = 0},
#'     \item auxiliary moment constraints \eqn{\sum d_i (X_i - \mu_x)/D_i = 0},
#'     \item and the design-based linkage between \eqn{\lambda_W} and the
#'       nonrespondent total: \eqn{T_0/(1-W) - \lambda_W \sum d_i / D_i = 0},
#'       where \eqn{T_0 = N_{\mathrm{pop}} - \sum d_i} on the analysis scale.
#'   }
#'
#'   When all design weights are equal and \eqn{N_{\mathrm{pop}}} and the respondent
#'   count match the simple random sampling setup, this system reduces to the
#'   Qin, Leung, and Shao (2002) equations (6)-(10).
#'
#' @keywords internal
el_build_equation_system_survey <- function(family, missingness_model_matrix, auxiliary_matrix,
                                            respondent_weights, N_pop, n_resp_weighted, mu_x_scaled) {
  force(family)
  force(missingness_model_matrix)
  force(auxiliary_matrix)
  force(respondent_weights)
  force(N_pop)
  force(n_resp_weighted)
  force(mu_x_scaled)
  K_beta <- ncol(missingness_model_matrix)
  K_aux <- if (is.null(auxiliary_matrix) || ncol(auxiliary_matrix) == 0) 0 else ncol(auxiliary_matrix)
# Center auxiliaries once
  X_centered <- if (K_aux > 0) sweep(auxiliary_matrix, 2, mu_x_scaled, "-") else matrix(nrow = nrow(missingness_model_matrix), ncol = 0)
# Design-based nonrespondent total on analysis scale
  T0 <- N_pop - n_resp_weighted
  ETA_CAP <- get_eta_cap()
  denom_floor <- nmar_get_el_denom_floor()
  function(params) {
    expected_len <- K_beta + 1L + 1L + K_aux
    if (length(params) != expected_len) {
      stop("Internal error: parameter length mismatch in survey EL equation system.", call. = FALSE)
    }
    beta_vec <- params[1:K_beta]
    z <- params[K_beta + 1L]
    lambda_W <- params[K_beta + 2L]
    lambda_x <- if (K_aux > 0) params[(K_beta + 3L):expected_len] else numeric(0)
    W <- plogis(z)
    W <- min(max(W, 1e-12), 1 - 1e-12)
    W_bounded <- W
    eta_raw <- as.vector(missingness_model_matrix %*% beta_vec)
    eta_i <- pmax(pmin(eta_raw, ETA_CAP), -ETA_CAP)
    w_i <- family$linkinv(eta_i)
    mu_eta_i <- family$mu.eta(eta_i)
# Stable score w.r.t eta
    w_i_clipped <- pmin(pmax(w_i, 1e-12), 1 - 1e-12)
    if (!is.null(family$name) && identical(family$name, "logit")) {
      s_eta_i <- 1 - w_i_clipped
    } else if (!is.null(family$name) && identical(family$name, "probit")) {
      log_phi <- stats::dnorm(eta_i, log = TRUE)
      log_Phi <- stats::pnorm(eta_i, log.p = TRUE)
      s_eta_i <- exp(log_phi - log_Phi)
    } else if (!is.null(family$score_eta)) {
      s_eta_i <- family$score_eta(eta_i, 1)
    } else {
      s_eta_i <- mu_eta_i / w_i_clipped
    }
# Denominator Di = 1 + lambda_W (w_i - W) + (X_centered %*% lambda_x)
    Xc_lambda <- if (K_aux > 0) as.vector(X_centered %*% lambda_x) else 0
    dpack <- el_denominator(lambda_W, W_bounded, Xc_lambda, w_i, denom_floor)
    inv_denominator <- dpack$inv
# beta block: s_eta - lambda_W * mu_eta / D_i
    beta_eq_term <- s_eta_i - lambda_W * mu_eta_i * inv_denominator
    eq_betas <- shared_weighted_Xty(missingness_model_matrix, respondent_weights, beta_eq_term)
# W equation (constraint sum): sum d_i (w_i - W) / D_i = 0
    eq_W_constraint <- as.numeric(crossprod(respondent_weights * inv_denominator, (w_i - W_bounded)))
# Auxiliary constraints: sum d_i (X - mu_x) / D_i = 0
    if (K_aux > 0) {
      eq_constraints <- shared_weighted_Xty(X_centered, respondent_weights, inv_denominator)
    } else {
      eq_constraints <- numeric(0)
    }
# W-lambda_W linkage: T0/(1-W) - lambda_W * sum d_i / D_i = 0
    sum_d_over_D <- as.numeric(crossprod(respondent_weights, inv_denominator))
    eq_W_link <- (T0 / (1 - W_bounded)) - lambda_W * sum_d_over_D
    c(as.vector(eq_betas), eq_W_constraint, as.vector(eq_constraints), eq_W_link)
  }
}
