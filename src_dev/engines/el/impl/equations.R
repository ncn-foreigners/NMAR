#' Empirical likelihood estimating equations for SRS
#'
#' Returns a function that evaluates the stacked EL system for
#' \eqn{\theta = (\beta, z, \lambda_x)} with \eqn{z = \operatorname{logit}(W)}.
#' Blocks correspond to:
#' \enumerate{
#' \item missingness model score equations in \eqn{\beta},
#' \item the response-rate equation in \eqn{W},
#' \item auxiliary moment constraints in \eqn{\lambda_x}.
#' }
#'
#' When no auxiliaries are present the last block is omitted. The system
#' matches QLS equations 7-10. We cap \eqn{\eta}, clip \eqn{w_i} in ratios,
#' and guard \eqn{D_i} away from zero to ensure numerical stability.
#'
#' \strong{Guarding policy:}
#' \itemize{
#' \item Cap \eqn{\eta}:
#' \code{eta <- pmax(pmin(eta, get_eta_cap()), -get_eta_cap())}.
#' \item Compute \code{w <- family$linkinv(eta)} and clip to
#' \code{[1e-12, 1 - 1e-12]} when used in ratios.
#' \item Denominator floor:
#' \code{Di <- pmax(Di_raw, nmar_get_el_denom_floor())}. In the Jacobian,
#' terms that depend on \code{d(1/Di)/d(.)} are multiplied by
#' \code{active = 1(Di_raw > floor)} to match the clamped equations.
#' }
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

  X_centered <- if (K_aux > 0) {
    sweep(auxiliary_matrix, 2, mu_x_scaled, "-")
    } else {
      matrix(nrow = nrow(missingness_model_matrix), ncol = 0)
    }

  C_const <- (N_pop / n_resp_weighted) - 1
  ETA_CAP <- get_eta_cap()

  function(params) {
    beta_vec <- params[1:K_beta]
    z <- params[K_beta + 1]
    W <- plogis(z)
    W <- min(max(W, 1e-12), 1 - 1e-12)
    lambda_x <- if (K_aux > 0) params[(K_beta + 2):length(params)] else numeric(0)
    W_bounded <- W
# QLS eq. 10: lambda_W = (N/n - 1) / (1 - W)
    lambda_W <- el_lambda_W(C_const, W_bounded)
    eta_raw <- as.vector(missingness_model_matrix %*% beta_vec)
    eta_state <- el_core_eta_state(family, eta_raw, ETA_CAP)
    w_i <- eta_state$w
    mu_eta_i <- eta_state$mu_eta
    s_eta_i <- eta_state$s_eta
# QLS eq. 5: Di = 1 + lambda_W * (w_i - W) + (Xc %*% lambda_x)
    Xc_lambda <- if (K_aux > 0) as.vector(X_centered %*% lambda_x) else 0
    dpack <- el_denominator(lambda_W, W_bounded, Xc_lambda, w_i, nmar_get_el_denom_floor())
    inv_denominator <- dpack$inv
# beta block (QLS Eq. 9): s_eta(eta) - lambda_W * mu.eta(eta) / Di
    beta_eq_term <- s_eta_i - lambda_W * mu_eta_i * inv_denominator
    eq_betas <- shared_weighted_Xty(missingness_model_matrix, respondent_weights, beta_eq_term)
# W equation (QLS eq. 8)
    eq_W <- as.numeric(crossprod(respondent_weights * inv_denominator, (w_i - W_bounded)))
# Auxiliary constraints (QLS eq. 7)
    if (K_aux > 0) {
      eq_constraints <- shared_weighted_Xty(X_centered, respondent_weights, inv_denominator)
    } else {
      eq_constraints <- numeric(0)
    }
    c(as.vector(eq_betas), eq_W, as.vector(eq_constraints))
  }
}

#' Empirical likelihood equations for survey designs
#'
#' Returns a function that evaluates the stacked EL system for survey designs
#' using design weights. Unknowns are
#' \eqn{\theta = (\beta, z, \lambda_W, \lambda_x)} with \eqn{z = \operatorname{logit}(W)}.
#' Blocks correspond to:
#' \itemize{
#' \item response-model score equations in \eqn{\beta},
#' \item the response-rate equation in \eqn{W} based on \eqn{\sum d_i (w_i - W)/D_i = 0},
#' \item auxiliary moment constraints \eqn{\sum d_i (X_i - \mu_x)/D_i = 0},
#' \item and the design-based linkage between \eqn{\lambda_W} and the
#' nonrespondent total: \eqn{T_0/(1-W) - \lambda_W \sum d_i / D_i = 0},
#' where \eqn{T_0 = N_{\mathrm{pop}} - \sum d_i} on the analysis scale.
#' }
#'
#' When all design weights are equal and \eqn{N_{\mathrm{pop}}} and the respondent
#' count match the simple random sampling setup, this system reduces to the QLS equations 6-10.
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

  X_centered <- if (K_aux > 0) {
    sweep(auxiliary_matrix, 2, mu_x_scaled, "-")
    } else {
      matrix(nrow = nrow(missingness_model_matrix), ncol = 0)
    }

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
    eta_state <- el_core_eta_state(family, eta_raw, ETA_CAP)
    w_i <- eta_state$w
    mu_eta_i <- eta_state$mu_eta
    s_eta_i <- eta_state$s_eta
# Denominator Di = 1 + lambda_W (w_i - W) + (X_centered %*% lambda_x)
    Xc_lambda <- if (K_aux > 0) as.vector(X_centered %*% lambda_x) else 0
    dpack <- el_denominator(lambda_W, W_bounded, Xc_lambda, w_i, denom_floor)
    inv_denominator <- dpack$inv
# beta block: s_eta - lambda_W * mu_eta / D_i
    beta_eq_term <- s_eta_i - lambda_W * mu_eta_i * inv_denominator
    eq_betas <- shared_weighted_Xty(missingness_model_matrix, respondent_weights, beta_eq_term)
# W constraint: sum d_i (w_i - W) / D_i = 0
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
