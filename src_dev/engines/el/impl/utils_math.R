#' Core computations
#'
#' Computes the capped linear predictor, response probabilities, derivatives,
#' and stable scores with respect to the linear predictor for a given family.
#' Centralizes the numerically delicate pieces (capping, clipping, score derivatives)
#' to be reused in EL equations and jacobian.
#'
#' @param family Response family.
#' @param eta_raw Numeric vector of unconstrained linear predictors.
#' @param eta_cap Scalar cap applied symmetrically to \code{eta_raw}.
#'
#' @return A list with components:
#' \describe{
#' \item{\code{eta}}{Capped linear predictor.}
#' \item{\code{w}}{Mean function \code{family$linkinv(eta)}.}
#' \item{\code{w_clipped}}{\code{w} clipped to \code{[1e-12, 1-1e-12]} for use in ratios.}
#' \item{\code{mu_eta}}{Derivative \code{family$mu.eta(eta)}.}
#' \item{\code{d2mu}}{Second derivative \code{family$d2mu.deta2(eta)}
#' when available, otherwise \code{NULL}.}
#' \item{\code{s_eta}}{Score with respect to \code{eta}.}
#' \item{\code{ds_eta_deta}}{Derivative of \code{s_eta} with respect to \code{eta}
#' when \code{d2mu} is available, otherwise \code{NULL}.}
#' }
#'
#' @keywords internal
el_core_eta_state <- function(family, eta_raw, eta_cap) {
  eta <- pmax(pmin(eta_raw, eta_cap), -eta_cap)

  w <- family$linkinv(eta)
  mu_eta <- family$mu.eta(eta)
  d2mu <- if (!is.null(family$d2mu.deta2)) family$d2mu.deta2(eta) else NULL

# Clip probabilities for ratios
  w_clipped <- pmin(pmax(w, 1e-12), 1 - 1e-12)

  is_logit <- !is.null(family$name) && identical(family$name, "logit")
  is_probit <- !is.null(family$name) && identical(family$name, "probit")

# Score wrt eta
  mills_ratio <- NULL
  if (is_logit) {
# For logit, s_eta = mu_eta / w = 1 - w
    s_eta <- 1 - w_clipped
  } else if (is_probit) {
# For probit, s_eta = phi/Phi (Mills ratio)
    log_phi <- stats::dnorm(eta, log = TRUE)
    log_Phi <- stats::pnorm(eta, log.p = TRUE)
    mills_ratio <- exp(log_phi - log_Phi)
    s_eta <- mills_ratio
  } else if (!is.null(family$score_eta)) {
# Fallback to family-provided score
    s_eta <- family$score_eta(eta, 1)
  } else {
# Generic fallback
    s_eta <- mu_eta / w_clipped
  }

# Derivative of score wrt eta
  if (!is.null(d2mu)) {
    if (is_logit) {
# s_eta = 1 - w gives d/deta(s_eta) = -dw/deta = -mu_eta
      ds_eta_deta <- -mu_eta
    } else if (is_probit) {
# s_eta = r = phi/Phi gives r' = -eta * r - r^2
      if (is.null(mills_ratio)) {
        log_phi <- stats::dnorm(eta, log = TRUE)
        log_Phi <- stats::pnorm(eta, log.p = TRUE)
        mills_ratio <- exp(log_phi - log_Phi)
      }
      ds_eta_deta <- -(eta * mills_ratio + mills_ratio^2)
    } else {
      ds_eta_deta <- (d2mu * w_clipped - mu_eta^2) / (w_clipped^2)
    }
  } else {
    ds_eta_deta <- NULL
  }

  list(
    eta = eta,
    w = w,
    w_clipped = w_clipped,
    mu_eta = mu_eta,
    d2mu = d2mu,
    s_eta = s_eta,
    ds_eta_deta = ds_eta_deta
  )
}

#' Compute lambda_W
#'
#' @param C_const numeric scalar: (N_pop / sum(d_resp) - 1)
#' @param W numeric scalar in (0,1)
#'
#' @keywords internal
el_lambda_W <- function(C_const, W) {
  W <- min(max(as.numeric(W)[1], 1e-12), 1 - 1e-12)
  C_const / (1 - W)
}

#' Compute denominator
#'
#' @param lambda_W numeric scalar
#' @param W numeric scalar in (0,1)
#' @param Xc_lambda numeric vector (X_centered \%*\% lambda_x) or 0
#' @param p_i numeric vector of response probabilities
#' @param floor numeric scalar > 0, denominator floor
#' @return list with denom, active, inv, inv_sq
#'
#' @keywords internal
el_denominator <- function(lambda_W, W, Xc_lambda, p_i, floor) {
  denom <- 1 + lambda_W * (p_i - W)
  if (length(Xc_lambda) > 1L || (is.numeric(Xc_lambda) && Xc_lambda[1] != 0)) {
    denom <- denom + as.numeric(Xc_lambda)
  }
  active <- as.numeric(denom > floor)
  denom_guard <- pmax(denom, floor)
  inv <- 1 / denom_guard
  list(denom = denom_guard, active = active, inv = inv, inv_sq = inv * inv)
}

#' Compute probability masses
#'
#' @param weights numeric respondent base weights (d_i)
#' @param denom numeric denominators Di after floor guard
#' @param floor numeric small positive guard
#' @param trim_cap numeric cap (>0) or Inf
#' @return list with mass_untrim, mass_trimmed, prob_mass, trimmed_fraction
#'
#' @keywords internal
el_masses <- function(weights, denom, floor, trim_cap) {
  mass_untrim <- as.numeric(weights) / as.numeric(denom)
  nn <- enforce_nonneg_weights(mass_untrim)
  if (!nn$ok) stop(nn$message, call. = FALSE)
  mass_untrim <- nn$weights
  trim_res <- trim_weights(mass_untrim, cap = trim_cap)
  mass_trim <- trim_res$weights
  total <- sum(mass_trim)
  prob_mass <- if (total > 0) mass_trim / total else rep(NA_real_, length(mass_trim))
  list(mass_untrim = mass_untrim,
       mass_trimmed = mass_trim,
       prob_mass = prob_mass,
       trimmed_fraction = trim_res$trimmed_fraction)
}

#' Compute the mean
#' @keywords internal
el_mean <- function(prob_mass, y) {
  sum(prob_mass * as.numeric(y))
}
