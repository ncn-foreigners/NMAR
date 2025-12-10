#' Response model families
#' @name nmar_response_families
#' @description Small family objects used by NMAR engines to compute response
#'   probabilities and derivatives in a link-agnostic way (logit/probit).
#' @details A family object is a list with components:
#'   - `name`: character identifier (e.g., "logit", "probit").
#'   - `linkinv(eta)`: mean function (Bernoulli), returns `p = P(R=1 | eta)`.
#'   - `mu.eta(eta)`: derivative `dp/deta`.
#'   - `d2mu.deta2(eta)`: second derivative `d^2p/deta^2` (used by analytic Jacobian).
#'   - `score_eta(eta, delta)`: score of the Bernoulli log-likelihood w.r.t. `eta`.
#'     For respondents (delta=1), this reduces to `mu.eta(eta)/linkinv(eta)` and is
#'     used in the EL estimating equations. Engines may ignore `delta`.
#'   These functions allow engines to switch between links without changing the
#'   estimating system. See Qin, Leung and Shao (2002) for the EL system.
#' @section Stability:
#'   - Logit: numerically well behaved for a wide range of eta; clip `p` away
#'     from 0/1 when forming ratios.
#'   - Probit: compute `phi/Phi` via a stable log-ratio in tails to avoid 0/0.
#' @name response_families
#' @keywords internal
#' @noRd
NULL

#' Logit family functions (link and derivatives)
#' @keywords internal
logit_family <- function() {
  list(
    name = "logit",
    linkinv = function(eta) stats::plogis(eta),
    mu.eta = function(eta) {
      p <- stats::plogis(eta)
      p * (1 - p)
    },
    d2mu.deta2 = function(eta) {
      p <- stats::plogis(eta)
      p * (1 - p) * (1 - 2 * p)
    },
    score_eta = function(eta, delta) {
      if (missing(delta)) delta <- 1
      delta <- validate_delta(delta)
      p <- stats::plogis(eta)
      delta - p
    }
  )
}

#' Probit family functions (link and derivatives)
#' @keywords internal
probit_family <- function() {
  list(
    name = "probit",
    linkinv = function(eta) stats::pnorm(eta),
    mu.eta = function(eta) stats::dnorm(eta),
    d2mu.deta2 = function(eta) -eta * stats::dnorm(eta),
    score_eta = function(eta, delta) {
      if (missing(delta)) delta <- 1
      delta <- validate_delta(delta)
      n <- length(eta)
      if (!n) return(numeric(0))

# Use log-domain ratios in tails for stability: phi/Phi and phi/(1-Phi)
      log_phi <- stats::dnorm(eta, log = TRUE)
      log_Phi <- stats::pnorm(eta, log.p = TRUE)
      log_tail <- stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE)

      score <- numeric(n)
      idx1 <- delta == 1
      if (any(idx1)) {
        score[idx1] <- exp(log_phi[idx1] - log_Phi[idx1])
      }
      idx0 <- delta == 0
      if (any(idx0)) {
        score[idx0] <- -exp(log_phi[idx0] - log_tail[idx0])
      }
      idx_other <- !(idx1 | idx0)
      if (any(idx_other)) {
# Fall back to the general expression with guarded p
        p_other <- stats::pnorm(eta[idx_other])
        p_other <- pmin(pmax(p_other, .Machine$double.eps), 1 - .Machine$double.eps)
        phi_other <- exp(log_phi[idx_other])
        score[idx_other] <- phi_other * ((delta[idx_other] - p_other) / (p_other * (1 - p_other)))
      }
      score
    }
  )
}

validate_delta <- function(delta) {
  delta <- as.numeric(delta)
  if (any(!is.finite(delta))) {
    stop("`delta` must be finite.", call. = FALSE)
  }
  if (any(delta < 0 | delta > 1)) {
    stop("`delta` must lie in [0, 1] for Bernoulli response models.", call. = FALSE)
  }
  delta
}
