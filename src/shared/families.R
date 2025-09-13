#' Response model families
#' @name nmar_response_families
#' @description Small family objects used by NMAR engines to compute response
#'   probabilities and derivatives in a link‑agnostic way (logit/probit).
#' @details A family object is a list with components:
#'   - `name`: character identifier (e.g., "logit", "probit").
#'   - `linkinv(eta)`: mean function (Bernoulli), returns `p = P(R=1 | eta)`.
#'   - `mu.eta(eta)`: derivative `dp/deta`.
#'   - `d2mu.deta2(eta)`: second derivative `d^2p/deta^2` (used by analytic Jacobian).
#'   - `score_eta(eta, delta)`: score of the Bernoulli log‑likelihood w.r.t. `eta`.
#'     For respondents (delta=1), this reduces to `mu.eta(eta)/linkinv(eta)` and is
#'     used in the EL estimating equations. Engines may ignore `delta`.
#'   These functions allow engines to switch between links without changing the
#'   estimating system. See Qin, Leung and Shao (2002) for the EL system.
#' @section Stability:
#'   - Logit: numerically well behaved for a wide range of eta; clip `p` away
#'     from 0/1 when forming ratios.
#'   - Probit: compute `phi/Phi` via a stable log‑ratio in tails to avoid 0/0.
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
      # Score wrt eta for log-likelihood of response probability: d/deta log p(eta)
      # This equals mu.eta / p for any Bernoulli mean model.
      p <- stats::plogis(eta)
      p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
      m <- p * (1 - p)
      m / p
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
      # Score wrt eta for log-likelihood of response probability: d/deta log p(eta)
      # Compute via stable log form for extreme tails.
      # score = phi(eta) / Phi(eta)
      log_phi <- stats::dnorm(eta, log = TRUE)
      log_Phi <- stats::pnorm(eta, log.p = TRUE)
      # When p is extremely close to 1, use (1-Phi) tail to avoid 0/0; fallback to direct ratio.
      s <- exp(log_phi - log_Phi)
      # Guard against under/overflow
      s[!is.finite(s)] <- {
        p <- stats::pnorm(eta)
        p <- pmin(pmax(p, 1e-300), 1 - 1e-300)
        stats::dnorm(eta) / p
      }[!is.finite(s)]
      s
    }
  )
}
