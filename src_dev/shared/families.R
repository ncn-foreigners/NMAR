#' Response model families
#'
#' Small family objects used by NMAR engines to compute response probabilities
#' and derivatives for Bernoulli response models with different links (logit,
#' probit).
#'
#' A family object is a list with components:
#' \describe{
#'   \item{\code{name}}{Character identifier, e.g. \code{"logit"}.}
#'   \item{\code{linkinv(eta)}}{Mean function returning \code{p = P(R = 1 | eta)}.}
#'   \item{\code{mu.eta(eta)}}{First derivative \code{dp/deta}.}
#'   \item{\code{d2mu.deta2(eta)}}{Second derivative \code{d^2p/deta^2} (optional).}
#'   \item{\code{score_eta(eta, delta)}}{Score with respect to \code{eta} for a
#'     Bernoulli log-likelihood. For respondents (\code{delta = 1}), this reduces
#'     to \code{mu.eta(eta) / linkinv(eta)}.}
#' }
#'
#' These functions let engines switch links without changing their estimating
#' equations. The EL engine follows Qin, Leung, and Shao (2002).
#'
#' @section Numerical stability:
#'   Probit score calculations use log-domain ratios in the tails to compute the
#'   Mills ratio \code{phi/Phi} without underflow. Engines should clip
#'   probabilities away from 0 and 1 when forming ratios.
#'
#' @name nmar_response_families
#' @keywords internal
#' @noRd
NULL

#' Construct a logit response family bundle
#' @return A list with components \code{name}, \code{linkinv}, \code{mu.eta},
#'   \code{d2mu.deta2}, and \code{score_eta}.
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

#' Construct a probit response family bundle
#' @return A list with components \code{name}, \code{linkinv}, \code{mu.eta},
#'   \code{d2mu.deta2}, and \code{score_eta}.
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

# Use log-domain ratios in tails for stability: phi/Phi and phi/(1-Phi).
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

# Internal helper: validate Bernoulli response indicators.
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
