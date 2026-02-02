#' Construct logit response family
#' @return A list with components \code{name}, \code{linkinv}, \code{mu.eta},
#' \code{d2mu.deta2}, and \code{score_eta}.
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

#' Construct probit response family
#' @return A list with components \code{name}, \code{linkinv}, \code{mu.eta},
#' \code{d2mu.deta2}, and \code{score_eta}.
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

# log-domain ratios in tails for stability
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

#' @keywords internal
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
