#' Stabilized M2(gamma)/M1(gamma) ratio estimator from residuals (optionally weighted)
#'
#' Computes `sum(w * eps_hat * exp(gamma * eps_hat)) / sum(w * exp(gamma * eps_hat))`
#' using a max-shift to avoid overflow/underflow. When `weights` is `NULL`, equal
#' weights are used. When `weights` is supplied, the stabilization shift is based
#' on `log(weights) + gamma * eps_hat` so that zero-weight observations do not
#' affect the reference value.
#'
#' @keywords internal
#' @noRd
il_m2_over_m1_ratio_weighted <- function(eps_hat, gamma, weights = NULL) {
  eps_hat <- as.numeric(eps_hat)

  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || !is.finite(gamma)) {
    stop("`gamma` must be a finite numeric scalar.", call. = FALSE)
  }

  if (length(eps_hat) == 0L) {
    stop("`eps_hat` must have length >= 1.", call. = FALSE)
  }
  if (anyNA(eps_hat) || any(!is.finite(eps_hat))) {
    stop("`eps_hat` must be finite (no NA/Inf).", call. = FALSE)
  }

  z <- gamma * eps_hat
  if (anyNA(z)) stop("Internal error: z has NA values.", call. = FALSE)

  if (is.null(weights)) {
    a <- max(z)
    if (!is.finite(a)) {
      stop(
        "Numerical overflow: max(gamma * eps_hat) is not finite. ",
        "This should not happen with finite inputs; check data scaling.",
        call. = FALSE
      )
    }

    wexp <- exp(z - a)
    denom <- sum(wexp)
    if (!is.finite(denom) || denom <= 0) {
      stop(
        "Numerical error computing stabilized weights for M2/M1 ratio. ",
        "Try standardizing covariates or check for extreme gamma.",
        call. = FALSE
      )
    }
    return(sum(eps_hat * wexp) / denom)
  }

  weights <- as.numeric(weights)
  if (length(weights) != length(eps_hat)) {
    stop("`weights` must have the same length as `eps_hat`.", call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite and nonnegative (no NA/Inf, no negatives).", call. = FALSE)
  }
  sum_w <- sum(weights)
  if (!is.finite(sum_w) || sum_w <= 0) {
    stop("`weights` must sum to a positive number.", call. = FALSE)
  }

# Zero weights must not affect the log-sum-exp shift
# Shifting on log(weights) avoids underflow when weight mass is concentrated
  pos <- weights > 0
  if (!any(pos)) stop("`weights` must have at least one positive entry.", call. = FALSE)
  lw <- log(weights[pos]) + z[pos]
  a <- max(lw)
  if (!is.finite(a)) {
    stop(
      "Numerical overflow: max(log(weights) + gamma * eps_hat) is not finite. ",
      "This should not happen with finite inputs; check data scaling.",
      call. = FALSE
    )
  }

  wexp <- exp(lw - a)
  denom <- sum(wexp)
  if (!is.finite(denom) || denom <= 0) {
    stop(
      "Numerical error computing stabilized weights for M2/M1 ratio. ",
      "Try standardizing covariates or check for extreme gamma.",
      call. = FALSE
    )
  }

  sum(eps_hat[pos] * wexp) / denom
}

#' Stabilized log M1(gamma) estimator from residuals (optionally weighted)
#'
#' Computes `log( sum(w * exp(gamma * eps_hat)) / sum(w) )` using a max-shift to
#' avoid overflow/underflow. When `weights` is `NULL`, equal weights are used.
#' When `weights` is supplied, the stabilization shift is based on
#' `log(weights) + gamma * eps_hat` so that zero-weight observations do not affect
#' the reference value.
#'
#' @keywords internal
#' @noRd
il_log_m1_hat_weighted <- function(eps_hat, gamma, weights = NULL) {
  eps_hat <- as.numeric(eps_hat)

  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || !is.finite(gamma)) {
    stop("`gamma` must be a finite numeric scalar.", call. = FALSE)
  }

  if (length(eps_hat) == 0L) {
    stop("`eps_hat` must have length >= 1.", call. = FALSE)
  }
  if (anyNA(eps_hat) || any(!is.finite(eps_hat))) {
    stop("`eps_hat` must be finite (no NA/Inf).", call. = FALSE)
  }

  z <- gamma * eps_hat
  if (anyNA(z)) stop("Internal error: z has NA values.", call. = FALSE)

  if (is.null(weights)) {
    a <- max(z)
    if (!is.finite(a)) {
      stop(
        "Numerical overflow: max(gamma * eps_hat) is not finite. ",
        "This should not happen with finite inputs; check data scaling.",
        call. = FALSE
      )
    }

    wexp <- exp(z - a)
    sw <- sum(wexp)
    if (!is.finite(sw) || sw <= 0) {
      stop(
        "Numerical error computing stabilized weights for log M1(gamma). ",
        "Try standardizing covariates or check for extreme gamma.",
        call. = FALSE
      )
    }
    return(as.numeric(a + log(sw / length(z))))
  }

  weights <- as.numeric(weights)
  if (length(weights) != length(eps_hat)) {
    stop("`weights` must have the same length as `eps_hat`.", call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite and nonnegative (no NA/Inf, no negatives).", call. = FALSE)
  }
  sum_w <- sum(weights)
  if (!is.finite(sum_w) || sum_w <= 0) {
    stop("`weights` must sum to a positive number.", call. = FALSE)
  }
  pos <- weights > 0
  if (!any(pos)) stop("`weights` must have at least one positive entry.", call. = FALSE)

  lw <- log(weights[pos]) + z[pos]
  a <- max(lw)
  if (!is.finite(a)) {
    stop(
      "Numerical overflow: max(log(weights) + gamma * eps_hat) is not finite. ",
      "This should not happen with finite inputs; check data scaling.",
      call. = FALSE
    )
  }
  sw <- sum(exp(lw - a))
  if (!is.finite(sw) || sw <= 0) {
    stop(
      "Numerical error computing stabilized weights for log M1(gamma). ",
      "Try standardizing covariates or check for extreme gamma.",
      call. = FALSE
    )
  }

  as.numeric((a + log(sw)) - log(sum_w))
}

#' Convenience wrappers (unweighted)
#'
#' @keywords internal
#' @noRd
il_m2_over_m1_ratio <- function(eps_hat, gamma) {
  il_m2_over_m1_ratio_weighted(eps_hat = eps_hat, gamma = gamma, weights = NULL)
}

#' @keywords internal
#' @noRd
il_log_m1_hat <- function(eps_hat, gamma) {
  il_log_m1_hat_weighted(eps_hat = eps_hat, gamma = gamma, weights = NULL)
}
