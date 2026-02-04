#' Enforce nonnegativity of weights
#'
#' Softly enforces nonnegativity of a numeric weight vector. Large negative
#' values are treated as errors, while small negative values
#' are truncated to zero.
#'
#' Values below \code{-tol} are treated as clearly negative. Values in
#' \code{[-tol, 0)} are clipped to zero.
#'
#' @param weights numeric vector of weights.
#' @param tol numeric tolerance below which negative values are treated as
#' numerical noise and clipped to zero.
#' @return A list with components:
#' \describe{
#' \item{\code{ok}}{logical; \code{TRUE} if no clearly negative weights
#' were found.}
#' \item{\code{message}}{character; diagnostic message when \code{ok} is
#' \code{FALSE}, otherwise \code{NULL}.}
#' \item{\code{weights}}{numeric vector of adjusted weights (original if
#' \code{ok} is \code{FALSE}, otherwise with small negatives clipped to
#' zero).}
#' }
#'
#' @keywords internal
enforce_nonneg_weights <- function(weights, tol = 1e-8) {
  if (!is.numeric(weights)) stop("`weights` must be numeric.", call. = FALSE)
  if (any(!is.finite(weights))) stop("`weights` must be finite (no NA/NaN/Inf).", call. = FALSE)
  min_w <- suppressWarnings(min(weights))
  if (is.finite(min_w) && min_w < -tol) {
    return(list(
      ok = FALSE,
      message = sprintf("Negative weights detected beyond tolerance (min = %.6f)", min_w),
      weights = weights
    ))
  }
  w <- as.numeric(weights)
  w[w < 0] <- 0
  list(ok = TRUE, message = NULL, weights = w)
}

#' Trim weights by capping and proportional redistribution
#'
#' Applies a cap to a nonnegative weight vector and, when feasible, redistributes
#' excess mass across the remaining positive entries so that the total sum is
#' preserved. When the requested cap is too tight to preserve the total mass,
#' all positive entries are set to the cap and the total sum decreases.
#'
#' Zero weights remain zero. Only entries that are positive after nonnegativity
#' enforcement can absorb redistributed mass.
#'
#' Internally, a simple water-filling style algorithm is used on the positive
#' weights: the largest weights are successively saturated at the cap and the
#' remaining weights are rescaled by a common factor chosen to maintain the
#' total sum.
#'
#' @param weights numeric vector of weights.
#' @param cap positive numeric scalar; maximum allowed weight, or \code{Inf} to disable trimming.
#' @param tol numeric tolerance used when testing whether a rescaling step respects the cap.
#' @param warn_tol numeric tolerance used when testing whether the total sum has been preserved.
#' @return A list with components:
#' \describe{
#' \item{\code{weights}}{numeric vector of trimmed weights.}
#' \item{\code{trimmed_fraction}}{fraction of entries at or very close to
#' the cap (within \code{tol}).}
#' \item{\code{preserved_sum}}{logical; \code{TRUE} if the total sum of
#' weights is preserved to within \code{warn_tol}.}
#' \item{\code{total_before}}{numeric; sum of the original weights.}
#' \item{\code{total_after}}{numeric; sum of the trimmed weights.}
#' }
#'
#' @keywords internal
trim_weights <- function(weights, cap, tol = 1e-12, warn_tol = 1e-8) {
  if (!is.numeric(weights)) stop("`weights` must be numeric.", call. = FALSE)
  if (length(cap) != 1L || !is.numeric(cap) || is.na(cap) || cap <= 0)
    stop("`cap` must be a single positive number or +Inf.", call. = FALSE)
  if (any(!is.finite(weights)))
    stop("`weights` must be finite (no NA/NaN/Inf).", call. = FALSE)

  n <- length(weights)
  if (n == 0L) {
    return(list(
      weights = weights,
      trimmed_fraction = 0,
      preserved_sum = TRUE,
      total_before = 0,
      total_after = 0
    ))
  }

  nn <- enforce_nonneg_weights(weights)
  if (!nn$ok) stop(nn$message, call. = FALSE)
  w <- nn$weights

  if (is.infinite(cap)) {
    tot <- sum(w)
    return(list(
      weights = w,
      trimmed_fraction = 0,
      preserved_sum = TRUE,
      total_before = tot,
      total_after = tot
    ))
  }

  total <- sum(w)
  pos_idx <- which(w > 0)
  mpos <- length(pos_idx)

  if (mpos == 0L) {
    return(list(
      weights = w,
      trimmed_fraction = 0,
      preserved_sum = (abs(total) <= warn_tol),
      total_before = total,
      total_after = sum(w)
    ))
  }

# Only entries that are originally positive can absorb mass.
# If total mass exceeds the maximum that can be assigned to these positions
# under the cap, the total cannot be preserved.
  if (total > mpos * cap + warn_tol) {
    out <- numeric(n)
    out[pos_idx] <- cap
    warning(sprintf(
      "Cannot preserve total mass: sum(weights)=%.6f exceeds %d * cap=%.6f.",
      total, mpos, mpos * cap
    ), call. = FALSE)
    return(list(
      weights = out,
      trimmed_fraction = mean(out >= cap - tol),
      preserved_sum = FALSE,
      total_before = total,
      total_after = sum(out)
    ))
  }

# Water-filling on positive subset: search over the number of the
# capped entries and a common scaling factor for the remaining ones so
# that the total sum is preserved and no weight exceeds the cap
  v <- w[pos_idx]
# Threshold t_i is the scaling factor at which weight i hits the cap
  t <- cap / v
  ord <- order(t) # increasing thresholds
  v <- v[ord]; t <- t[ord]
  m <- length(v)
  suf <- rev(cumsum(rev(v)))

  s <- NA_real_
  k_cap <- 0L
  for (k in 0:(m - 1L)) {
    denom <- suf[k + 1L]
    s_k <- (total - k * cap) / denom
    if (s_k <= t[k + 1L] + tol) {
# In feasible cases, preserving the total requires a scaling factor >= 1
      s <- max(s_k, 1)
      k_cap <- k
      break
    }
  }

  if (is.na(s)) {
    out_pos <- rep(cap, m)
  } else {
    if (k_cap == 0L) {
      out_pos <- c(s * v)
    } else if (k_cap < m) {
      out_pos <- c(rep(cap, k_cap), s * v[(k_cap + 1L):m])
    } else {
      out_pos <- rep(cap, m)
    }
  }

  result <- numeric(n)
  result[pos_idx[ord]] <- pmin(out_pos, cap)
  total_after <- sum(result)
  preserved <- abs(total_after - total) <= warn_tol

  list(
    weights = result,
    trimmed_fraction = mean(result >= cap - tol),
    preserved_sum = preserved,
    total_before = total,
    total_after = total_after
  )
}
