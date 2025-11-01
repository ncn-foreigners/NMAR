#' Enforce (near-)nonnegativity of weights
#' @keywords internal
enforce_nonneg_weights <- function(weights, tol = 1e-8) {
  if (!is.numeric(weights)) stop("`weights` must be numeric.", call. = FALSE)
  if (any(!is.finite(weights))) stop("`weights` must be finite (no NA/NaN/Inf).", call. = FALSE)
  min_w <- suppressWarnings(min(weights))
  if (is.finite(min_w) && min_w < -tol) {
    return(list(
      ok = FALSE,
      message = sprintf("Negative EL weights produced (min = %.6f)", min_w),
      weights = weights
    ))
  }
  w <- as.numeric(weights)
  w[w < 0] <- 0
  list(ok = TRUE, message = NULL, weights = w)
}

#' Trim weights by capping and proportional redistribution (water-filling)
#' @keywords internal
trim_weights <- function(weights, cap, tol = 1e-12, warn_tol = 1e-8) {
# Validation
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

# Nonnegativity (soft enforcement)
  nn <- enforce_nonneg_weights(weights)
  if (!nn$ok) stop(nn$message, call. = FALSE)
  w <- nn$weights

# No cap -> identity
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

# Feasibility: only original-positive entries can absorb mass
  if (total > mpos * cap + warn_tol) {
    out <- numeric(n)
    out[pos_idx] <- cap
    warning(sprintf(
      "EL mass trimming: cannot preserve total mass: sum(weights)=%.6f exceeds %d * cap=%.6f.",
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

# Water-filling on positive subset
  v <- w[pos_idx]
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
