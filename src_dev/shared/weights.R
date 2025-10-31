#' Weight utilities
#' @keywords internal
trim_weights <- function(weights, cap) {
  if (is.infinite(cap)) {
    return(list(weights = weights, trimmed_fraction = 0))
  }
# Guard against NA/NaN values
  if (any(!is.finite(weights))) {
    stop("Weights contain non-finite values (NA, NaN, or Inf)", call. = FALSE)
  }
  final_weights <- weights
  original_total <- sum(weights)
  eligible_for_redistribution <- rep(TRUE, length(weights))
  while (TRUE) {
    weights_to_cap <- (final_weights > cap) & eligible_for_redistribution
    if (!any(weights_to_cap, na.rm = TRUE)) break
    excess <- sum(final_weights[weights_to_cap] - cap)
    final_weights[weights_to_cap] <- cap
    eligible_for_redistribution[weights_to_cap] <- FALSE
    eligible_weights <- final_weights[eligible_for_redistribution]
    if (sum(eligible_weights) < 1e-8) break
    redistribution <- excess * (eligible_weights / sum(eligible_weights))
    final_weights[eligible_for_redistribution] <- eligible_weights + redistribution
  }
  total_after <- sum(final_weights)
  if (abs(total_after - original_total) > 1e-8) {
    warning(
      sprintf(
        "Trimming reduced total weight from %.6f to %.6f; consider using a larger cap.",
        original_total, total_after
      ),
      call. = FALSE
    )
  }
  list(weights = final_weights, trimmed_fraction = mean(weights > cap))
}

#' @keywords internal
enforce_nonneg_weights <- function(weights, tol = 1e-8) {
  min_w <- suppressWarnings(min(weights, na.rm = TRUE))
  if (is.finite(min_w) && min_w < -tol) {
    return(list(ok = FALSE, message = sprintf("Negative EL weights produced (min = %.6f)", min_w), weights = weights))
  }
  weights[weights < 0] <- 0
  list(ok = TRUE, message = NULL, weights = weights)
}
