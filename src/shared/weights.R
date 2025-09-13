#' Weight utilities
#' @keywords internal
trim_weights <- function(weights, cap) {
  if (is.infinite(cap)) {
    return(list(weights = weights, trimmed_fraction = 0))
  }
  final_weights <- weights
  eligible_for_redistribution <- rep(TRUE, length(weights))
  while (TRUE) {
    weights_to_cap <- (final_weights > cap) & eligible_for_redistribution
    if (!any(weights_to_cap)) break
    excess <- sum(final_weights[weights_to_cap] - cap)
    final_weights[weights_to_cap] <- cap
    eligible_for_redistribution[weights_to_cap] <- FALSE
    eligible_weights <- final_weights[eligible_for_redistribution]
    if (sum(eligible_weights) < 1e-8) break
    redistribution <- excess * (eligible_weights / sum(eligible_weights))
    final_weights[eligible_for_redistribution] <- eligible_weights + redistribution
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
