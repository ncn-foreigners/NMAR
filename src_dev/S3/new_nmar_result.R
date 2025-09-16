#' Construct Result Object
#' @keywords internal
new_nmar_result <- function(y_hat, se, weights, coefficients, vcov, converged, class) {

  result <- list(
    y_hat = y_hat,
    se = se,
    weights = weights,
    coefficients = coefficients,
    vcov = vcov,
    converged = converged
  )

  structure(result, class = c(class, "nmar_result"))
}
