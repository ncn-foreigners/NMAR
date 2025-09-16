#' @keywords internal
new_nmar_result_exptilt <- function(y_hat, se, weights, coefficients, vcov, converged, class) {


  # if (is.null(data_info$method)) {
  #   data_info$method <- "Exptilt"
  # }


  result <- new_nmar_result(
    y_hat = y_hat,
    se = se,
    weights = weights,
    coefficients = coefficients,
    vcov = vcov,
    converged = converged,
    class = "nmar_result_exptilt"
  )


  return(result)
}
