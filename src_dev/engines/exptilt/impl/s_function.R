#' @exportS3Method NULL
s_function.nmar_exptilt <- function(model, delta, x, theta = model$theta) {


  x_mat <- as.matrix(x)

  x_aug <- cbind(1, x_mat, model$y_1[1:nrow(x_mat)])

  eta <- as.vector(x_aug %*% theta)
  pi_vals <- model$family$linkinv(eta)

  pi_deriv <- model$family$mu.eta(eta) * x_aug

  numerator <- NULL
  denominator <- NULL


  if (delta == 1) {
    numerator <- pi_deriv
    denominator <- pi_vals
  } else if (delta == 0) {
    numerator <- -pi_deriv
    denominator <- 1 - pi_vals
  }

  result <- numerator / denominator
  result[is.nan(result)] <- 0
  result
}
