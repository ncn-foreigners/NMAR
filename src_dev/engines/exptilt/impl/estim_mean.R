#' @exportS3Method NULL
estim_mean.nmar_exptilt <- function(model) {


  x_mat <- as.matrix(model$x_1[, model$cols_delta])
  x_aug <- cbind(1, x_mat, model$x_1[, model$col_y])
  x_aug <- apply(x_aug, 2, as.numeric) # Convert each column to numeric

  theta <- unname(model$theta) # Remove names


# browser()
  eta <- as.vector(x_aug %*% theta)

  probabilities <- model$family$linkinv(eta)

  numerator <- sum(model$y_1 * model$design_weights / probabilities)
  denominator <- sum(model$design_weights / probabilities)

  return(numerator / denominator)
}
