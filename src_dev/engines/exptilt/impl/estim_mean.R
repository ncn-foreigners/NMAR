#' @exportS3Method NULL
estim_mean.nmar_exptilt <- function(model) {


  x_mat <- as.matrix(model$data_1[, model$cols_delta])
  x_aug <- cbind(1, x_mat, model$data_1[, model$col_y])
  x_aug <- apply(x_aug, 2, as.numeric)

  theta <- unname(model$theta) # Remove names


# browser()
  eta <- as.vector(x_aug %*% theta)

  probabilities <- model$family$linkinv(eta)
  probabilities <- pmax(probabilities, .Machine$double.eps)

  numerator <- sum(model$y_1
                   * model$design_weights[model$respondent_mask]
                   / probabilities)
  denominator <- sum(
    1
    * model$design_weights[model$respondent_mask]
    / probabilities)

  test <- numerator / denominator
# browser()
  return(numerator / denominator)
}
