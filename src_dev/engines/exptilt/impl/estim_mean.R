
#' @exportS3Method NULL
estim_mean.nmar_exptilt <- function(model) {


  x_mat <- as.matrix(model$x_1[, model$cols_delta])
  x_aug <- cbind(1, x_mat, model$x_1[, model$col_y])

  eta <- as.vector(x_aug %*% model$theta)

  probabilities <- model$family$linkinv(eta)

  numerator <- sum(model$y_1 * model$respondent_weights / probabilities)
  denominator <- sum(model$respondent_weights / probabilities)

  return(numerator / denominator)
}
