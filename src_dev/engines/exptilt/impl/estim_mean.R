#' @exportS3Method NULL
estim_mean.nmar_exptilt <- function(model) {
  x_mat <- as.matrix(model$x_1[, model$cols_delta])
  x_aug <- cbind(1, x_mat, model$x_1[, model$col_y])

  eta <- as.vector(x_aug %*% model$theta)

  probabilities <- model$family$linkinv(eta)

  if (!is.null(model$design_weights) && !is.null(model$respondent_mask)) {
    resp_w_local <- model$design_weights[model$respondent_mask]
  } else {
    resp_w_local <- rep(1, length(model$y_1))
  }
  numerator <- sum(model$y_1 * resp_w_local / probabilities)
  denominator <- sum(resp_w_local / probabilities)

  return(numerator / denominator)
}
