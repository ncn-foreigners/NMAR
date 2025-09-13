# estim_mean.nmar_exptilt <- function(model){
#   numerator<-sum(model$y_1/pi_func(model, model$x_1[,model$cols_delta], model$x_1[,model$col_y], func = "reg"))
#   denominator<-sum(1/pi_func(model, model$x_1[,model$cols_delta], model$x_1[,model$col_y], func = "reg"))
#
#   return(sum(numerator/denominator))
# }
#
#' @exportS3Method NULL
estim_mean.nmar_exptilt <- function(model) {
  family <- if (model$prob_model_type == "logit") {
    logit_family()
  } else if (model$prob_model_type == "probit") {
    probit_family()
  }

  x_mat <- as.matrix(model$x_1[, model$cols_delta])
  x_aug <- cbind(1, x_mat, model$x_1[, model$col_y])

  eta <- as.vector(x_aug %*% model$theta)

  probabilities <- family$linkinv(eta)

  numerator <- sum(model$y_1 / probabilities)
  denominator <- sum(1 / probabilities)

  return(numerator / denominator)
}
