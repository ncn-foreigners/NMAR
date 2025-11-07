#' @exportS3Method NULL
estim_mean.nmar_exptilt <- function(model) {
# Build X_delta from declared delta columns intersected with current design
  col_delta <- intersect(model$cols_delta, colnames(model$x_1))
  x_mat <- as.matrix(model$x_1[, col_delta, drop = FALSE])
  x_aug <- cbind(1, x_mat, model$x_1[, model$col_y, drop = FALSE])
  x_aug <- apply(x_aug, 2, as.numeric) # Convert each column to numeric

# Align theta to the augmented design order when names are available
  aug_names <- c("(Intercept)", col_delta, model$col_y)
  theta_raw <- model$theta
  theta <- as.numeric(theta_raw)
  if (!is.null(names(theta_raw))) {
    theta_named <- as.numeric(theta_raw[aug_names])
    if (all(!is.na(theta_named)) && length(theta_named) == ncol(x_aug)) {
      theta <- theta_named
    } else {
# Fallback: keep raw ordering but coerce length to match design
      if (length(theta) != ncol(x_aug)) {
        if (length(theta) > ncol(x_aug)) theta <- theta[seq_len(ncol(x_aug))] else theta <- c(theta, rep(0, ncol(x_aug) - length(theta)))
      }
    }
  } else {
    if (length(theta) != ncol(x_aug)) {
      if (length(theta) > ncol(x_aug)) theta <- theta[seq_len(ncol(x_aug))] else theta <- c(theta, rep(0, ncol(x_aug) - length(theta)))
    }
  }


# browser()
  eta <- as.vector(x_aug %*% theta)

  probabilities <- model$family$linkinv(eta)
  probabilities <- pmax(probabilities, .Machine$double.eps)

  numerator <- sum(model$y_1 * model$design_weights / probabilities)
  denominator <- sum(model$design_weights / probabilities)

  return(numerator / denominator)
}
