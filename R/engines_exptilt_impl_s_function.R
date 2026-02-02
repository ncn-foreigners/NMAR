#' @exportS3Method NULL
s_function.nmar_exptilt <- function(model, delta, x, theta = model$theta) {

  theta_numeric <- as.numeric(theta)

  if (delta == 1) {
# Observed case: no expansion needed
    X_full <- cbind(1, x)
    X_full <- apply(X_full, 2, as.numeric)

    eta_full <- as.vector(X_full %*% theta_numeric)

    pi_val_full <- model$family$linkinv(eta_full)
    pi_deriv_full <- model$family$mu.eta(eta_full)
    pi_val_safe <- pmax(pi_val_full, .Machine$double.eps)

    result_matrix <- (pi_deriv_full / pi_val_safe) * X_full
    return(result_matrix)

  } else if (delta == 0) {
# Unobserved case: use outer product to avoid expansion
    x_mat <- as.matrix(x)
    y_vec <- as.vector(model$y_1)
    n_x0 <- nrow(x_mat)
    n_y1 <- length(y_vec)
    p <- length(theta_numeric)
    p_delta <- ncol(x_mat)

# theta = [intercept, theta_delta_1, ..., theta_delta_p, theta_y]
    theta_intercept <- theta_numeric[1]
    theta_delta <- if (p_delta > 0) theta_numeric[2:(1 + p_delta)] else numeric(0)
    theta_y <- theta_numeric[p]

# Compute eta[i,j]
    if (p_delta > 0) {
      x_contrib <- theta_intercept + as.vector(x_mat %*% theta_delta)
    } else {
      x_contrib <- rep(theta_intercept, n_x0)
    }
    eta_matrix <- outer(x_contrib, theta_y * y_vec, "+") # n_x0 × n_y1

# Apply link functions
    pi_val_matrix <- model$family$linkinv(eta_matrix)
    pi_deriv_matrix <- model$family$mu.eta(eta_matrix)
    one_minus_pi_safe <- pmax(1 - pi_val_matrix, .Machine$double.eps)

# Score factor: -pi_deriv / (1 - pi)
    score_factor <- -pi_deriv_matrix / one_minus_pi_safe # n_x0 × n_y1

# expanded output: (n_x0 * n_y1) × p matrix
    n_expanded <- n_x0 * n_y1
    result_matrix <- matrix(0, nrow = n_expanded, ncol = p)

# Parameter 1: Intercept - just the score_factor vectorized
    result_matrix[, 1] <- as.vector(t(score_factor))

# Parameters 2 to (1+p_delta): Delta coefficients
    if (p_delta > 0) {
      for (k in 1:p_delta) {
        score_k <- score_factor * x_mat[, k]
        result_matrix[, 1 + k] <- as.vector(t(score_k))
      }
    }

# Y coefficient
    score_y <- sweep(score_factor, 2, y_vec, "*")
    result_matrix[, p] <- as.vector(t(score_y))

    return(result_matrix)
  } else {
    stop("delta must be 0 or 1 for vectorized calculation.")
  }
}
