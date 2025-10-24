#' @exportS3Method NULL
s_function.nmar_exptilt <- function(model, delta, x, theta = model$theta) {

  theta_numeric <- as.numeric(theta)
  SAFE_THRESHOLD <- 1e-3

  if (delta == 1) {
# Observed case: no expansion needed
    X_full <- cbind(1, x)
    X_full <- apply(X_full, 2, as.numeric)

    eta_full <- as.vector(X_full %*% theta_numeric)

    pi_val_full <- model$family$linkinv(eta_full)
    pi_deriv_full <- model$family$mu.eta(eta_full)
    pi_val_safe <- pmin(pmax(pi_val_full, SAFE_THRESHOLD), 1 - SAFE_THRESHOLD)

    result_matrix <- (pi_deriv_full / pi_val_safe) * X_full
    return(result_matrix)

  } else if (delta == 0) {
# Unobserved case: use outer product to avoid expansion
# Goal: compute score for each (x[i], y[j]) pair
# score[row, param] where row = (i-1)*n_y1 + j (expanded index)

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

# Compute eta[i,j] using outer product
    if (p_delta > 0) {
      x_contrib <- theta_intercept + as.vector(x_mat %*% theta_delta)
    } else {
      x_contrib <- rep(theta_intercept, n_x0)
    }
    eta_matrix <- outer(x_contrib, theta_y * y_vec, "+") # n_x0 × n_y1

# Apply link functions
    pi_val_matrix <- model$family$linkinv(eta_matrix)
    pi_deriv_matrix <- model$family$mu.eta(eta_matrix)
    pi_val_safe <- pmin(pmax(pi_val_matrix, SAFE_THRESHOLD), 1 - SAFE_THRESHOLD)

# Score factor: -pi_deriv / (1 - pi)
    score_factor <- -pi_deriv_matrix / (1 - pi_val_safe) # n_x0 × n_y1

# Now create expanded output: (n_x0 * n_y1) × p matrix
# Row ordering: all j for i=1, then all j for i=2, etc.
    n_expanded <- n_x0 * n_y1
    result_matrix <- matrix(0, nrow = n_expanded, ncol = p)

# Parameter 1: Intercept - just the score_factor vectorized
    result_matrix[, 1] <- as.vector(t(score_factor))

# Parameters 2 to (1+p_delta): Delta coefficients
    if (p_delta > 0) {
      for (k in 1:p_delta) {
# score[i,j,k] = score_factor[i,j] * x[i,k]
# Multiply each row i of score_factor by x[i,k]
        score_k <- score_factor * x_mat[, k] # R recycles x_mat[,k] across columns
        result_matrix[, 1 + k] <- as.vector(t(score_k))
      }
    }

# Last parameter: Y coefficient
# score[i,j,p] = score_factor[i,j] * y[j]
# Multiply each column j of score_factor by y[j]
    score_y <- sweep(score_factor, 2, y_vec, "*")
    result_matrix[, p] <- as.vector(t(score_y))

    return(result_matrix)
  } else {
    stop("delta must be 0 or 1 for vectorized calculation.")
  }
}
