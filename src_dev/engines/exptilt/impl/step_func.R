#' @exportS3Method NULL

step_func <- function(model, theta, O_matrix_nieobs) {

  n_x1 <- nrow(model$x_1)
  n_x0 <- nrow(model$x_0)
  n_y1 <- length(model$y_1)

  theta_numeric <- as.numeric(theta)
  p <- length(theta_numeric)

  inv_C <- 1 * model$design_weights / as.vector(model$C_matrix_nieobs)
  inv_C[!is.finite(inv_C)] <- 0
  W_numerator_matrix <- O_matrix_nieobs * model$f_matrix_nieobs

  common_term <- sweep(W_numerator_matrix, MARGIN = 2, STATS = inv_C, FUN = "*")
# browser()
  X_obs_for_s_func <- cbind(model$x_1[, model$cols_delta], model$y_1[1:n_x1])

  s_values_obs <- s_function(model, 1, X_obs_for_s_func, theta_numeric)

  s_values_obs <- s_values_obs * model$design_weights[1:n_x1]

# Unobserved contribution
  s_values_unobs_expanded <- s_function(model, 0, model$x_0[, model$cols_delta, drop = FALSE], theta_numeric)

# Optimize: avoid transpose and multiple copies
# common_term is n_x0 × n_y1
# We need it in row-major order matching s_values_unobs_expanded
# s_values_unobs_expanded rows correspond to: all y for x[1], all y for x[2], ...
# which matches t(common_term) vectorized

  common_term_vec <- as.vector(t(common_term)) # length n_x0 * n_y1

# Weight the scores
  s_values_unobs_weighted <- common_term_vec * s_values_unobs_expanded

# Reshape and sum
# Each column k is the score for parameter k at all (i,j) pairs
# Need to sum within each i (across all j), then sum across i

  numerators <- matrix(0, nrow = n_x0, ncol = p)
  for (k in 1:p) {
    score_matrix_k <- matrix(s_values_unobs_weighted[, k], nrow = n_y1, ncol = n_x0)
    numerators[, k] <- colSums(score_matrix_k)
  }

  denominator <- rowSums(common_term)

# cat('Dim of numerators:', dim(numerators), '\n')
# cat('Dim of denominators:', length(denominator), '\n')
# cat('Dim of s_values_obs:', dim(s_values_obs), '\n')
# cat('Dim of s_values_unobs:', dim(numerators), '\n')
# cat('Dim of inv_C:', length(inv_C), '\n')
# cat('Dim of common_term:', dim(common_term), '\n')
# cat('Dim of O_matrix_nieobs:', dim(O_matrix_nieobs), '\n')
# cat('Dim of model$f_matrix_nieobs:', dim(model$f_matrix_nieobs), '\n')

  result_nieobs <- colSums(numerators / denominator)
  result_obs <- colSums(s_values_obs)

  result_total <- result_nieobs + result_obs

  return(result_total)
}
