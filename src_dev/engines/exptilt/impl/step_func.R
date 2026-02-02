#' @exportS3Method NULL

step_func <- function(model, theta, O_matrix_nieobs) {
  n_x1 <- nrow(model$data_1)
  n_x0 <- nrow(model$data_0)
  n_y1 <- length(model$y_1)

  theta_numeric <- as.numeric(theta)
  p <- length(theta_numeric)

  C_safe <- pmax(abs(model$C_matrix_nieobs), .Machine$double.eps)
  inv_C <- 1 / as.vector(C_safe)
  W_numerator_matrix <- O_matrix_nieobs * model$f_matrix_nieobs * model$design_weights[!model$respondent_mask]

  common_term <- sweep(W_numerator_matrix, MARGIN = 2, STATS = inv_C, FUN = "*")
# browser()
  X_obs_for_s_func <- cbind(model$data_1[, model$cols_delta], model$y_1[1:n_x1])

  s_values_obs <- s_function(model, 1, X_obs_for_s_func, theta_numeric)

  s_values_obs <- s_values_obs * model$design_weights[model$respondent_mask]

# Unobserved contribution
  s_values_unobs_expanded <- s_function(model, 0, model$data_0[, model$cols_delta, drop = FALSE], theta_numeric)


  common_term_vec <- as.vector(t(common_term)) # length n_x0 * n_y1

  s_values_unobs_weighted <- common_term_vec * s_values_unobs_expanded


  numerators <- matrix(0, nrow = n_x0, ncol = p)
  for (k in 1:p) {
    score_matrix_k <- matrix(s_values_unobs_weighted[, k], nrow = n_y1, ncol = n_x0)
    numerators[, k] <- colSums(score_matrix_k)
  }
  numerators <- numerators * model$design_weights[!model$respondent_mask]

  denominator <- rowSums(common_term)
  denominator_safe <- pmax(denominator, .Machine$double.eps)

  result_nieobs <- colSums(numerators / denominator_safe)
  result_obs <- colSums(s_values_obs)

  result_total <- result_nieobs + result_obs
# browser()
# browser()
  return(result_total)
}
