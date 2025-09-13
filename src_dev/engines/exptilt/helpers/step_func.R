#' @exportS3Method NULL
step_func <- function(model,theta, O_matrix_nieobs) {
  s_values_obs <- s_function(model,1, model$x_1[,model$cols_delta],theta)
  s_values_unobs <- s_function(model,0, model$x_1[,model$cols_delta],theta)

  inv_C <- 1 / as.vector(model$C_matrix_nieobs)
  common_term <- O_matrix_nieobs * model$f_matrix_nieobs * rep(inv_C, each = nrow(O_matrix_nieobs))

  denominator <- rowSums(common_term)

  numerators <- common_term %*% s_values_unobs

  result_nieobs <- colSums(numerators / denominator)

  result_obs <- colSums(s_values_obs)
  test_weights_todoremove <- common_term / denominator
  # browser()
  return(result_nieobs + result_obs)
}


