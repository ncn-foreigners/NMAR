estim_var.nmar_exptilt <- function(model){
  s_values_unobs <- s_function(model,0, model$x_1[,model$cols_delta],theta)

  inv_C <- 1 / as.vector(model$C_matrix_nieobs)
  common_term <- model$O_matrix_nieobs * model$f_matrix_nieobs * rep(inv_C, each = nrow(model$O_matrix_nieobs))

  denominator <- rowSums(common_term)


  weights <- common_term / denominator #TODO test if dim is matching
  p = pi_func(model,model$x_1[,model$cols_delta],'reg',model$theta)

  #it is ok if Var is close to 0
  cat("Mean of p:", mean(p), "Variance of p:", var(p), "\n")


  #density num of coefs refers to density f.e beta, intercept, sigma
  # S1=matrix(0, nrow=nrow(model$x_1), ncol=model$density_num_of_coefs)
  S1 <- t(sapply(1:nrow(model$x_1), function(i) {
    model$density_fun_gradient(
      model$y_1[i],
      model$x_for_y_obs[i, , drop = FALSE]
    )
  }))

  F11 = calculate_fisher_information(model$y_1,model$x_1,model$density_num_of_coefs,model$x_for_y_obs,model$density_fun_hess)
  browser()
}

