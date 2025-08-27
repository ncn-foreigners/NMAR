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

  #it is OK if values are close to 0
  F11 = calculate_fisher_information(model$y_1,model$x_1,model$density_num_of_coefs,model$x_for_y_obs,model$density_fun_hess)
  # browser()

  calculate_FI21 <- function() {


    n_unobs <- nrow(model$x_0)

    n_obs <- nrow(model$x_1)


    num_phi_params <- length(model$theta)
    num_gamma_params <- model$density_num_of_coefs

    FI21 <- matrix(0, nrow = num_phi_params, ncol = num_gamma_params)

    #TODO: optimize
    for (i in 1:n_unobs) {

      x_i_delta <- model$x_0[i, model$cols_delta, drop = FALSE]
      x_i_gamma <- model$x_for_y_unobs[i, , drop = FALSE]


      x_i_delta_rep <- x_i_delta[rep(1, n_obs), , drop = FALSE]

      s_ij_matrix <- s_function(model, delta = 0, x = x_i_delta_rep, theta = model$theta)


      s1_ij_matrix <- t(sapply(1:n_obs, function(j) {
        model$density_fun_gradient(model$y_1[j], x_i_gamma)
      }))


      w_i <- weights[i, ]


      # browser()
      s_bar_0i <- colSums(w_i * s_ij_matrix)


      s_dev_matrix <- s_ij_matrix - matrix(s_bar_0i, nrow = n_obs, ncol = num_phi_params, byrow = TRUE)


      # t(A) %*% (w * B) oblicza sum_j(w_j * A_j^T * B_j)
      cov_i <- t(s_dev_matrix) %*% (w_i * s1_ij_matrix)


      FI21 <- FI21 + cov_i
    }

    return(-FI21)
  }
  FI21 <- calculate_FI21()
  browser()

}

