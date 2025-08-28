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

  z_function <- function(model, x, theta = model$theta) {

    pi_vals <- pi_func(model, x, func = "reg", theta = theta)


    pi_deriv <- pi_func(model, x, func = "deriv", theta = theta)


    denominator <- pi_vals * (1 - pi_vals)


    denominator[denominator < 1e-9] <- 1 # zero case

    result <- pi_deriv / as.vector(denominator)
    result[is.nan(result)] <- 0

    return(result)
  }

  calculate_FI22 <- function(model, weights) {

    n_unobs <- nrow(model$x_0)
    n_obs <- nrow(model$x_1)
    num_phi_params <- length(model$theta)


    FI22 <- matrix(0, nrow = num_phi_params, ncol = num_phi_params)

    #TODO optimize
    for (i in 1:n_unobs) {

      x_i_delta <- model$x_0[i, model$cols_delta, drop = FALSE]
      x_i_delta_rep <- x_i_delta[rep(1, n_obs), , drop = FALSE]


      #TODO verify - this code might be repeated
      s_ij_matrix <- s_function.nmar_exptilt(model, delta = 0, x = x_i_delta_rep, theta = model$theta)
      w_i <- weights[i, ]
      s_bar_0i <- colSums(w_i * s_ij_matrix)


      z_ij_matrix <- z_function(model, x = x_i_delta_rep, theta = model$theta)
      z_bar_0i <- colSums(w_i * z_ij_matrix)

      outer_prod_i <- s_bar_0i %*% t(z_bar_0i)
      FI22 <- FI22 + outer_prod_i
      # browser()
    }

    return(-FI22)
  }
  #TODO - CHECK Below. Fi22 Too big and K seems too low comparing to author
  FI22 <- calculate_FI22(model, weights)
  K <- FI21 %*% solve(F11)
  browser()


  calculate_B <- function(model, esty, FI22) {

    # Pobierz zmienne dla respondentów
    x_obs_delta <- model$x_1[, model$cols_delta, drop = FALSE]
    y_obs <- model$y_1

    # Oblicz p i u_i dla respondentów
    p_obs <- pi_func(model, x_obs_delta, func = "reg", theta = model$theta)
    u_i <- y_obs - esty

    # Oblicz pochodną pi po phi dla respondentów
    pi_deriv_obs <- pi_func(model, x_obs_delta, func = "deriv", theta = model$theta)

    # Oblicz pierwszą część formuły B (suma)
    # W artykule autora `B1`
    sum_term <- t(u_i / p_obs) %*% pi_deriv_obs

    # B = suma %*% odwrotność(FI22)
    B <- sum_term %*% solve(FI22)

    return(B)
  }

  calculate_S2 <- function(model, weights) {

    n_total <- nrow(model$x)
    num_phi_params <- length(model$theta)

    S2 <- matrix(0, nrow = n_total, ncol = num_phi_params)

    # Identyfikatory respondentów i nierespondentów w pełnym zbiorze danych
    respondent_indices <- which(!is.na(model$x[, model$col_y]))
    non_respondent_indices <- which(is.na(model$x[, model$col_y]))

    # --- 1. Obliczenia dla respondentów ---
    x_obs_delta <- model$x_1[, model$cols_delta, drop = FALSE]
    S2[respondent_indices, ] <- s_function.nmar_exptilt(model, delta = 1, x = x_obs_delta, theta = model$theta)

    # --- 2. Obliczenia dla nierespondentów ---
    n_unobs <- nrow(model$x_0)
    n_obs <- nrow(model$x_1)

    # W pętli, ponieważ każda `s_bar_0i` jest specyficzna dla `x_i` nierespondenta
    for (i in 1:n_unobs) {
      x_i_delta <- model$x_0[i, model$cols_delta, drop = FALSE]
      x_i_delta_rep <- x_i_delta[rep(1, n_obs), , drop = FALSE]

      s_ij_matrix <- s_function.nmar_exptilt(model, delta = 0, x = x_i_delta_rep, theta = model$theta)
      w_i <- weights[i, ]
      s_bar_0i <- colSums(w_i * s_ij_matrix)

      # Przypisz do odpowiedniego wiersza w macierzy S2
      S2[non_respondent_indices[i], ] <- s_bar_0i
    }

    return(S2)
  }



  #TODO below 2 lines repeated - refactor
  p_obs <- pi_func(model, model$x_1[, model$cols_delta, drop = FALSE], func = "reg", theta = model$theta)
  esty <- sum(model$y_1 / p_obs) / sum(1 / p_obs)


  B <- calculate_B(model, esty, FI22)


  S2 <- calculate_S2(model, weights)


  n_total <- nrow(model$x)
  eta <- numeric(n_total)


  respondent_indices <- which(!is.na(model$x[, model$col_y]))


  u_i_obs <- model$y_1 - esty


  correction_term_obs <- S2[respondent_indices, ] - S1 %*% t(K)


  eta[respondent_indices] <- (u_i_obs / p_obs) - B %*% t(correction_term_obs)


  non_respondent_indices <- which(is.na(model$x[, model$col_y]))
  correction_term_unobs <- S2[non_respondent_indices, ]

  eta[non_respondent_indices] <- -B %*% t(correction_term_unobs)


  tau <- sum(1 / p_obs)


  var_est <- n_total * var(eta) / (tau^2)

  #TODO below for test only
  # cat("Estimated Mean (esty):", esty, "\n")
  cat("Estimated Variance (var_est):", var_est, "\n")

  # browser()
  return(var_est)
}

