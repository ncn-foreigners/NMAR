calculate_fisher_information <- function(y_1,x_1,density_num_of_coefs,x_for_y_obs,density_fun_hess) {
  n_obs <- nrow(x_1)
  n_params <- density_num_of_coefs

  # Inicjalizacja macierzy informacji Fishera
  FI <- matrix(0, nrow = n_params, ncol = n_params)
  # param_names <- names(bbmle::coef(model$density_result$model))
  # rownames(FI) <- colnames(FI) <- param_names

  # Dla każdej obserwacji oblicz i sumuj oczekiwany hesjan
  for (i in 1:n_obs) {
    x_row <- x_for_y_obs[i, , drop = FALSE]
    y_value <- y_1[i]

    # Oblicz hesjan dla tej obserwacji
    hessian <- density_fun_hess(y_value, x_row)

    if (!is.null(hessian)) {
      # Informacja Fishera to minus oczekiwana wartość hesjanu
      FI <- FI - hessian
    }
  }

  return(FI)
}
