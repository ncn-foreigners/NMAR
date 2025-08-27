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

    # Liczba nierespondentów (wiersze macierzy wag)
    n_unobs <- nrow(model$x_0)
    # Liczba respondentów (kolumny macierzy wag)
    n_obs <- nrow(model$x_1)

    # Wymiary macierzy FI21
    # Wiersze: liczba parametrów modelu odpowiedzi (phi)
    num_phi_params <- length(model$theta)
    # Kolumny: liczba parametrów modelu gęstości (gamma)
    num_gamma_params <- model$density_num_of_coefs

    # 1. Inicjalizacja macierzy FI21 zerami
    FI21 <- matrix(0, nrow = num_phi_params, ncol = num_gamma_params)

    # Pobranie wag - zakładam, że są już policzone i dostępne jako `model$weights`
    # Wymiar `model$weights` to [n_unobs, n_obs]
    # weights <- model$weights # Użyj już obliczonych wag

    # Pętla po każdym NIERESPONDENCIE (i)
    for (i in 1:n_unobs) {
      # Covariates nierespondenta `i`
      x_i_delta <- model$x_0[i, model$cols_delta, drop = FALSE]
      x_i_gamma <- model$x_for_y_unobs[i, , drop = FALSE]

      # --- Obliczenia dla nierespondenta `i` uśredniane po respondentach `j` ---

      # 2. Oblicz macierz s_ij: wartości funkcji score s(phi) dla x_i oraz każdego y_j respondenta
      # Replikujemy x_i, aby dopasować wymiar do wektora y_1
      x_i_delta_rep <- x_i_delta[rep(1, n_obs), , drop = FALSE]

      # s_ij_matrix będzie miała wymiar [n_obs, num_phi_params]
      # Każdy wiersz `j` to s(phi | delta=0, x_i, y_j)
      s_ij_matrix <- s_function(model, delta = 0, x = x_i_delta_rep, theta = model$theta)

      # 3. Oblicz macierz s1_ij: wartości funkcji score s1(gamma) dla x_i oraz każdego y_j respondenta
      # s1_ij_matrix będzie miała wymiar [n_obs, num_gamma_params]
      # Każdy wiersz `j` to gradient log(f1(y_j | x_i, gamma))
      s1_ij_matrix <- t(sapply(1:n_obs, function(j) {
        model$density_fun_gradient(model$y_1[j], x_i_gamma)
      }))

      # Wagi dla bieżącego nierespondenta `i`
      w_i <- weights[i, ]

      # 4. Oblicz s_bar_0i: warunkowa wartość oczekiwana s(phi) dla nierespondenta `i`
      # To jest wektor o długości num_phi_params
      # browser()
      s_bar_0i <- colSums(w_i * s_ij_matrix)

      # 5. Oblicz odchylenia s_ij od ich średniej s_bar_0i
      s_dev_matrix <- s_ij_matrix - matrix(s_bar_0i, nrow = n_obs, ncol = num_phi_params, byrow = TRUE)

      # 6. Oblicz macierz kowariancji dla nierespondenta `i`
      # To jest suma ważonych iloczynów zewnętrznych wektorów odchyleń i wektorów s1
      # Używamy mnożenia macierzy do efektywnego obliczenia tej sumy:
      # t(A) %*% (w * B) oblicza sum_j(w_j * A_j^T * B_j)
      cov_i <- t(s_dev_matrix) %*% (w_i * s1_ij_matrix)

      # 7. Dodaj wynik do sumy całkowitej
      FI21 <- FI21 + cov_i
    }

    # 8. Zgodnie z formułą, na końcu mnożymy przez -1
    return(-FI21)
  }
  FI21 <- calculate_FI21()
  browser()

}

