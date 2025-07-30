
s_function <- function(model,delta,x,theta=model$theta) {
  # stopifnot(
  #   !any(is.na(delta)),
  #   delta %in% c(0, 1),
  #   length(theta) == ncol(x) + 1 + 1 # +1 for intercept, +1 for y
  # )
  x_mat <- as.matrix(x)
  y_vec <- as.vector(model$y_1)

  pi_vals <- pi_func(model, x_mat,  func = "reg",theta=theta)
  pi_deriv <- pi_func(model,x_mat, func = "deriv",theta=theta)


  # numerator <- pi_deriv *(delta - pi_vals)
  # denominator <- pi_vals * (1 - pi_vals)
  numerator <- NULL
  denominator <- NULL

  #optimization
  if(delta==1){
    numerator <- pi_deriv
    denominator <- pi_vals
  }
  else if(delta==0){
    numerator <- -pi_deriv
    denominator <- 1 - pi_vals
  }



  result <- numerator / denominator
  # result <- sweep(pi_deriv, MARGIN = 1, STATS = 1 / pi_vals, FUN = "*") #CZAT
  result[is.nan(result)] <- 0  # Obsługa przypadków gdy pi = 0 lub 1
  result
}
