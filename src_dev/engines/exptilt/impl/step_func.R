#' #' @exportS3Method NULL
#'
step_func <- function(model, theta, O_matrix_nieobs) {

  # Get dimensions
  n_x1 <- nrow(model$x_1)
  n_x0 <- nrow(model$x_0)
  n_y1 <- length(model$y_1)
  p <- length(theta)

  # Calculate common_term FIRST
  inv_C <- 1 * model$respondent_weights / as.vector(model$C_matrix_nieobs)
  inv_C[!is.finite(inv_C)] <- 0

  # common_term <- O_matrix_nieobs * model$f_matrix_nieobs * rep(inv_C, each = nrow(O_matrix_nieobs)) * model$respondent_weights
  W_numerator_matrix <- O_matrix_nieobs * model$f_matrix_nieobs

  # Now, correctly scale each column j by the corresponding inv_C[j] value
  # MARGIN = 2 specifies that the operation should be applied column-wise.
  common_term <- sweep(W_numerator_matrix, MARGIN = 2, STATS = inv_C, FUN = "*")

  # Calculate s_values_obs using i,j approach
  s_values_obs <- matrix(0, nrow = n_x1, ncol = p)
  for (i in 1:n_x1) {
    s_val_ij <- s_function(model, 1, model$x_1[, model$cols_delta], theta, i = i, j = i)
    s_values_obs[i, ] <- s_val_ij * model$respondent_weights[i]
  }

  # Calculate s_values_unobs using i,j approach for common_term multiplication
  numerators <- matrix(0, nrow = n_x0, ncol = p)
  for (i in 1:n_x0) {
    for (j in 1:n_y1) {
      s_val_ij <- s_function(model, 0, model$x_0[, model$cols_delta], theta, i = i, j = j)
      # s_val_ij<-c(1,1)
      numerators[i, ] <- numerators[i, ] + common_term[i, j] * s_val_ij
    }
  }

  denominator <- rowSums(common_term)
  # Avoid division by values close to 0
  # denominator[denominator < 0.01] <- 0.01
  cat('Dim of numerators:', dim(numerators), '\n')
  cat('Dim of denominators:', length(denominator), '\n')
  cat('Dim of s_values_obs:', dim(s_values_obs), '\n')
  cat('Dim of s_values_unobs:', dim(numerators), '\n')
  cat('Dim of inv_C:', length(inv_C), '\n')
  cat('Dim of common_term:', dim(common_term), '\n')
  cat('Dim of O_matrix_nieobs:', dim(O_matrix_nieobs), '\n')
  cat('Dim of model$f_matrix_nieobs:', dim(model$f_matrix_nieobs), '\n')
  result_nieobs <- colSums(numerators / denominator)
  result_obs <- colSums(s_values_obs)
  test_weights_todoremove <- numerators / denominator
  cat('----')
  cat(theta)
  cat("\n")
  cat(result_obs)
  cat("\n")
  cat(result_nieobs)
  cat("\n")
  browser()
  return(result_nieobs + result_obs)
}



#
# step_func <- function(model, theta, O_matrix_nieobs) {
#   # Get dimensions
#   n_x1 <- nrow(model$x_1)
#   n_x0 <- nrow(model$x_0)
#   n_y1 <- length(model$y_1)
#   p <- length(theta)
#   # theta=c(0.66,-0.36)
#   # Define numerical stability constants
#   EPS <- .Machine$double.eps
#   MIN_DENOM <- 1e-10  # More conservative than 0.01
#
#   # Calculate common_term with better numerical handling
#   inv_C <- model$respondent_weights / as.vector(model$C_matrix_nieobs)
#   inv_C[!is.finite(inv_C) | inv_C < 0] <- 0  # Also handle negative weights
#
#   # Use more stable computation for common_term
#   common_term <- O_matrix_nieobs * model$f_matrix_nieobs *
#     rep(inv_C, each = nrow(O_matrix_nieobs)) *
#     model$respondent_weights
#
#   # Ensure common_term is non-negative and finite
#   common_term[common_term < 0 | !is.finite(common_term)] <- 0
#
#   # Calculate s_values_obs - vectorized approach
#   s_values_obs <- matrix(0, nrow = n_x1, ncol = p)
#   for (i in 1:n_x1) {
#     s_val_ij <- s_function(model, 1, model$x_1[, model$cols_delta], theta, i = i, j = i)
#     # Handle potential numerical issues in s_function output
#     s_val_ij[!is.finite(s_val_ij)] <- 0
#     s_values_obs[i, ] <- s_val_ij #* model$respondent_weights[i]
#   }
#
#   # Calculate s_values_unobs with improved stability
#   numerators <- matrix(0, nrow = n_x0, ncol = p)
#
#   for (i in 1:n_x0) {
#     row_total <- 0
#     for (j in 1:n_y1) {
#       if (common_term[i, j] > EPS) {  # Skip negligible terms
#         s_val_ij <- s_function(model, 0, model$x_0[, model$cols_delta], theta, i = i, j = j)
#         # s_val_ij[!is.finite(s_val_ij)] <- 0
#         # s_val_ij<-c(1,1)
#         numerators[i,] <- numerators[i, ] + common_term[i, j] * s_val_ij
#         row_total <- row_total + common_term[i, j]
#       }
#     }
#     # Store denominator for this row
#     if (i == 1) denominators <- numeric(n_x0)
#     denominators[i] <- row_total
#   }
#
#   # Safe denominator handling
#   denominators <- pmax(denominators, MIN_DENOM)
#
#   # Stable division
#   result_nieobs <- colSums(numerators / denominators)
#   result_obs <- colSums(s_values_obs)
#
#   # Final result with sanity checks
#   final_result <- result_nieobs + result_obs
#   # final_result<-abs(result_nieobs + result_obs)/min(abs(result_nieobs),abs(result_obs))
#   # final_result[!is.finite(final_result)] <- 0
#
#   # Optional: Remove debugging output for production
#
#     cat('Dim of numerators:', dim(numerators), '\n')
#     cat('Dim of denominators:', length(denominators), '\n')
#     cat('Dim of s_values_obs:', dim(s_values_obs), '\n')
#     cat('Dim of s_values_unobs:', dim(numerators), '\n')
#     cat('Dim of inv_C:', length(inv_C), '\n')
#     cat('Dim of common_term:', dim(common_term), '\n')
#     cat('Dim of O_matrix_nieobs:', dim(O_matrix_nieobs), '\n')
#     cat('Dim of model$f_matrix_nieobs:', dim(model$f_matrix_nieobs), '\n')
#
#
#
#     # final_result[1]<-final_result[1]+abs(theta[1])#abs?
#     cat('-----------')
#     cat("Theta:", theta, "\n")
#     cat("Result:", final_result, "\n")
#     print(result_obs)
#     print(result_nieobs)
#     browser()
#
#   return(final_result)
# }
