generate_Odds <- function(model) {
  # stopifnot(
  #   any(is.na(theta)) == FALSE,
  #   any(is.na(x)) == FALSE,
  #   any(is.na(y)) == FALSE,
  #   is.matrix(x) || is.data.frame(x),
  #   is.vector(y) || is.matrix(y),
  #   length(theta) == ncol(x) + 1 + 1 # +1 for intercept, +1 for y
  # )
  x_mat <- as.matrix(model$x_0[,model$cols_delta,drop=FALSE])
  y_vec <- as.vector(model$y_1)


  i_indices <- rep(1:nrow(x_mat), each = length(y_vec))
  j_indices <- rep(1:length(y_vec), times = nrow(x_mat))


  x_expanded <- x_mat[i_indices, , drop = FALSE]
  y_expanded <- y_vec[j_indices]


  x_aug <- cbind(1, x_expanded, y_expanded)  #+Intercept, y
  eta <- x_aug %*% model$theta


  if(model$prob_model_type == "logit") {
    odds <- exp(-eta)
  } else if(model$prob_model_type == "probit") {
    odds <- pnorm(-eta)
  }


  matrix(odds, nrow = nrow(x_mat), ncol = length(y_vec), byrow = FALSE)
}
