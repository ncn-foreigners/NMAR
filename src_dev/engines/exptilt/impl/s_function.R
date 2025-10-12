#' @exportS3Method NULL
s_function.nmar_exptilt <- function(model, delta, x, theta = model$theta) {

  theta_numeric <- as.numeric(theta)

  EPS <- .Machine$double.eps
  SAFE_THRESHOLD <- 1e-3

  if (delta == 1) {
    X_full_raw <- cbind(1, x)
  } else if (delta == 0) {
# browser()
    n_rows <- nrow(x)
    n_y1 <- length(model$y_1)

    X_expanded <- as.matrix(x)[rep(1:n_rows, each = n_y1), ]

    Y_expanded <- rep(as.numeric(model$y_1), n_rows)

    X_full_raw <- cbind(1, X_expanded, Y_expanded)
  } else {
    stop("delta must be 0 or 1 for vectorized calculation.")
  }

  X_full <- apply(X_full_raw, 2, as.numeric)

  eta_full <- as.vector(X_full %*% theta_numeric)

  pi_val_full <- model$family$linkinv(eta_full)
  pi_deriv_full <- model$family$mu.eta(eta_full)

  pi_val_safe_full <- pmin(pmax(pi_val_full, SAFE_THRESHOLD), 1 - SAFE_THRESHOLD)

  if (delta == 1) {
    Numerator <- pi_deriv_full * X_full
    Denominator <- pi_val_safe_full
  } else if (delta == 0) {
    Numerator <- -pi_deriv_full * X_full
    Denominator <- 1 - pi_val_safe_full
  }

  result_matrix <- Numerator / Denominator
# browser()
  return(result_matrix)
}
