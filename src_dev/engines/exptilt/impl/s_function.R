

#' @exportS3Method NULL
s_function.nmar_exptilt <- function(model, delta, x, theta = model$theta, i = NULL, j = NULL) {
  # Define machine precision constant
  EPS <- .Machine$double.eps
  # Safe threshold to avoid extreme values
  SAFE_THRESHOLD <- 1e-3

  # If i and j are provided, calculate for specific indices
  if (!is.null(i) && !is.null(j)) {
    # Use specific row i from x and element j from y_1
    x_vec <- c(1, as.numeric(x[i, ]), model$y_1[j])

    eta <- as.vector(sum(x_vec * theta))
    pi_val <- model$family$linkinv(eta)
    pi_deriv <- model$family$mu.eta(eta)

    # Safe clamping of pi_val to prevent division by zero or near-zero
    pi_val_safe <- pmin(pmax(pi_val, SAFE_THRESHOLD), 1 - SAFE_THRESHOLD)

    if (delta == 1) {
      numerator <- pi_deriv * x_vec
      denominator <- pi_val_safe
    } else if (delta == 0) {
      numerator <- -pi_deriv * x_vec
      denominator <- 1 - pi_val_safe
    } else {
      stop("delta must be either 0 or 1")
    }


    # result <- (delta - pi_val_safe)*x_vec
    result <- numerator / denominator
    # result <- (delta-pi_val_safe)*pi_deriv
    # Use more precise zero detection
    # result[abs(result) < EPS | is.nan(result)] <- 0
    # cat('s function outtput')
    # print(as.numeric(result))
    # browser()
    return(as.numeric(result))
  }
  # Original behavior when no i,j provided
  else {
    x_mat <- as.matrix(x)
    x_aug <- cbind(1, x_mat, model$y_1[1:nrow(x_mat)])

    eta <- as.vector(x_aug %*% theta)
    pi_vals <- model$family$linkinv(eta)
    pi_deriv <- model$family$mu.eta(eta) * x_aug

    # Improved clamping with smooth boundaries
    pi_vals_safe <- pmin(pmax(pi_vals, SAFE_THRESHOLD), 1 - SAFE_THRESHOLD)

    if (delta == 1) {
      numerator <- pi_deriv
      denominator <- pi_vals_safe
    } else if (delta == 0) {
      numerator <- -pi_deriv
      denominator <- 1 - pi_vals_safe
    } else {
      stop("delta must be either 0 or 1")
    }

    result <- numerator / denominator
    # Better handling of numerical issues
    result[abs(result) < EPS | is.nan(result) | is.infinite(result)] <- 0

    return(result)
  }
}
