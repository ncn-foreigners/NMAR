#' @exportS3Method NULL
s_function.nmar_exptilt <- function(model, delta, x, theta = model$theta, y = NULL, ...) {
  # Score contributions are computed row-wise as
  #   d(log L) / d phi using a vectorized numerator/denominator form, branching on delta:
  #   - delta = 1:  result = mu_eta(eta)[1,x,y] / pi(eta)
  #   - delta = 0:  result = -mu_eta(eta)[1,x,y] / (1-pi(eta))
  #   where eta = [1,x,y]^T theta, pi(eta)=linkinv(eta), and mu_eta(eta) = d pi/d eta.
  # We accept an explicit y so the caller (step_func) can evaluate
  #   s(phi; delta, x, y) on the respondent y-grid while holding a
  #   nonrespondent x fixed (E-step with delta=0).
  #
  # changes:
  # - Previously, this helper implicitly paired y with the rows of x coming
  #   from respondents, i.e., used y_j alongside x_{1j}. When the E-step
  #   needs s(phi; 0, x_{1i}, y_j) (hold nonrespondent x_{1i} fixed and
  #   average over the respondent grid y_j), callers could not supply y_j
  #   independent of x. That led downstream code to evaluate the wrong score
  #   argument for delta=0 (using x_{1j} instead of x_{1i}).
  # - We now expose an explicit y argument so callers can form cartesian
  #   grid (x_{1i}, y_j) correctly for the discrete E-step, while keeping all
  #   computations vectorized

  x_mat <- as.matrix(x)
  n_rows <- nrow(x_mat)
  if (!n_rows) return(matrix(0, nrow = 0, ncol = length(theta)))

  # Pair y with rows of x
  if (is.null(y)) {
    y_vec <- model$y_1[seq_len(n_rows)]
  } else {
    y_vec <- as.numeric(y)
    if (length(y_vec) == 1L && n_rows > 1L) y_vec <- rep(y_vec, n_rows)
    if (length(y_vec) != n_rows) stop("Length of `y` must match `nrow(x)`.", call. = FALSE)
  }

  # delta as a row-wise vector (allow scalar)
  delta_vec <- delta
  if (length(delta_vec) == 1L && n_rows > 1L) delta_vec <- rep(delta_vec, n_rows)
  if (length(delta_vec) != n_rows) stop("Length of `delta` must match `nrow(x)`.", call. = FALSE)

  # Design and linear predictor
  design_mat <- cbind(1, x_mat, y_vec)
  storage.mode(design_mat) <- "numeric"
  if (!is.null(names(theta))) colnames(design_mat) <- names(theta)
  eta <- as.vector(design_mat %*% as.numeric(theta))

  # Mean and derivative wrt eta
  p <- model$family$linkinv(eta)
  mu_eta <- model$family$mu.eta(eta)

  # Guard against division by 0
  eps <- .Machine$double.eps
  p <- pmin(pmax(p, eps), 1 - eps)

  # pi'(eta),[1,x,y]
  pi_deriv_times_design <- design_mat * as.numeric(mu_eta)

  # Branch on delta using the optimal numerator/denominator form
  is_one <- (delta_vec == 1)
  is_zero <- !is_one
  # Initialise result matrix
  result <- matrix(0, nrow = n_rows, ncol = ncol(design_mat))
  if (any(is_one)) {
    result[is_one, ] <- pi_deriv_times_design[is_one, , drop = FALSE] / p[is_one]
  }
  if (any(is_zero)) {
    result[is_zero, ] <- -pi_deriv_times_design[is_zero, , drop = FALSE] / (1 - p[is_zero])
  }

  result[!is.finite(result)] <- 0
  # Preserve column names for downstream comparisons/tests
  if (!is.null(colnames(design_mat))) colnames(result) <- colnames(design_mat)
  result
}
