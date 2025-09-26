#' @exportS3Method NULL
s_function.nmar_exptilt <- function(model, delta, x, theta = model$theta, y = NULL) {
  # We compute the score contributions d(log-likelihood)/d phi for each row by
  # evaluating the Bernoulli score with respect to eta (using the shared family
  # helpers) and multiplying by the design vector [1, x1, y]. This mirrors
  # s(phi; delta, x1, y) in Riddles et al.

  x_mat <- as.matrix(x)
  n_rows <- nrow(x_mat)
  if (!n_rows) {
    return(matrix(0, nrow = 0, ncol = length(theta)))
  }

  if (is.null(y)) {
    # Previous implementation silently pulled the first n_rows respondent y's.
    # Retain that fallback for existing callers while allowing explicit y input
    # for the EM E-step (nonrespondent / respondent pairing)
    y_vec <- model$y_1[seq_len(n_rows)]
  } else {
    y_vec <- as.numeric(y)
    if (length(y_vec) == 1L && n_rows > 1L) {
      y_vec <- rep(y_vec, n_rows)
    }
    if (length(y_vec) != n_rows) {
      stop("Length of `y` must match the number of rows in `x`.", call. = FALSE)
    }
  }

  delta_vec <- delta
  if (length(delta_vec) == 1L && n_rows > 1L) {
    delta_vec <- rep(delta_vec, n_rows)
  }
  if (length(delta_vec) != n_rows) {
    stop("Length of `delta` must match the number of rows in `x`.", call. = FALSE)
  }

  design_mat <- cbind(1, x_mat, y_vec)
  storage.mode(design_mat) <- "numeric"
  theta_num <- as.numeric(theta)
  if (length(theta_num) != ncol(design_mat)) {
    stop("Length of `theta` does not match the design vector in s_function.", call. = FALSE)
  }
  if (!is.null(names(theta))) {
    colnames(design_mat) <- names(theta)
  }

  eta <- as.vector(design_mat %*% theta_num)
  score_eta <- model$family$score_eta(eta, delta_vec)

  # Gradient with respect to phi is score_eta * d eta / d phi, and d eta / d phi equals design_mat
  result <- score_eta * design_mat
  result[!is.finite(result)] <- 0
  result
}
