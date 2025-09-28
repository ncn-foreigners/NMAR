#' @exportS3Method NULL
generate_Odds.nmar_exptilt <- function(model, ...) {
  x_mat <- as.matrix(model$x_0[, model$cols_delta, drop = FALSE])
  y_vec <- as.vector(model$y_1)


  i_indices <- rep(1:nrow(x_mat), each = length(y_vec))
  j_indices <- rep(1:length(y_vec), times = nrow(x_mat))


  x_expanded <- x_mat[i_indices, , drop = FALSE]
  y_expanded <- y_vec[j_indices]


  x_aug <- cbind(1, x_expanded, y_expanded) #+Intercept, y
  eta <- x_aug %*% model$theta

  # Conditional odds O(x, y) = Pr(delta = 0 | x, y) / Pr(delta = 1 | x, y)
  if (model$prob_model_type == "logit") {
    # For the logit link, O = exp(-eta) exactly, but we clamp p for stability so
    # the implementation is consistent across links
    p <- stats::plogis(eta)
    p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
    odds <- (1 - p) / p
  } else if (model$prob_model_type == "probit") {
    # Use log-space to avoid catastrophic cancellation in tails:
    #   log O = log(1 - Phi(eta)) - log Phi(eta)
    log_p <- stats::pnorm(eta, log.p = TRUE)
    log_tail <- stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE)
    log_odds <- log_tail - log_p
    odds <- exp(log_odds)
  } else {
    stop("Unsupported prob_model_type for odds generation.", call. = FALSE)
  }

  matrix(odds, nrow = nrow(x_mat), ncol = length(y_vec), byrow = FALSE)
}
