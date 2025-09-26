#' @exportS3Method NULL
step_func <- function(model, theta, O_matrix_nieobs) {
  # Solve the mean score equation S2(phi) = 0 by combining respondent scores
  # with the fractional imputation terms for nonrespondents (see eq. 13 in the
  # paper). We explicitly build s(phi; delta, x1, y) so that each nonrespondent i
  # uses their covariates x1i while averaging over respondent outcomes yj with
  # weights w_ij proportional to O(x1i, yj; phi) * f1(yj | x1i; gamma_hat) / C(yj; gamma_hat)

  x_obs_delta <- model$x_1[, model$cols_delta, drop = FALSE]
  y_obs <- model$y_1
  n_obs_rows <- if (is.null(dim(x_obs_delta))) length(x_obs_delta) else nrow(x_obs_delta)
  respondent_weights <- model$respondent_weights
  if (is.null(respondent_weights) || !length(respondent_weights)) {
    respondent_weights <- rep(1, n_obs_rows)
  }
  x_unobs_delta <- model$x_0[, model$cols_delta, drop = FALSE]
  n_non_rows <- if (is.null(dim(x_unobs_delta))) length(x_unobs_delta) else nrow(x_unobs_delta)
  nonrespondent_weights <- model$nonrespondent_weights
  if (is.null(nonrespondent_weights) || !length(nonrespondent_weights)) {
    nonrespondent_weights <- rep(1, n_non_rows)
  }

  # Respondent contribution: sum_j d_j * s(phi; delta = 1, x1j, yj)
  s_obs <- s_function(model, delta = 1, x = x_obs_delta, theta = theta, y = y_obs)
  if (!is.null(s_obs) && length(respondent_weights)) {
    s_obs <- s_obs * respondent_weights
  }
  score_obs <- colSums(s_obs)

  n_nonrespondents <- n_non_rows
  if (n_nonrespondents == 0L) {
    return(as.numeric(score_obs))
  }

  # Fractional weights for each nonrespondent i over respondent outcomes j
  # common_term[i, j] = O_ij * f_ij / C_j
  inv_C <- as.vector(1 / model$C_matrix_nieobs)
  inv_C[!is.finite(inv_C)] <- 0
  common_term <- O_matrix_nieobs * model$f_matrix_nieobs
  common_term <- sweep(common_term, 2, inv_C, `*`)

  weight_denominator <- rowSums(common_term)
  weight_denominator[weight_denominator == 0] <- 1
  fractional_weights <- common_term / weight_denominator

  # Nonrespondent contribution: sum_i d_i * sum_j w_ij * s(phi; delta = 0, x1i, yj)
  score_nonresp <- numeric(length(theta))

  for (i in seq_len(n_nonrespondents)) {
    x_i_delta <- model$x_0[i, model$cols_delta, drop = FALSE]
    x_i_rep <- matrix(rep(x_i_delta, each = length(y_obs)), nrow = length(y_obs), byrow = TRUE)

    s_ij <- s_function(model, delta = 0, x = x_i_rep, theta = theta, y = y_obs)
    score_nonresp <- score_nonresp + nonrespondent_weights[i] * as.numeric(crossprod(fractional_weights[i, ], s_ij))
  }

  as.numeric(score_obs + score_nonresp)
}
