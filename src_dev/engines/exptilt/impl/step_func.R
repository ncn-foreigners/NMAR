#' @exportS3Method NULL
step_func <- function(model, theta, O_matrix_nieobs) {
  # Previous implementation computed
  #   s_values_unobs <- s_function(delta = 0, x = X_obs, theta)
  # and then multiplied by the fractional weight matrix O_{ij}f_{ij} / C_j).
  # That evaluates s(phi; 0, x_{1j}, y_j) at RESPONDENT covariates x_{1j},
  # not at the NONRESPONDENT covariates x_{1i}. In the ET E-step (discrete case),
  # the nonrespondent contribution must be
  #   E_0[s(phi; 0, x_{1i}, Y)] = sum_j w_{ij}s(phi; 0, x_{1i}, y_j),
  # i.e., hold x_{1i} fixed and average over the respondent grid y_j.
  #
  # Below we keep the vectorized structure but compute
  # s(phi; 0, x_{1i}, y_j) by replicating x_{1i} across the y_j grid and
  # calling s_function with an explicit y argument. This matches the papers
  # E-step, is algebraically correct, and preserves the performance of the
  # original vectorization.
  # We solve the mean score equation S_2(phi)=0 by combining respondent scores
  # with fractional-imputation terms for nonrespondents (eq. 13).
  # We explicitly build s(phi; delta, x_1, y) so that each nonrespondent i
  # uses their covariates x_{1i} while averaging over respondent outcomes y_j with
  # weights w_{ij} ~ O(x_{1i},y_j; phi)f_1(y_j | x_{1i}; gamma_har) / C(y_j; gamma_hat).

  x_obs_delta <- model$x_1[, model$cols_delta, drop = FALSE]
  y_obs <- model$y_1
  n_obs_rows <- if (is.null(dim(x_obs_delta))) length(x_obs_delta) else nrow(x_obs_delta)
  # Derive weights
  if (!is.null(model$design_weights) && !is.null(model$respondent_mask)) {
    respondent_weights <- model$design_weights[model$respondent_mask]
  } else {
    respondent_weights <- rep(1, n_obs_rows)
  }
  x_unobs_delta <- model$x_0[, model$cols_delta, drop = FALSE]
  n_non_rows <- if (is.null(dim(x_unobs_delta))) length(x_unobs_delta) else nrow(x_unobs_delta)
  if (!is.null(model$design_weights) && !is.null(model$respondent_mask)) {
    nonrespondent_weights <- model$design_weights[!model$respondent_mask]
  } else {
    nonrespondent_weights <- rep(1, n_non_rows)
  }

  # Respondent contribution: sum_j d_j s(\phi; 1, x_{1j}y_j)
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

  # Nonrespondent contribution (vectorised):
  # Build replicated design for all (i,j) pairs: i repeats, j varies fastest
  if (n_nonrespondents > 0) {
    x_rep <- model$x_0[, model$cols_delta, drop = FALSE]
    x_rep <- x_rep[rep(seq_len(n_nonrespondents), each = length(y_obs)), , drop = FALSE]
    y_rep <- rep(y_obs, times = n_nonrespondents)

    s_ij <- s_function(model, delta = 0, x = x_rep, theta = theta, y = y_rep)
    # Weight by fractional weights: flatten by row-major order (i blocks)
    w_flat <- as.vector(t(fractional_weights))
    s_ij_weighted <- s_ij * w_flat
    # Sum within each i block of length length(y_obs)
    grp <- rep(seq_len(n_nonrespondents), each = length(y_obs))
    S_non_by_i <- rowsum(s_ij_weighted, group = grp, reorder = FALSE)
    # Apply nonrespondent design weights and sum across i
    score_nonresp <- as.numeric(colSums(S_non_by_i * nonrespondent_weights))
  } else {
    score_nonresp <- rep(0, length(theta))
  }

  # Final score
  as.numeric(score_obs + score_nonresp)
}
