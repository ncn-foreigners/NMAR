test_that("s_function matches analytic score for logit", {
  model <- list(
    family = nmar:::logit_family(),
    y_1 = c(0.2, -0.1, 0.5)
  )
  class(model) <- "nmar_exptilt"

  x <- matrix(c(0.3, -0.4, 0.7, 1.2, -0.9, 0.1), ncol = 2)
  colnames(x) <- c("x1", "x2")
  theta <- c(`(Intercept)` = 0.1, x1 = -0.2, x2 = 0.3, y = 0.4)
  delta <- c(1, 0, 1)

  score_pkg <- nmar:::s_function.nmar_exptilt(model, delta = delta, x = x, theta = theta, y = model$y_1)

  design_mat <- cbind(1, x, model$y_1)
  eta <- as.vector(design_mat %*% theta)
  p <- stats::plogis(eta)
  score_manual <- (delta - p) * design_mat

  colnames(score_manual) <- colnames(score_pkg)
  expect_equal(score_pkg, score_manual, tolerance = 1e-10)
})

test_that("s_function matches analytic score for probit", {
  model <- list(
    family = nmar:::probit_family(),
    y_1 = c(-0.3, 1.1)
  )
  class(model) <- "nmar_exptilt"

  x <- matrix(c(0.5, -0.8, 0.2, 0.6), ncol = 2)
  colnames(x) <- c("x1", "x2")
  theta <- c(`(Intercept)` = -0.1, x1 = 0.4, x2 = -0.2, y = 0.3)
  delta <- c(1, 0)

  score_pkg <- nmar:::s_function.nmar_exptilt(model, delta = delta, x = x, theta = theta, y = model$y_1)

  design_mat <- cbind(1, x, model$y_1)
  eta <- as.vector(design_mat %*% theta)
  phi <- stats::dnorm(eta)
  Phi <- stats::pnorm(eta)
  Phi <- pmin(pmax(Phi, .Machine$double.eps), 1 - .Machine$double.eps)
  odds_ratio <- delta
  odds_ratio[delta == 1] <- phi[delta == 1] / Phi[delta == 1]
  odds_ratio[delta == 0] <- -phi[delta == 0] / (1 - Phi[delta == 0])
  score_manual <- odds_ratio * design_mat

  colnames(score_manual) <- colnames(score_pkg)
  expect_equal(score_pkg, score_manual, tolerance = 1e-10)
})

test_that("fractional weights sum to one", {
  model <- list(
    cols_delta = c("x1", "x2"),
    col_y = "y",
    y_1 = c(0.1, -0.2, 0.5),
    x_1 = as.matrix(data.frame(x1 = c(0.2, -0.3, 0.4), x2 = c(0.5, -0.6, 0.7))),
    x_0 = as.matrix(data.frame(x1 = c(-0.4, 0.8), x2 = c(0.1, -0.2))),
    x_for_y_obs = as.matrix(data.frame(x1 = c(0.2, -0.3, 0.4), x2 = c(0.5, -0.6, 0.7))),
    x_for_y_unobs = as.matrix(data.frame(x1 = c(-0.4, 0.8), x2 = c(0.1, -0.2))),
    theta = c(`(Intercept)` = 0.3, x1 = -0.2, x2 = 0.1, y = 0.4),
    family = nmar:::logit_family(),
    respondent_weights = rep(1, 3),
    nonrespondent_weights = rep(1, 2),
    y_dens = "normal"
  )
  class(model) <- "nmar_exptilt"

  model$prob_model_type <- "logit"
  model$density_fun <- function(y, x) stats::dnorm(y, mean = 0, sd = 1)
  model$O_matrix_nieobs <- nmar:::generate_Odds.nmar_exptilt(model)
  n_non <- nrow(model$x_for_y_unobs)
  n_resp <- length(model$y_1)
  model$f_matrix_nieobs <- matrix(0, nrow = n_non, ncol = n_resp)
  for (i_idx in seq_len(n_non)) {
    for (j_idx in seq_len(n_resp)) {
      model$f_matrix_nieobs[i_idx, j_idx] <- model$density_fun(model$y_1[j_idx], model$x_for_y_unobs[i_idx, , drop = FALSE])
    }
  }
  C_vec <- numeric(n_resp)
  for (j_idx in seq_len(n_resp)) {
    dens_vals <- numeric(nrow(model$x_for_y_obs))
    for (l in seq_len(nrow(model$x_for_y_obs))) {
      dens_vals[l] <- model$density_fun(model$y_1[j_idx], model$x_for_y_obs[l, , drop = FALSE])
    }
    C_vec[j_idx] <- sum(model$respondent_weights * dens_vals)
  }
  model$C_matrix_nieobs <- matrix(C_vec, ncol = 1)

  inv_C <- as.vector(1 / model$C_matrix_nieobs)
  weights_raw <- model$O_matrix_nieobs * model$f_matrix_nieobs
  weights_raw <- sweep(weights_raw, 2, inv_C, `*`)
  normalized <- sweep(weights_raw, 1, rowSums(weights_raw), `/`)

  expect_equal(rowSums(normalized), rep(1, nrow(model$x_0)), tolerance = 1e-8)
})

test_that("generate_Odds matches closed form", {
  model <- list(
    x_0 = as.matrix(data.frame(x1 = c(-0.4, 0.6), x2 = c(0.3, -0.5))),
    y_1 = c(0.2, -1.1, 0.7),
    cols_delta = c("x1", "x2"),
    theta = c(`(Intercept)` = 0.2, x1 = -0.3, x2 = 0.4, y = -0.5)
  )
  class(model) <- "nmar_exptilt"

  model$prob_model_type <- "logit"
  odds_logit <- nmar:::generate_Odds.nmar_exptilt(model)
  n_non <- nrow(model$x_0)
  n_resp <- length(model$y_1)
  i_indices <- rep(seq_len(n_non), each = n_resp)
  j_indices <- rep(seq_len(n_resp), times = n_non)
  odds_vec <- numeric(length(i_indices))
  for (idx in seq_along(i_indices)) {
    i <- i_indices[idx]
    j <- j_indices[idx]
    x_aug <- c(1, model$x_0[i, ], model$y_1[j])
    eta <- sum(x_aug * model$theta)
    odds_vec[idx] <- exp(-eta)
  }
  odds_manual <- matrix(odds_vec, nrow = n_non, ncol = n_resp, byrow = FALSE)
  expect_equal(odds_logit, odds_manual, tolerance = 1e-10)

  model$prob_model_type <- "probit"
  odds_probit <- nmar:::generate_Odds.nmar_exptilt(model)
  odds_vec_probit <- numeric(length(i_indices))
  for (idx in seq_along(i_indices)) {
    i <- i_indices[idx]
    j <- j_indices[idx]
    x_aug <- c(1, model$x_0[i, ], model$y_1[j])
    eta <- sum(x_aug * model$theta)
    log_p <- stats::pnorm(eta, log.p = TRUE)
    log_tail <- stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE)
    odds_vec_probit[idx] <- exp(log_tail - log_p)
  }
  odds_manual_probit <- matrix(odds_vec_probit, nrow = n_non, ncol = n_resp, byrow = FALSE)
  expect_equal(odds_probit, odds_manual_probit, tolerance = 1e-10)
})
