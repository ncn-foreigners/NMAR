test_that("analytic Jacobian matches numeric Jacobian at solution (standardize=FALSE)", {
  skip_if_not_installed("numDeriv")
  set.seed(3401)
  N <- 200
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-1 + 0.5 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_

  # Fit once to get parameter estimates
  fit <- nmar:::el.data.frame(df, Y_miss ~ X,
    response_predictors = NULL,
    auxiliary_means = c(X = 0), standardize = FALSE,
    variance_method = "delta"
  )
  expect_true(fit$converged)

  # Reconstruct inputs used by builders
  parsed <- nmar:::prepare_el_inputs(Y_miss ~ X, df, NULL)
  dat2 <- parsed$data
  fmls <- parsed$formula_list
  resp_var <- all.vars(fmls$response)[1]
  obs_idx <- which(dat2[[resp_var]] == 1)
  resp_df <- dat2[obs_idx, ]
  Z_un <- model.matrix(update(fmls$response, NULL ~ .), data = resp_df)
  X_un <- model.matrix(fmls$auxiliary, data = resp_df)
  aux_means <- c(X = 0)
  sc <- nmar:::validate_and_apply_nmar_scaling(FALSE, !is.null(fmls$auxiliary), Z_un, if (is.null(fmls$auxiliary)) matrix(nrow = nrow(Z_un), ncol = 0) else X_un, if (is.null(fmls$auxiliary)) NULL else aux_means)

  Z <- sc$response_model_matrix_scaled
  Xc <- sc$auxiliary_matrix_scaled
  mu_x <- sc$mu_x_scaled
  n_resp_wt <- length(obs_idx)
  N_pop <- nrow(dat2)
  wts <- rep(1, length(obs_idx))

  eq_fun <- nmar:::build_equation_system(nmar:::logit_family(), Z, Xc, wts, N_pop, n_resp_wt, mu_x)
  jac_fun <- nmar:::build_el_jacobian(nmar:::logit_family(), Z, Xc, wts, N_pop, n_resp_wt, mu_x)

  beta_hat <- fit$model$coefficients # unscaled since standardize=FALSE
  z <- stats::qlogis(fit$model$nuisance$W_hat)
  lambda_hat <- if (!is.null(fit$model$nuisance$lambda_x)) as.numeric(fit$model$nuisance$lambda_x) else numeric(0)
  theta <- c(as.numeric(beta_hat), z, lambda_hat)

  J_num <- numDeriv::jacobian(eq_fun, theta)
  J_ana <- jac_fun(theta)

  denom <- max(1e-8, norm(J_num, type = "F"))
  rel_diff <- norm(J_ana - J_num, type = "F") / denom
  expect_lt(rel_diff, 1e-4)
})
