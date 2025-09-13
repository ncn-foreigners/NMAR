test_that("constraint sums are near zero at solution (no trimming)", {
  set.seed(1234)
  N <- 300
  X <- rnorm(N)
  Y <- 1 + 0.5 * X + rnorm(N)
  p <- plogis(-0.5 + 0.6 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[!R] <- NA_real_

  fit <- nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = el_engine(auxiliary_means = c(X = 0), trim_cap = Inf, variance_method = "delta", standardize = FALSE)

  )
  expect_true(fit$converged)
  # Reconstruct components to compute raw constraint sums from stored diagnostics inputs
  parsed <- nmar:::prepare_el_inputs(Y_miss ~ X, df, NULL)
  dat2 <- parsed$data
  fmls <- parsed$formula_list
  resp_var <- all.vars(fmls$response)[1]
  obs_idx <- which(dat2[[resp_var]] == 1)
  resp_df <- dat2[obs_idx, ]
  Z_un <- model.matrix(update(fmls$response, NULL ~ .), data = resp_df)
  X_un <- model.matrix(fmls$auxiliary, data = resp_df)
  sc <- nmar:::validate_and_apply_nmar_scaling(FALSE, !is.null(fmls$auxiliary), Z_un, if (is.null(fmls$auxiliary)) matrix(nrow = nrow(Z_un), ncol = 0) else X_un, if (is.null(fmls$auxiliary)) NULL else c(X = 0))
  Z <- sc$response_model_matrix_scaled
  Xc <- sc$auxiliary_matrix_scaled
  mu_x <- sc$mu_x_scaled
  n_resp_wt <- length(obs_idx)
  N_pop <- nrow(dat2)
  wts <- rep(1, length(obs_idx))

  fam <- nmar:::logit_family()
  eq_fun <- nmar:::build_equation_system(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
  beta_hat <- fit$coefficients$response_model
  z <- stats::qlogis(fit$coefficients$nuisance$W_hat)
  lambda_hat <- if (!is.null(fit$coefficients$nuisance$lambda_x)) as.numeric(fit$coefficients$nuisance$lambda_x) else numeric(0)
  theta <- c(as.numeric(beta_hat), z, lambda_hat)
  # Equations close to zero
  res <- as.numeric(eq_fun(theta))
  expect_lt(max(abs(res)), 1e-5)
})
