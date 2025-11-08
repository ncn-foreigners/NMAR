test_that("constraint sums are near zero at solution (no trimming)", {
  df <- make_iid_nmar(n = 300, alpha = 0.6, seed = 1234)

  fit <- nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = make_engine(auxiliary_means = c(X = 0), trim_cap = Inf, variance_method = "none", standardize = FALSE)
  )
  expect_true(fit$converged)
# Jacobian quality reported
  diag <- fit$diagnostics
  expect_true(is.finite(diag$jacobian_condition_number) || is.na(diag$jacobian_condition_number))
  engine <- make_engine(
    auxiliary_means = c(X = 0),
    trim_cap = Inf,
    variance_method = "none",
    standardize = FALSE
  )
  runtime <- build_el_runtime(Y_miss ~ X, df, engine)
  dat2 <- runtime$data
  fmls <- runtime$internal_formula
  resp_var <- runtime$response_var
  obs_idx <- runtime$observed_indices
  resp_df <- dat2[obs_idx, ]
  Z_un <- model.matrix(update(fmls$response, NULL ~ .), data = resp_df)
  X_un <- model.matrix(fmls$auxiliary, data = resp_df)
  sc <- NMAR:::validate_and_apply_nmar_scaling(FALSE, !is.null(fmls$auxiliary), Z_un, if (is.null(fmls$auxiliary)) matrix(nrow = nrow(Z_un), ncol = 0) else X_un, if (is.null(fmls$auxiliary)) NULL else c(X = 0))
  Z <- sc$response_model_matrix_scaled
  Xc <- sc$auxiliary_matrix_scaled
  mu_x <- sc$mu_x_scaled
  n_resp_wt <- length(obs_idx)
  N_pop <- nrow(dat2)
  wts <- rep(1, length(obs_idx))

  fam <- NMAR:::logit_family()
  eq_fun <- NMAR:::el_build_equation_system(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
  beta_hat <- fit$model$coefficients
  z <- stats::qlogis(fit$model$nuisance$W_hat)
  lambda_hat <- if (!is.null(fit$model$nuisance$lambda_x)) as.numeric(fit$model$nuisance$lambda_x) else numeric(0)
  theta <- c(as.numeric(beta_hat), z, lambda_hat)
# Equations close to zero
  res <- as.numeric(eq_fun(theta))
  expect_lt(max(abs(res)), 1e-5)
})
