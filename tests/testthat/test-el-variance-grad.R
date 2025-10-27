test_that("analytic vs numeric gradient for g are close when untrimmed", {
  set.seed(1234)
  n <- 80
  Y_full <- rbinom(n, 1, 0.5)
  R <- rbinom(n, 1, 0.8)
  Y_miss <- Y_full
  Y_miss[R == 0] <- NA
  df <- data.frame(Y = Y_miss, X = rnorm(n))
# Build inputs
  parsed <- NMAR:::prepare_el_inputs(Y ~ X, df)
  fmls <- parsed$formula_list
  delta_var <- all.vars(fmls$response)[1]
  resp_idx <- which(parsed$data[[delta_var]] == 1)
  resp_df <- parsed$data[resp_idx, ]
  Z_un <- model.matrix(update(fmls$response, NULL ~ .), data = resp_df)
  X_un <- matrix(nrow = nrow(resp_df), ncol = 0)
  sc <- NMAR:::validate_and_apply_nmar_scaling(TRUE, FALSE, Z_un, X_un, NULL, weights = rep(1, nrow(resp_df)))
  Z <- sc$response_model_matrix_scaled
  Xc <- sc$auxiliary_matrix_scaled
  mu_x <- sc$mu_x_scaled
  wts <- rep(1, nrow(resp_df))
  N_pop <- nrow(parsed$data)
  fam <- NMAR:::logit_family()
  eq_fun <- NMAR:::el_build_equation_system(fam, Z, Xc, wts, N_pop, sum(wts), mu_x)
  jac_fun <- NMAR:::build_el_jacobian(fam, Z, Xc, wts, N_pop, sum(wts), mu_x)
# Solve to get estimates
  init <- c(rep(0, ncol(Z)), 0) # no auxiliaries
  sol <- nleqslv::nleqslv(x = init, fn = eq_fun, jac = jac_fun, method = "Newton")
  expect_lte(sol$termcd, 2)
  est <- sol$x
# Post-solution useful quantities
  eta <- as.vector(Z %*% est[1:ncol(Z)])
  w_i <- fam$linkinv(eta)
  W <- plogis(est[ncol(Z) + 1])
  lambda_W <- ((N_pop / sum(wts)) - 1) / (1 - W)
  denom <- 1 + lambda_W * (w_i - W)
  denom <- pmax(as.numeric(denom), 1e-8)
# Analytic gradient
  g_anal <- NMAR:::el_grad_g_analytic(fam, Z, Xc, mu_x, wts, eta, w_i, W, denom, lambda_W, resp_df$Y, sum(wts), N_pop)
# Numeric gradient via mean function builder
  g_fn <- NMAR:::el_build_mean_fn(fam, Z, Xc, mu_x, wts, N_pop, sum(wts), Inf, resp_df$Y)
  g_num <- NMAR:::grad_numeric(est, g_fn)
# Attach names to numeric gradient to match analytic ordering for comparison
  names(g_num) <- names(g_anal)
  expect_equal(as.numeric(g_anal), as.numeric(g_num), tolerance = 1e-4)
})
