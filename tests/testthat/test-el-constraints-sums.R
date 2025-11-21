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
# Reconstruct components to compute raw constraint sums from stored diagnostics inputs
  spec <- el_prepare_inputs(
    formula = Y_miss ~ X,
    data = df,
    weights = NULL,
    n_total = NULL
  )
  dat2 <- spec$analysis_data
  obs_idx <- which(spec$respondent_mask)
  resp_df <- dat2[obs_idx, ]
  Z_un <- spec$missingness_design
  X_un <- spec$aux_design_full[spec$respondent_mask, , drop = FALSE]
  aux_means <- if (ncol(X_un) > 0) c(X = 0) else NULL
  aux_mat <- if (ncol(X_un) > 0) X_un else matrix(nrow = nrow(Z_un), ncol = 0)
  sc <- validate_and_apply_nmar_scaling(FALSE, ncol(aux_mat) > 0, Z_un, aux_mat, aux_means)
  Z <- sc$response_model_matrix_scaled
  Xc <- sc$auxiliary_matrix_scaled
  mu_x <- sc$mu_x_scaled
  n_resp_wt <- length(obs_idx)
  N_pop <- nrow(dat2)
  wts <- rep(1, length(obs_idx))

  fam <- logit_family()
  eq_fun <- el_build_equation_system(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
  beta_hat <- fit$model$coefficients
  z <- stats::qlogis(fit$model$nuisance$W_hat)
  lambda_hat <- if (!is.null(fit$model$nuisance$lambda_x)) as.numeric(fit$model$nuisance$lambda_x) else numeric(0)
  theta <- c(as.numeric(beta_hat), z, lambda_hat)
# Equations close to zero
  res <- as.numeric(eq_fun(theta))
  expect_lt(max(abs(res)), 1e-5)
})
