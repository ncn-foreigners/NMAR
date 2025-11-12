test_that("analytic Jacobian matches numeric Jacobian at solution (logit and probit)", {
  skip_if_not_installed("numDeriv")
  df <- make_iid_nmar(n = 200, alpha = 0.5, seed = 3401)
# Fit once (logit default) to get a stable solution near the root
  fit <- NMAR:::el.data.frame(df, Y_miss ~ X,
    auxiliary_means = c(X = 0), standardize = TRUE,
    variance_method = "none"
  )
  expect_type(fit$converged, "logical")

  parsed <- NMAR:::el_prepare_inputs(Y_miss ~ X, df)
  dat2 <- parsed$data
  resp_var <- parsed$delta_name
  obs_idx <- which(dat2[[resp_var]] == 1)
  resp_df <- dat2[obs_idx, ]
  Z_un <- parsed$response_matrix
  X_un <- parsed$auxiliary_matrix
  aux_means <- c(X = 0)
  sc <- NMAR:::validate_and_apply_nmar_scaling(FALSE, parsed$has_aux, Z_un, if (parsed$has_aux) X_un else matrix(nrow = nrow(Z_un), ncol = 0), if (parsed$has_aux) aux_means else NULL)
  Z <- sc$response_model_matrix_scaled
  Xc <- sc$auxiliary_matrix_scaled
  mu_x <- sc$mu_x_scaled
  n_resp_wt <- length(obs_idx)
  N_pop <- nrow(dat2)
  wts <- rep(1, length(obs_idx))

  beta_hat <- fit$model$coefficients # unscaled since standardize=FALSE
  z <- stats::qlogis(fit$model$nuisance$W_hat)
  lambda_hat <- if (!is.null(fit$model$nuisance$lambda_x)) as.numeric(fit$model$nuisance$lambda_x) else numeric(0)
  theta <- c(as.numeric(beta_hat), z, lambda_hat)

  for (fam in list(NMAR:::logit_family(), NMAR:::probit_family())) {
    eq_fun <- NMAR:::el_build_equation_system(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
    jac_fun <- NMAR:::el_build_jacobian(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
    J_num <- numDeriv::jacobian(eq_fun, theta)
    J_ana <- jac_fun(theta)
    denom <- max(1e-8, norm(J_num, type = "F"))
    rel_diff <- norm(J_ana - J_num, type = "F") / denom
    expect_lt(rel_diff, 1e-4)
  }
})
