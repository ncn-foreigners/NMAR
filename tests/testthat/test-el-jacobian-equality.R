test_that("analytic Jacobian matches numeric Jacobian at solution (logit and probit)", {
  skip_if_not_installed("numDeriv")
  df <- make_iid_nmar(n = 200, alpha = 0.5, seed = 3401)
# Fit once (logit default) to get a stable solution near the root
  fit <- el.data.frame(df, Y_miss ~ X,
    auxiliary_means = c(X = 0), standardize = TRUE,
    variance_method = "none"
  )
  expect_type(fit$converged, "logical")

  spec <- el_build_input_spec(
    formula = Y_miss ~ X,
    data = df,
    weights_full = NULL,
    population_total = NULL,
    is_survey = FALSE,
    design_object = NULL,
    auxiliary_means = NULL
  )
  dat2 <- if (inherits(spec$analysis_object, "survey.design")) spec$analysis_object$variables else spec$analysis_object
  obs_idx <- spec$respondent_indices
  resp_df <- dat2[obs_idx, ]
  Z_un <- spec$missingness_design
  X_un <- spec$auxiliary_design_full[spec$respondent_mask, , drop = FALSE]
  aux_means <- if (ncol(X_un) > 0) setNames(rep(0, ncol(X_un)), colnames(X_un)) else NULL
  aux_mat <- if (ncol(X_un) > 0) X_un else matrix(nrow = nrow(Z_un), ncol = 0)
  sc <- validate_and_apply_nmar_scaling(FALSE, ncol(aux_mat) > 0, Z_un, aux_mat, aux_means)
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

  for (fam in list(logit_family(), probit_family())) {
    eq_fun <- el_build_equation_system(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
    jac_fun <- el_build_jacobian(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
    J_num <- numDeriv::jacobian(eq_fun, theta)
    J_ana <- jac_fun(theta)
    denom <- max(1e-8, norm(J_num, type = "F"))
    rel_diff <- norm(J_ana - J_num, type = "F") / denom
    expect_lt(rel_diff, 1e-4)
  }
})
