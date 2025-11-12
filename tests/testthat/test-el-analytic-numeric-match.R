test_that("analytic vs numeric Jacobian agree across links and scaling", {
  skip_if_not_installed("numDeriv")
  set.seed(2025)
# Generate data with auxiliaries and extra response predictor
  n <- 250
  X1 <- rnorm(n); X2 <- rnorm(n); Z <- rnorm(n)
  Y <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
  p <- plogis(-0.5 + 0.6 * scale(Y)[, 1] + 0.4 * Z)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X1 = X1, X2 = X2, Z = Z)
  df$Y_miss[!R] <- NA_real_

# Helper to check one configuration
  check_one <- function(fam, standardize) {
    parsed <- NMAR:::el_prepare_inputs(Y_miss ~ X1 + X2 | Z, df)
    dat2 <- parsed$data
    resp_var <- parsed$delta_name
    obs_idx <- which(dat2[[resp_var]] == 1)
    resp_df <- dat2[obs_idx, , drop = FALSE]
    Z_un <- parsed$response_matrix
    X_un <- parsed$auxiliary_matrix
    aux_means <- setNames(rep(0, ncol(X_un)), colnames(X_un))
    sc <- NMAR:::validate_and_apply_nmar_scaling(standardize, TRUE, Z_un, X_un, aux_means, weights = rep(1, nrow(resp_df)))
    Z <- sc$response_model_matrix_scaled
    Xc <- sc$auxiliary_matrix_scaled
    mu_x <- sc$mu_x_scaled
    n_resp_wt <- nrow(resp_df)
    N_pop <- nrow(dat2)
    wts <- rep(1, n_resp_wt)
    eq_fun <- NMAR:::el_build_equation_system(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
    jac_fun <- NMAR:::el_build_jacobian(fam, Z, Xc, wts, N_pop, n_resp_wt, mu_x)
# Solve to the root using a better init and controls; allow Broyden fallback
    W0 <- sum(wts) / N_pop
    W0 <- min(max(W0, 1e-12), 1 - 1e-12)
    z0 <- stats::qlogis(W0)
    init <- c(rep(0, ncol(Z)), z0, rep(0, ncol(Xc)))
    sol <- nleqslv::nleqslv(x = init, fn = eq_fun, jac = jac_fun, method = "Newton", control = list(ftol = 1e-10, xtol = 1e-10, maxit = 200))
    if (sol$termcd > 2) {
      sol <- nleqslv::nleqslv(x = init, fn = eq_fun, method = "Broyden", control = list(ftol = 1e-10, xtol = 1e-10, maxit = 200))
    }
    expect_lte(sol$termcd, 3)
    theta <- sol$x
# Compare analytic vs numeric J at the solution
    J_num <- numDeriv::jacobian(eq_fun, theta)
    J_ana <- jac_fun(theta)
    denom <- max(1e-8, norm(J_num, type = "F"))
    rel_diff <- norm(J_ana - J_num, type = "F") / denom
    expect_lt(rel_diff, 1e-4)
  }

  check_one(NMAR:::logit_family(), TRUE)
  check_one(NMAR:::logit_family(), FALSE)
  check_one(NMAR:::probit_family(), TRUE)
  check_one(NMAR:::probit_family(), FALSE)
})
