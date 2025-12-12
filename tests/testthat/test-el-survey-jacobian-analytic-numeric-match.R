test_that("survey analytic Jacobian matches numeric Jacobian (with auxiliaries)", {
  skip_if_not_installed("survey")
  skip_if_not_installed("numDeriv")

  set.seed(4242)
  n <- 180
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 0.2 + 0.6 * X + rnorm(n)

  p <- plogis(-0.4 + 0.7 * scale(Y)[, 1] + 0.2 * Z)
  R <- runif(n) < p
  if (all(R)) R[sample.int(n, 1)] <- FALSE
  if (!any(R)) R[sample.int(n, 1)] <- TRUE

  w <- runif(n, 0.5, 2)

  df <- data.frame(Y_miss = Y, X = X, Z = Z, w = w)
  df$Y_miss[!R] <- NA_real_

  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  weights_full <- as.numeric(weights(des))

  spec <- el_prepare_inputs(
    formula = Y_miss ~ X | Z,
    data = des$variables,
    weights = weights_full,
    n_total = NULL,
    design_object = des
  )

  auxiliary_summary <- el_resolve_auxiliaries(
    aux_design_full = spec$aux_design_full,
    respondent_mask = spec$respondent_mask,
    auxiliary_means = NULL,
    weights_full = weights_full
  )

  sc <- validate_and_apply_nmar_scaling(
    standardize = FALSE,
    has_aux = ncol(auxiliary_summary$auxiliary_design) > 0,
    response_model_matrix_unscaled = spec$missingness_design,
    aux_matrix_unscaled = auxiliary_summary$auxiliary_design,
    mu_x_unscaled = auxiliary_summary$means,
    weights = spec$respondent_weights
  )

  Zs <- sc$response_model_matrix_scaled
  Xs <- sc$auxiliary_matrix_scaled
  mu_x <- sc$mu_x_scaled

  fam <- logit_family()
  n_resp_weighted <- sum(spec$respondent_weights)
  N_pop <- spec$N_pop

  eq_fun <- el_build_equation_system_survey(
    fam, Zs, Xs, spec$respondent_weights, N_pop, n_resp_weighted, mu_x
  )
  jac_fun <- el_build_jacobian_survey(
    fam, Zs, Xs, spec$respondent_weights, N_pop, n_resp_weighted, mu_x
  )

  W0 <- min(max(n_resp_weighted / N_pop, 1e-12), 1 - 1e-12)
  z0 <- stats::qlogis(W0)

  theta <- c(rep(0, ncol(Zs)), z0, 0.1, rep(0, ncol(Xs)))

  J_num <- numDeriv::jacobian(eq_fun, theta)
  J_ana <- jac_fun(theta)

  denom <- max(1e-8, norm(J_num, type = "F"))
  rel_diff <- norm(J_ana - J_num, type = "F") / denom
  expect_lt(rel_diff, 1e-4)
})

test_that("survey link Jacobian wrt z respects W clamping", {
  skip_if_not_installed("survey")

  set.seed(9876)
  n <- 80
  y <- rnorm(n)
  r <- rep(TRUE, n)
  r[sample.int(n, n %/% 2)] <- FALSE
  w <- runif(n, 0.5, 2)
  df <- data.frame(y_miss = ifelse(r, y, NA_real_), w = w)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  weights_full <- as.numeric(weights(des))

  spec <- el_prepare_inputs(
    formula = y_miss ~ 1,
    data = des$variables,
    weights = weights_full,
    n_total = NULL,
    design_object = des
  )

  Zs <- spec$missingness_design
  Xs <- matrix(nrow = nrow(Zs), ncol = 0)
  mu_x <- numeric(0)

  fam <- logit_family()
  n_resp_weighted <- sum(spec$respondent_weights)
  N_pop <- spec$N_pop

  eq_fun <- el_build_equation_system_survey(
    fam, Zs, Xs, spec$respondent_weights, N_pop, n_resp_weighted, mu_x
  )
  jac_fun <- el_build_jacobian_survey(
    fam, Zs, Xs, spec$respondent_weights, N_pop, n_resp_weighted, mu_x
  )

# Force W to be clamped at the upper bound; derivative wrt z must be 0.
  z_big <- 40
  theta <- c(rep(0, ncol(Zs)), z_big, 0.1)
  idx_z <- ncol(Zs) + 1L
  idx_link <- length(theta)

  J_ana <- jac_fun(theta)
  d_ana <- as.numeric(J_ana[idx_link, idx_z])

  eps <- 1e-4
  theta_p <- theta; theta_p[idx_z] <- theta_p[idx_z] + eps
  theta_m <- theta; theta_m[idx_z] <- theta_m[idx_z] - eps
  d_num <- (eq_fun(theta_p)[idx_link] - eq_fun(theta_m)[idx_link]) / (2 * eps)

  expect_equal(d_num, 0, tolerance = 1e-10)
  expect_equal(d_ana, d_num, tolerance = 1e-8)
})
