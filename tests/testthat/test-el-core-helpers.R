test_that("el_run_solver solves toy system (Newton with analytic Jacobian)", {
  skip_on_cran()
  f <- function(x) c(x[1] - 1, x[2] - 2, x[3])
  j <- function(x) diag(3)
  init <- c(0, 0, 0)
  ctrl <- list(maxit = 50)
  out <- NMAR:::el_run_solver(
    equation_system_func = f,
    analytical_jac_func = j,
    init = init,
    final_control = ctrl,
    top_args = list(),
    solver_method = "auto",
    use_solver_jac = TRUE,
    K_beta = 2, K_aux = 0,
    respondent_weights = c(1, 1, 1),
    N_pop = 3
  )
  expect_lte(out$solution$termcd, 2)
  expect_equal(out$solution$x, c(1, 2, 0), tolerance = 1e-8)
})

test_that("el_post_solution returns sane weights and mean in trivial case", {
  skip_on_cran()
  set.seed(123)
  n <- 5
  df <- data.frame(y = c(1, 0, 1, 0, 1), x = rnorm(n))
  X <- model.matrix(~x, data = df)
# No auxiliaries
  aux <- matrix(nrow = n, ncol = 0)
  mu_x <- NULL
  wts <- rep(1, n)
  fam <- NMAR:::logit_family()
  est <- c(0, 0, 0) # beta=(0,0), z=0 => W=0.5
  out <- NMAR:::el_post_solution(
    estimates = est,
    response_model_matrix_scaled = X,
    response_model_matrix_unscaled = X,
    auxiliary_matrix_scaled = aux,
    mu_x_scaled = mu_x,
    respondent_data = df,
    outcome_var = "y",
    family = fam,
    N_pop = sum(wts),
    respondent_weights = wts,
    K_beta = ncol(X),
    K_aux = 0,
    nmar_scaling_recipe = NULL,
    standardize = FALSE,
    trim_cap = Inf
  )
  expect_false(out$error)
# With lambda_W=0, denominators are 1, so weights equal base weights
  expect_equal(out$weights, wts)
  expect_equal(out$y_hat, mean(df$y))
})

## Removed: variance Jacobian selection; runtime is analytic-only

test_that("el_build_start applies user-supplied beta/W/lambda on the correct scale", {
  set.seed(1618)
  n <- 25
  df <- data.frame(y = rnorm(n), x = rnorm(n))
  Z_un <- model.matrix(~ y + x, data = df)
  X_un <- model.matrix(~ x - 1, data = df)
  aux_means <- c(x = mean(df$x))
  sc <- NMAR:::validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = TRUE,
    response_model_matrix_unscaled = Z_un,
    auxiliary_matrix_unscaled = X_un,
    mu_x_unscaled = aux_means,
    weights = rep(1, n)
  )
  start <- list(
    beta = c("(Intercept)" = 0.25, y = -0.15, x = 0.35),
    W = 0.65,
    lambda = c(x = 0.5)
  )
  res <- NMAR:::el_build_start(
    response_model_matrix_scaled = sc$response_model_matrix_scaled,
    auxiliary_matrix_scaled = sc$auxiliary_matrix_scaled,
    nmar_scaling_recipe = sc$nmar_scaling_recipe,
    start = start,
    N_pop = 2 * n,
    respondent_weights = rep(1, n)
  )
  expected_beta <- NMAR:::scale_coefficients(
    start$beta,
    sc$nmar_scaling_recipe,
    colnames(sc$response_model_matrix_scaled)
  )
  expected_lambda <- NMAR:::scale_aux_multipliers(
    start$lambda,
    sc$nmar_scaling_recipe,
    colnames(sc$auxiliary_matrix_scaled)
  )
  expect_equal(res$init_beta, expected_beta, tolerance = 1e-12)
  expect_equal(res$init_lambda, expected_lambda, tolerance = 1e-12)
  expect_equal(res$init_z, qlogis(start$W), tolerance = 1e-12)
})

test_that("el_build_start prefers start$z over start$W when both supplied", {
  set.seed(2718)
  n <- 10
  df <- data.frame(y = rnorm(n), x = rnorm(n))
  Z_un <- model.matrix(~ y + x, data = df)
  sc <- NMAR:::validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = FALSE,
    response_model_matrix_unscaled = Z_un,
    auxiliary_matrix_unscaled = matrix(nrow = n, ncol = 0),
    mu_x_unscaled = NULL,
    weights = rep(1, n)
  )
  start <- list(
    beta = c("(Intercept)" = 0.1, y = 0.2, x = -0.3),
    W = 0.8,
    z = 0.25
  )
  res <- NMAR:::el_build_start(
    response_model_matrix_scaled = sc$response_model_matrix_scaled,
    auxiliary_matrix_scaled = sc$auxiliary_matrix_scaled,
    nmar_scaling_recipe = sc$nmar_scaling_recipe,
    start = start,
    N_pop = 2 * n,
    respondent_weights = rep(1, n)
  )
  expect_equal(res$init_z, 0.25, tolerance = 1e-12)
})
