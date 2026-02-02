test_that("el_run_solver solves toy system (Newton with analytic Jacobian)", {
  skip_on_cran()
  f <- function(x) c(x[1] - 1, x[2] - 2, x[3])
  j <- function(x) diag(3)
  init <- c(0, 0, 0)
  ctrl <- list(maxit = 50)
  out <- el_run_solver(
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
  fam <- logit_family()
  est <- c(0, 0, 0) # beta=(0,0), z=0 => W=0.5
  out <- el_post_solution(
    estimates = est,
    missingness_model_matrix_scaled = X,
    missingness_model_matrix_unscaled = X,
    auxiliary_matrix_scaled = aux,
    mu_x_scaled = mu_x,
    response_outcome = df$y,
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
