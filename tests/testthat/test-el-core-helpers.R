test_that("el_run_solver solves toy system with chosen method", {
  skip_on_cran()
# Toy system: root at (1,2,0)
  f <- function(x) c(x[1] - 1, x[2] - 2, x[3])
  j <- function(x) diag(3)
  init <- c(0, 0, 0)
  ctrl <- list(maxit = 50)
# Newton
  out_n <- NMAR:::el_run_solver(
    equation_system_func = f,
    analytical_jac_func = j,
    init = init,
    final_control = ctrl,
    top_args = list(),
    solver_method = "newton",
    use_solver_jac = TRUE,
    K_beta = 2, K_aux = 0,
    respondent_weights = c(1, 1, 1),
    N_pop = 3
  )
  expect_equal(out_n$method, "Newton")
  expect_lte(out_n$solution$termcd, 2)
  expect_equal(out_n$solution$x, c(1, 2, 0), tolerance = 1e-8)
# Broyden
  out_b <- NMAR:::el_run_solver(
    equation_system_func = f,
    analytical_jac_func = j,
    init = init,
    final_control = ctrl,
    top_args = list(),
    solver_method = "broyden",
    use_solver_jac = FALSE,
    K_beta = 2, K_aux = 0,
    respondent_weights = c(1, 1, 1),
    N_pop = 3
  )
  expect_equal(out_b$method, "Broyden")
  expect_lte(out_b$solution$termcd, 2)
  expect_equal(out_b$solution$x, c(1, 2, 0), tolerance = 1e-6)
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

test_that("el_select_variance_jacobian chooses numeric when analytic differs", {
  skip_on_cran()
# Toy system: zero at (1,2,0)
  f <- function(x) c(x[1] - 1, x[2] - 2, x[3])
  j_good <- function(x) diag(3)
  j_bad <- function(x) diag(3) + 0.1 # intentionally biased to trigger rel_diff gate
  at <- c(1, 2, 0)
  sel_auto <- NMAR:::el_select_variance_jacobian(
    equation_system_func = f,
    analytical_jac_func = j_bad,
    estimates = at,
    variance_jacobian = "auto"
  )
  expect_equal(sel_auto$A_source, "numeric")
  expect_equal(sel_auto$jacobian_auto_rule, "rel_diff_high")
# Forcing analytic should return analytic regardless
  sel_an <- NMAR:::el_select_variance_jacobian(
    equation_system_func = f,
    analytical_jac_func = j_good,
    estimates = at,
    variance_jacobian = "analytic"
  )
  expect_equal(sel_an$A_source, "analytic")
})
