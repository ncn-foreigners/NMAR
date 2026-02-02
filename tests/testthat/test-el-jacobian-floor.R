test_that("EL Jacobian respects denominator floor (active mask)", {
  skip_on_cran()
  set.seed(42)

# Setup: choose parameters so many Di are below the floor
  n <- 60
  X <- cbind(`(Intercept)` = rep(1, n))
  aux <- matrix(nrow = n, ncol = 0)
  wts <- rep(1, n)
  N_pop <- ceiling(1.2 * n) # C_const > 0
  fam <- logit_family()
  n_resp_wt <- sum(wts)
  mu_x <- numeric(0)

  eq <- el_build_equation_system(
    family = fam,
    missingness_model_matrix = X,
    auxiliary_matrix = aux,
    respondent_weights = wts,
    N_pop = N_pop,
    n_resp_weighted = n_resp_wt,
    mu_x_scaled = mu_x
  )
  jac <- el_build_jacobian(
    family = fam,
    missingness_model_matrix = X,
    auxiliary_matrix = aux,
    respondent_weights = wts,
    N_pop = N_pop,
    n_resp_weighted = n_resp_wt,
    mu_x_scaled = mu_x
  )

# Pick W very close to 1 so lambda_W is large; pick beta << 0 so p_i << W
  beta <- -12
  W <- 0.99
  z <- qlogis(W)
  theta <- c(beta, z)

# Sanity: verify that many denominators would be below floor without clamping
# Reconstruct denominator components mirroring engine logic
  C_const <- (N_pop / n_resp_wt) - 1
  eta <- as.numeric(X %*% beta)
  p <- fam$linkinv(eta)
  lambda_W <- C_const / (1 - W)
  denom_raw <- 1 + lambda_W * (p - W) # no auxiliaries
  floor_val <- getOption("nmar.el_denom_floor", 1e-8)
  expect_true(mean(denom_raw <= floor_val) > 0.2)

# Analytic vs numeric Jacobian should still agree tightly thanks to 'active'
  fn <- function(t) as.numeric(eq(t))
  J_num <- numDeriv::jacobian(fn, x = theta)
  J_ana <- jac(theta)
  expect_equal(unname(J_ana), unname(J_num), tolerance = 1e-6)
})
