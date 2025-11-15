# Simple Jacobian check for EL equations with unified denominator floor

test_that("EL analytic Jacobian matches numeric (no aux, logit)", {
  skip_on_cran()
  set.seed(123)

# Small synthetic problem
  n <- 50
  y <- rnorm(n)
  X <- cbind(`(Intercept)` = 1, x1 = rnorm(n))

# Build objects the same way as the engine would for respondents only
  fam <- NMAR:::logit_family()
  resp_w <- rep(1, n)
  N_pop <- ceiling(1.2 * n) # ensure lambda_W > 0

  eq <- NMAR:::el_build_equation_system(
    family = fam,
    missingness_model_matrix = X,
    auxiliary_matrix = matrix(nrow = n, ncol = 0),
    respondent_weights = resp_w,
    N_pop = N_pop,
    n_resp_weighted = sum(resp_w),
    mu_x_scaled = numeric(0)
  )
  jac <- NMAR:::el_build_jacobian(
    family = fam,
    missingness_model_matrix = X,
    auxiliary_matrix = matrix(nrow = n, ncol = 0),
    respondent_weights = resp_w,
    N_pop = N_pop,
    n_resp_weighted = sum(resp_w),
    mu_x_scaled = numeric(0)
  )

# Random but reasonable parameters: beta near 0, z near logit(0.6)
  beta <- c(0.1, -0.2)
  z <- qlogis(0.6)
  theta <- c(beta, z)

# Numeric Jacobian via numDeriv
  fn <- function(t) as.numeric(eq(t))
  J_num <- numDeriv::jacobian(fn, x = theta)
  J_an <- jac(theta)

# Tolerance generous but should be tight when denominators are well above floor
# Compare numerically, ignoring dimnames
  expect_equal(unname(J_an), unname(J_num), tolerance = 1e-6)
})
