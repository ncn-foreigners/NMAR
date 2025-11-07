test_that("EL probit equations/Jacobian are stable in extreme tails", {
  skip_on_cran()
  set.seed(123)

# Construct a small respondents-only problem with probit link
  n <- 40
# Intercept-only design to force constant extreme eta
  X <- cbind(`(Intercept)` = rep(1, n))
  aux <- matrix(nrow = n, ncol = 0)
  wts <- rep(1, n)
  N_pop <- ceiling(1.2 * n) # ensure lambda_W > 0

  fam <- NMAR:::probit_family()
  n_resp_wt <- sum(wts)
  mu_x <- numeric(0)

  eq <- NMAR:::el_build_equation_system(
    family = fam,
    response_model_matrix = X,
    auxiliary_matrix = aux,
    respondent_weights = wts,
    N_pop = N_pop,
    n_resp_weighted = n_resp_wt,
    mu_x_scaled = mu_x
  )
  jac <- NMAR:::el_build_jacobian(
    family = fam,
    response_model_matrix = X,
    auxiliary_matrix = aux,
    respondent_weights = wts,
    N_pop = N_pop,
    n_resp_weighted = n_resp_wt,
    mu_x_scaled = mu_x
  )

# Parameters: very negative eta for all units; moderate W
  beta <- -12 # extreme tail
  z <- qlogis(0.6)
  theta <- c(beta, z)

# Equations finite and stable
  vals <- as.numeric(eq(theta))
  expect_true(all(is.finite(vals)))

# Analytic vs numeric Jacobian should agree tightly
  fn <- function(t) as.numeric(eq(t))
  J_num <- numDeriv::jacobian(fn, x = theta)
  J_ana <- jac(theta)
  expect_equal(unname(J_ana), unname(J_num), tolerance = 1e-6)
})
