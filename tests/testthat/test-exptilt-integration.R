## Optional integration smoke tests: skipped on CI by default because each run
## spins up a full EM fit. Set NMAR_RUN_INTEGRATION=1 locally to exercise them
test_that("exptilt converges on simple IID example", {
  skip_if(Sys.getenv("NMAR_RUN_INTEGRATION") != "1", "Set NMAR_RUN_INTEGRATION=1 to run integration tests.")
  set.seed(123)
  # Use a moderately sized synthetic sample so the analytic delta variance
  # has a stable Fisher information (avoids singular FI22 in tiny samples)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y_true <- 0.5 + 0.3 * x1 - 0.2 * x2 + rnorm(n, sd = 0.2)
  respond <- rbinom(n, 1, plogis(2 + 0.1 * y_true))
  if (all(respond == 1)) respond[sample.int(n, 1)] <- 0  # ensure at least one missing
  y_obs <- ifelse(respond == 1, y_true, NA)
  dat <- data.frame(y = y_obs, x1 = x1, x2 = x2)

  fit <- nmar::nmar(
    y ~ x1 + x2,
    data = dat,
    engine = nmar::exptilt_engine(
      y_dens = "normal",
      variance_method = "delta",
      min_iter = 1,
      max_iter = 8,
      tol_value = 1e-3
    )
  )
  expect_true(fit$converged)
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$std_error))

  # The fit surfaces successfully and produces finite estimates
  expect_true(is.finite(fit$std_error))
})

test_that("exptilt handles simple survey design", {
  skip_if(Sys.getenv("NMAR_RUN_INTEGRATION") != "1", "Set NMAR_RUN_INTEGRATION=1 to run integration tests.")
  skip_if_not_installed("survey")
  set.seed(321)
  # Larger survey sample keeps replicate fits numerically stable during
  # bootstrap variance calculations
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  w <- runif(n, 0.5, 2)
  y_true <- 1 + 0.4 * x1 + rnorm(n, sd = 0.3)
  respond <- rbinom(n, 1, plogis(2 + 0.1 * y_true))
  if (all(respond == 1)) respond[sample.int(n, 1)] <- 0
  y_obs <- ifelse(respond == 1, y_true, NA)
  dat <- data.frame(y = y_obs, x1 = x1, x2 = x2, w = w)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = dat)

  fit <- nmar::nmar(
    y ~ x1 + x2,
    data = des,
    engine = nmar::exptilt_engine(
      y_dens = "normal",
      variance_method = "bootstrap",
      bootstrap_reps = 10,
      min_iter = 1,
      max_iter = 8,
      tol_value = 1e-3
    )
  )
  expect_true(fit$converged)
  expect_true(is.finite(fit$estimate))
  expect_equal(fit$sample$n_total, n)
})
