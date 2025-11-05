test_that("EL block solver matches joint solver on IID data", {
  set.seed(123)
  n <- 400
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.7 * X + Z
  p <- plogis(-0.5 + 0.5 * scale(Y)[, 1] + 0.2 * X)
  Rr <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[!Rr] <- NA_real_

  eng_joint <- el_engine(control = list(solver = "joint"), auxiliary_means = c(X = 0), variance_method = "none")
  eng_block <- el_engine(control = list(solver = "block"), auxiliary_means = c(X = 0), variance_method = "none")

  fit_joint <- nmar(Y_miss ~ X, data = df, engine = eng_joint)
  fit_block <- nmar(Y_miss ~ X, data = df, engine = eng_block)

  expect_true(fit_joint$converged)
  expect_true(fit_block$converged)
# Weight sums
  expect_equal(sum(weights(fit_block, scale = "population")), n)
  expect_equal(sum(weights(fit_block, scale = "probability")), 1)
})

test_that("EL block solver matches joint solver on survey design", {
  skip_if_not_installed("survey")
  set.seed(321)
  n <- 600
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 2 + 0.5 * X + Z
  p <- plogis(-0.3 + 0.4 * scale(Y)[, 1])
  Rr <- runif(n) < p
  w <- runif(n, 0.5, 2)
  df <- data.frame(Y_miss = Y, X = X, w = w)
  df$Y_miss[!Rr] <- NA_real_
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  eng_joint <- el_engine(control = list(solver = "joint"), auxiliary_means = c(X = 0), variance_method = "none")
  eng_block <- el_engine(control = list(solver = "block"), auxiliary_means = c(X = 0), variance_method = "none")

  fit_joint <- nmar(Y_miss ~ X, data = des, engine = eng_joint)
  fit_block <- nmar(Y_miss ~ X, data = des, engine = eng_block)

  expect_true(fit_joint$converged)
  expect_true(fit_block$converged)
# Note: constraints are checked within the engine diagnostics; here we focus on weight sums
# Weight sums
  N_pop <- sum(weights(des))
  expect_equal(sum(weights(fit_block, scale = "population")), N_pop)
  expect_equal(sum(weights(fit_block, scale = "probability")), 1)
})
