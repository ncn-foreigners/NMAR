test_that("nleqslv top-level control (global/xscalm) is passed and recorded", {
  set.seed(101)
  n <- 400
  X <- rnorm(n)
  Y <- 1 + 0.4 * X + rnorm(n)
  p <- plogis(-0.5 + 0.3 * scale(Y)[, 1])
  R <- rbinom(n, 1, p)
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[R == 0] <- NA_real_
  eng <- NMAR::el_engine(auxiliary_means = c(X = 0), variance_method = 'none', control = list(global = "none", xscalm = "fixed", maxit = 5))
  fit <- NMAR::nmar(Y_miss ~ X, data = df, engine = eng)
  di <- fit$diagnostics
  expect_true(fit$converged)
  expect_identical(di$nleqslv_global, "none")
  expect_identical(di$nleqslv_xscalm, "fixed")
  expect_true(is.finite(di$solver_iterations))
  expect_lte(di$solver_iterations, 5)
})

test_that("invalid global/xscalm are coerced to defaults", {
  set.seed(102)
  n <- 200
  Y <- rnorm(n)
  R <- rbinom(n, 1, plogis(0))
  df <- data.frame(Y_miss = Y)
  df$Y_miss[R == 0] <- NA_real_
  eng <- NMAR::el_engine(variance_method = 'none', control = list(global = "badvalue", xscalm = "bad"))
  fit <- NMAR::nmar(Y_miss ~ 1, data = df, engine = eng)
  di <- fit$diagnostics
# Defaults when invalid: qline, auto
  expect_identical(di$nleqslv_global, "qline")
  expect_identical(di$nleqslv_xscalm, "auto")
})
