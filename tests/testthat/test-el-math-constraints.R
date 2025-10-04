test_that("estimating equations solved to tight tolerance (max residual small)", {
  set.seed(3402)
  N <- 300
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-1 + 0.5 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_

  fit <- NMAR:::el.data.frame(df, Y_miss ~ X,
    response_predictors = NULL,
    auxiliary_means = c(X = 0), standardize = FALSE,
    trim_cap = Inf, variance_method = "delta"
  )
  expect_true(fit$converged)
  expect_lt(fit$diagnostics$max_equation_residual, 1e-6)
  expect_true(is.finite(fit$diagnostics$min_denominator))
})
