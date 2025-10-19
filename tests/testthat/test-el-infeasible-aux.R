test_that("Infeasible/inconsistent auxiliaries trigger a warning and diagnostics (soft check)", {
  set.seed(7002)
  N <- 800
  X <- rnorm(N)
  Y <- 1 + 0.5 * X + rnorm(N)
  p <- plogis(-0.4 + 0.3 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[!R] <- NA_real_
# Supply a wildly inconsistent auxiliary mean to trigger infeasibility guard
  eng <- el_engine(auxiliary_means = c(X = 10), variance_method = "none", on_failure = "return")
  expect_warning(
    fit <- nmar(Y_miss ~ X, data = df, engine = eng),
    regexp = "Auxiliary means appear far from respondents' support"
  )
# Soft diagnostics fields must exist (values may be NA in edge draws)
  expect_true("aux_inconsistency_max_z" %in% names(fit$diagnostics))
  expect_true("aux_inconsistency_cols" %in% names(fit$diagnostics))
})
