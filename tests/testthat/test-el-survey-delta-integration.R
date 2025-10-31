test_that("Survey delta variance is finite in a benign SRS regime (integration)", {
  skip("EL delta variance disabled; use bootstrap for SEs")
  skip_on_cran()
  skip_if_not_installed("survey")
  if (Sys.getenv("NMAR_RUN_INTEGRATION", unset = "0") != "1") skip("integration-only")

  set.seed(7004)
  N <- 3000
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-0.4 + 0.3 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[!R] <- NA_real_
  des <- survey::svydesign(ids = ~1, weights = ~1, data = df)
  eng <- el_engine(auxiliary_means = c(X = mean(X)), variance_method = "delta", standardize = TRUE)
  fit <- nmar(Y_miss ~ X, data = des, engine = eng)
  expect_true(isTRUE(fit$converged))
  se <- as.numeric(fit$se)
  expect_true(is.finite(se))
})
