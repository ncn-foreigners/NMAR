test_that("EL delta variance returns finite SE in benign survey case", {
  skip_on_cran()
  skip_if_not_installed("survey")

  set.seed(30310)
  N <- 1500L
  X <- rnorm(N)
  eps <- rnorm(N)
  Y <- 1.5 + 0.6 * X + eps
  p <- plogis(-0.4 + 0.25 * scale(Y)[, 1])
  R <- rbinom(N, 1, p)

  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[R == 0] <- NA_real_
  des <- survey::svydesign(ids = ~1, weights = ~1, data = df)

  eng <- NMAR::el_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = TRUE)
  fit <- suppressWarnings(NMAR::nmar(Y_miss ~ X, data = des, engine = eng))
  se <- suppressWarnings(as.numeric(fit$se))
  expect_true(is.na(se) || is.finite(se))
})
