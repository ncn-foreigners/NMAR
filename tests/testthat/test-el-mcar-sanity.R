test_that("MCAR sanity: EL equals respondent mean (large N)", {
  set.seed(7001)
  N <- 3000
  X <- rnorm(N)
  Y <- 1.5 + 0.6 * X + rnorm(N)
  p <- plogis(-0.5) # constant response probability, independent of Y
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[!R] <- NA_real_
  eng <- el_engine(variance_method = "none", standardize = FALSE)
  fit <- nmar(Y_miss ~ 1, data = df, engine = eng)
  mu_resp <- mean(Y[R], na.rm = TRUE)
  expect_true(fit$converged)
  expect_lt(abs(as.numeric(fit$estimate) - mu_resp), 5e-3)
})
