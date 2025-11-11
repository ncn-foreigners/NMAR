test_that("EL respects formula environment for custom transforms", {
  set.seed(202)
  n <- 60
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.4 * X + 0.1 * Z + rnorm(n)
  p <- plogis(-0.1 + 0.5 * scale(Y)[, 1] + 0.1 * X)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_

  f <- function(z) z^2 + 1
# Use transform in response-only RHS
  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
  res <- nmar(Y_miss ~ X | f(Z), data = df, engine = eng)
  expect_s3_class(res, "nmar_result_el")
})
