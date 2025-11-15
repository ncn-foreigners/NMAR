test_that("EL errors early when auxiliaries contain NA (before solver)", {
  set.seed(123)
# Build a small dataset with NA in an auxiliary variable X
  n <- 50
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.5 * X + 0.3 * Z + rnorm(n)
  p <- plogis(-0.3 + 0.4 * scale(Y)[, 1] + 0.2 * Z)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_
# Inject a single NA in the auxiliary variable among respondents
  idx <- which(!is.na(df$Y_miss))[1]
  df$X[idx] <- NA_real_

  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")

# Expect a clear error about NA in the auxiliary covariate and not a row-mismatch error
  expect_error(
    nmar(Y_miss ~ X, data = df, engine = eng),
    regexp = "contains NA values",
    info = "Auxiliary NA should be caught before solver setup"
  )
})
