test_that("EL removes auxiliary intercept (+1) with a warning", {
  set.seed(101)
  n <- 80
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.5 * X + rnorm(n)
  p <- plogis(-0.3 + 0.4 * scale(Y)[, 1])
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_

  eng <- el_engine(variance_method = "none")

# Intercept in auxiliary part should be removed with a warning, not error
  expect_warning(
    {
      res1 <- nmar(Y_miss ~ X + 1, data = df, engine = eng, trace_level = 0)
      expect_s3_class(res1, "nmar_result_el")
    },
    regexp = "intercept",
    fixed = FALSE
  )

# Also fails when partitioned with response-only predictors
  expect_warning(
    {
      res2 <- nmar(Y_miss ~ X + 1 | Z, data = df, engine = eng, trace_level = 0)
      expect_s3_class(res2, "nmar_result_el")
    },
    regexp = "intercept",
    fixed = FALSE
  )
})
