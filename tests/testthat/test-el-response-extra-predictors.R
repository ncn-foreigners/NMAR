test_that("response_predictors can include non-auxiliary variables", {
  set.seed(4242)
  N <- 350
  X <- rnorm(N)
  Z <- rnorm(N)
  Y <- 1 + 0.7 * X + 0.3 * Z + rnorm(N)
  p <- plogis(-0.4 + 0.6 * scale(Y)[, 1] + 0.3 * Z)
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_

# Outcome RHS contains only X (auxiliary). Include Z as response-only predictor.
  res <- nmar:::el.data.frame(
    data = df,
    formula = Y_miss ~ X,
    response_predictors = c("Z"),
    auxiliary_means = c(X = 0),
    standardize = FALSE,
    trim_cap = Inf,
    on_failure = "return",
    variance_method = "delta"
  )

  expect_s3_class(res, "nmar_result_el")
  expect_true(isTRUE(res$converged))
  expect_true(is.numeric(res[['std_error']]))
  expect_true(is.na(res[['std_error']]) || res[['std_error']] >= 0)
# Coefficient vector should include Z as a response predictor
  expect_true(any(grepl("(^|\\b)Z(\\b|$)", names(res$model$coefficients))))
})
