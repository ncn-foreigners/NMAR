test_that("outcome transformations propagate through estimation", {
  set.seed(741)
  n <- 150
  X <- rnorm(n)
  Y_true <- exp(0.5 + 0.4 * X + rnorm(n, sd = 0.2))
  p <- plogis(-0.2 + 0.3 * scale(Y_true)[, 1] + 0.5 * X)
  respond <- runif(n) < p
  Y_miss <- Y_true
  Y_miss[!respond] <- NA_real_
  df <- data.frame(Y_miss = Y_miss, X = X)

  fit <- nmar(
    log(Y_miss) ~ X,
    data = df,
    engine = el_engine(standardize = FALSE, variance_method = "none", trim_cap = Inf)
  )
  expect_true(fit$converged)
  expect_identical(fit$estimate_name, "log(Y_miss)")
  obs <- !is.na(df$Y_miss)
  w <- weights(fit)
  expect_length(w, sum(obs))
  manual <- sum(w * log(df$Y_miss[obs]))
  expect_equal(manual, fit$y_hat, tolerance = 1e-8)
})

test_that("invalid outcome transforms are rejected", {
  set.seed(123)
  df <- data.frame(Y_miss = c(1, 2, NA, 4), X = 1:4)
  bad_engine <- el_engine(standardize = FALSE, variance_method = "none", trim_cap = Inf)
  expect_error(
    nmar(log(Y_miss - 5) ~ X, data = df, engine = bad_engine),
    "non-finite"
  )
})
