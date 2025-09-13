test_that("on_failure=return yields converged=FALSE on inconsistent aux means", {
  set.seed(3301)
  N <- 200
  X1 <- rnorm(N)
  X2 <- rnorm(N)
  Y <- 1 + 0.4 * X1 + rnorm(N)
  p <- plogis(-0.5 + 0.6 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X1 = X1, X2 = X2)
  df[!R, "Y_miss"] <- NA_real_
  bad_aux <- c(X1 = 10, X2 = -10)
  fit <- nmar:::el.data.frame(df, Y_miss ~ X1 + X2, response_predictors = NULL, auxiliary_means = bad_aux, on_failure = "return", variance_method = "delta")
  expect_false(fit$converged)
  expect_true(is.na(fit$y_hat))
})

test_that("trimming caps weights and sets trimmed_fraction > 0", {
  set.seed(3302)
  N <- 300
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-1 + 0.8 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_
  fit <- nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = el_engine(auxiliary_means = c(X = 0), trim_cap = 2, variance_method = "bootstrap", bootstrap_reps = 10, suppress_warnings = TRUE)
  )
  w <- weights(fit)
  expect_true(max(w) <= 2 + 1e-8)
  expect_true(attr(w, "trimmed_fraction") > 0)
})
