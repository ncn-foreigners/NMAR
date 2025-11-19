test_that("on_failure=return: inconsistent aux means yield a warning and diagnostics, not hard failure", {
  set.seed(3301)
  N <- 200
  X1 <- rnorm(N); X2 <- rnorm(N)
  Y <- 1 + 0.4 * X1 + rnorm(N)
  p <- plogis(-0.5 + 0.6 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X1 = X1, X2 = X2)
  df[!R, "Y_miss"] <- NA_real_
  bad_aux <- c(X1 = 10, X2 = -10)
  expect_warning(
    fit <- el.data.frame(df, Y_miss ~ X1 + X2,
                                 auxiliary_means = bad_aux, on_failure = "return", variance_method = "delta"),
    regexp = "Auxiliary means appear far from respondents' support"
  )
# Soft diagnostics fields must exist
  expect_true("aux_inconsistency_max_z" %in% names(fit$diagnostics))
  expect_true("aux_inconsistency_cols" %in% names(fit$diagnostics))
})

test_that("trimming caps weights and sets trimmed_fraction > 0", {
  set.seed(3302)
  df <- make_iid_nmar(n = 300, alpha = 0.8, seed = 3302)
  fit <- nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = make_engine(auxiliary_means = c(X = 0), trim_cap = 2, variance_method = "bootstrap", bootstrap_reps = 10)
  )
  w <- weights(fit)
  expect_true(max(w) <= 2 + 1e-8)
  expect_true(attr(w, "trimmed_fraction") > 0)
})
