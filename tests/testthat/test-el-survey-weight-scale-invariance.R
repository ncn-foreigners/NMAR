test_that("EL survey estimator is invariant to uniform rescaling of design weights", {
  skip_on_cran()
  skip_if_not_installed("survey")

  set.seed(2025)
  n <- 200
  X <- rnorm(n)
  Y <- 2 + 0.7 * X + rnorm(n)
# Missingness depends on Y (NMAR) to exercise the full machinery
  p <- plogis(-0.8 + 0.5 * scale(Y)[, 1])
  R <- rbinom(n, 1, p)
  Y_miss <- Y
  Y_miss[R == 0] <- NA_real_

# Base weights and a scaled copy
  w <- runif(n, 0.5, 2)
  cfac <- 7.3
  w2 <- cfac * w

  df1 <- data.frame(Y_miss = Y_miss, X = X, w = w)
  df2 <- data.frame(Y_miss = Y_miss, X = X, w2 = w2)

  des1 <- survey::svydesign(ids = ~1, weights = ~w, data = df1)
  des2 <- survey::svydesign(ids = ~1, weights = ~w2, data = df2)

  eng <- el_engine(variance_method = "none")

  fit1 <- nmar(Y_miss ~ X, data = des1, engine = eng)
  fit2 <- nmar(Y_miss ~ X, data = des2, engine = eng)

# Point estimate invariance (QLS Eq. 11 with weighted generalization)
  expect_equal(as.numeric(fit1$y_hat), as.numeric(fit2$y_hat), tolerance = 1e-8)

# Probability-scale weights should be identical
  p1 <- as.numeric(weights(fit1, scale = "probability"))
  p2 <- as.numeric(weights(fit2, scale = "probability"))
  expect_equal(p1, p2, tolerance = 1e-10)

# Population-scale weights should differ by the same cfac
  Wpop1 <- as.numeric(weights(fit1, scale = "population"))
  Wpop2 <- as.numeric(weights(fit2, scale = "population"))
# Compute realized scaling from reported population totals
  N1 <- fit1$sample$n_total
  N2 <- fit2$sample$n_total
  expect_true(is.finite(N1) && is.finite(N2) && N1 > 0 && N2 > 0)
  c_est <- as.numeric(N2 / N1)
  expect_equal(c_est, cfac, tolerance = 1e-12)
  expect_equal(Wpop2, c_est * Wpop1, tolerance = 1e-10)

# Constraint residuals scaled properly and negligible relative to weight scale
  diag1 <- fit1$diagnostics
  diag2 <- fit2$diagnostics
  sumw1 <- diag1$sum_respondent_weights
  sumw2 <- diag2$sum_respondent_weights
# Normalize by respondent-weight sum to remove scale
  norm_eqW1 <- abs(diag1$constraint_sum_W) / sumw1
  norm_eqW2 <- abs(diag2$constraint_sum_W) / sumw2
  expect_equal(norm_eqW1, norm_eqW2, tolerance = 1e-8)
  expect_lt(norm_eqW1, 1e-8)
  if (!is.null(diag1$constraint_sum_aux) && length(diag1$constraint_sum_aux) > 0) {
    a1 <- abs(diag1$constraint_sum_aux) / sumw1
    a2 <- abs(diag2$constraint_sum_aux) / sumw2
    expect_equal(a1, a2, tolerance = 1e-8)
    expect_true(all(a1 < 1e-8))
  }
})
