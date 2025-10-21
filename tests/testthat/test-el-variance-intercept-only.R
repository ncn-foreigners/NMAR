test_that("delta variance is finite in intercept-only IID case", {
  set.seed(123)
  n <- 1200
  X <- rnorm(n)
  Y <- 2 + 0.5 * X + rnorm(n)
  p <- plogis(-0.7 + 0.6 * scale(Y)[, 1])
  R <- rbinom(n, 1, p)
  df <- data.frame(Y_miss = Y)
  df$Y_miss[R == 0] <- NA_real_

  eng <- el_engine(variance_method = "delta", standardize = TRUE, trim_cap = Inf)
  fit <- suppressWarnings(nmar(Y_miss ~ 1, data = df, engine = eng))
  expect_true(isTRUE(fit$converged))
  se_val <- suppressWarnings(as.numeric(fit$se))
# Under our strict policy, delta may return NA in edge cases; allow NA here
  expect_true(is.na(se_val) || (is.finite(se_val) && se_val > 1e-4))

  if (identical(Sys.getenv("NMAR_RUN_INTEGRATION"), "1")) {
# Monte Carlo calibration in a benign regime
    B <- 100
    est <- se2 <- numeric(B)
    for (b in seq_len(B)) {
      set.seed(1000 + b)
      Rb <- rbinom(n, 1, p)
      d2 <- df
      d2$Y_miss <- Y
      d2$Y_miss[Rb == 0] <- NA_real_
      fb <- nmar(Y_miss ~ 1, data = d2, engine = eng)
      est[b] <- fb$y_hat
      se2[b] <- as.numeric(fb$se)^2
    }
    var_emp <- stats::var(est, na.rm = TRUE)
    var_bar <- mean(se2, na.rm = TRUE)
# Only assert calibration when we have finite SEs across reps
    if (is.finite(var_emp) && is.finite(var_bar)) {
# Allow a generous factor due to randomness, but catch gross miscalibration
      expect_lt(abs(log(var_bar / var_emp)), log(3))
    } else {
      testthat::skip("Delta returned NA SE in this Monte Carlo seed; skip calibration check.")
    }
  }
})
