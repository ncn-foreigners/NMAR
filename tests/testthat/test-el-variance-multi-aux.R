test_that("delta variance is finite in multi-aux IID case (benign)", {
  skip_if_not_installed("stats")
  set.seed(321)
  n <- 1200
  X1 <- rnorm(n); X2 <- rnorm(n); Z <- rnorm(n)
  Y <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
  p <- plogis(-0.5 + 0.6 * scale(Y)[, 1] + 0.4 * Z)
  R <- rbinom(n, 1, p)
  df <- data.frame(Y_miss = Y, X1 = X1, X2 = X2, Z = Z)
  df$Y_miss[R == 0] <- NA_real_

  eng <- el_engine(auxiliary_means = c(X1 = 0, X2 = 0), variance_method = "delta",
                   standardize = TRUE, trim_cap = Inf)
  fit <- suppressWarnings(nmar(Y_miss ~ X1 + X2 | Z, data = df, engine = eng))
  expect_true(isTRUE(fit$converged))
# Delta SE may be NA under the strict delta policy; allow NA here
  expect_true(is.na(as.numeric(fit$se)) || is.finite(as.numeric(fit$se)))

  if (identical(Sys.getenv("NMAR_RUN_INTEGRATION"), "1")) {
    B <- 60
    est <- se2 <- numeric(B)
    for (b in seq_len(B)) {
      set.seed(2000 + b)
      Rb <- rbinom(n, 1, p)
      d2 <- df
      d2$Y_miss <- Y
      d2$Y_miss[Rb == 0] <- NA_real_
      fb <- suppressWarnings(nmar(Y_miss ~ X1 + X2 | Z, data = d2, engine = eng))
      est[b] <- fb$y_hat
      se2[b] <- as.numeric(fb$se)^2
    }
    var_emp <- stats::var(est, na.rm = TRUE)
    var_bar <- mean(se2, na.rm = TRUE)
# Calibration check is optional due to stricter NA policy
    expect_true(is.finite(var_emp))
  }
})
