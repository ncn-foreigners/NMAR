test_that("EL delta variance calibrates in a benign IID regime (integration)", {
  skip_on_cran()
  if (Sys.getenv("NMAR_RUN_INTEGRATION", unset = "0") != "1") skip("integration-only")

  set.seed(20210)
  N <- 2000L
  X <- rnorm(N)
  Z <- rnorm(N)
  eps <- rnorm(N)
  Y <- 2 + 0.5 * X + eps
  mu_true <- mean(Y)
  p <- plogis(-0.4 + 0.3 * scale(Y)[, 1])

  B <- 100L
  est <- se_d <- rep(NA_real_, B)
  eng <- NMAR::el_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = TRUE)

  for (b in seq_len(B)) {
    set.seed(20210 + b)
    R <- rbinom(N, 1, p)
    df <- data.frame(Y_miss = Y, X = X)
    df$Y_miss[R == 0] <- NA_real_
    fit <- NMAR::nmar(Y_miss ~ X, data = df, engine = eng)
    est[b] <- as.numeric(fit$y_hat)
    se_d[b] <- suppressWarnings(as.numeric(fit$se))
  }
# Basic sanity: delta SEs mostly finite and non-trivial
  expect_gt(mean(is.finite(se_d)), 0.8)
  expect_gt(median(se_d, na.rm = TRUE), 1e-3)

  var_emp <- stats::var(est, na.rm = TRUE)
  var_bar <- mean(se_d^2, na.rm = TRUE)
  ratio <- var_emp / var_bar
# Rough calibration gate (first-order delta; allow generous bounds)
  expect_true(is.finite(ratio))
  expect_gt(ratio, 0.3)
  expect_lt(ratio, 3.0)
})
