test_that("EL supports response-only with no auxiliaries (y ~ 1 | z)", {
  set.seed(1001)
  n <- 120
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.6 * X + 0.3 * Z + rnorm(n)
  p <- plogis(-0.3 + 0.5 * scale(Y)[, 1] + 0.4 * Z)
  R <- rbinom(n, 1, p)
  df <- data.frame(Y_miss = Y, Z = Z)
  df$Y_miss[R == 0] <- NA_real_
  eng <- el_engine(auxiliary_means = NULL, variance_method = "none")
  fit <- nmar(Y_miss ~ 1 | Z, data = df, engine = eng)
  expect_s3_class(fit, "nmar_result_el")
  expect_true(isTRUE(fit$converged))
  expect_true(is.finite(as.numeric(fit$y_hat)) || is.na(as.numeric(fit$y_hat)))
})

test_that("EL allows overlap between aux and response sides", {
  set.seed(1002)
  n <- 100
  X <- rnorm(n)
  Y <- 2 + 0.5 * X + rnorm(n)
  p <- plogis(-0.4 + 0.5 * scale(Y)[, 1])
  R <- rbinom(n, 1, p)
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[R == 0] <- NA_real_
  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
  fit <- nmar(Y_miss ~ X | X, data = df, engine = eng)
  expect_s3_class(fit, "nmar_result_el")
  expect_true(isTRUE(fit$converged))
})

# test_that("ET rejects overlap and Y on RHS of response partition", {
#   set.seed(1003)
#   n <- 120
#   X <- rnorm(n); Z <- rnorm(n)
#   Y <- 1.2 + 0.4 * X + 0.2 * Z + rnorm(n)
#   p <- plogis(-0.2 + 0.6 * scale(Y)[, 1])
#   R <- rbinom(n, 1, p)
#   df <- data.frame(Y_miss = Y, X = X, Z = Z)
#   df$Y_miss[R == 0] <- NA_real_
#   eng <- exptilt_engine(y_dens = "normal", variance_method = "delta", standardize = FALSE)
# # Y on RHS should be rejected by ET (outcome in missingness not allowed by traits)
#   expect_error(nmar(Y_miss ~ X | Y_miss, data = df, engine = eng))
# # Overlap X on both sides should also be rejected for ET
#   expect_error(nmar(Y_miss ~ X | X, data = df, engine = eng))
# })

test_that("Engines accept parentheses and transforms in RHS partitions", {
  set.seed(1004)
  n <- 80
  X1 <- rnorm(n); X2 <- rnorm(n); Z1 <- rnorm(n)
  Y <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
  p <- plogis(-0.5 + 0.5 * scale(Y)[, 1])
  R <- rbinom(n, 1, p)
  df <- data.frame(Y_miss = Y, X1 = X1, X2 = X2, Z1 = Z1)
  df$Y_miss[R == 0] <- NA_real_
  if (!anyNA(df$Y_miss)) df$Y_miss[1] <- NA_real_
# EL with grouped aux and transformed response predictor
  eng_el <- el_engine(auxiliary_means = c(X1 = 0, X2 = 0), variance_method = "none")
  expect_error(
    nmar(Y_miss ~ (X1 + X2) | (Z1 + I(X1^2)), data = df, engine = eng_el),
    NA
  )
# ET with simple pattern (no overlap: response-only depends on Z1 and its transform)
  eng_et <- exptilt_engine(y_dens = "normal", variance_method = "none", standardize = FALSE)
  expect_error(
    nmar(Y_miss ~ (X1 + X2) | (Z1 + I(Z1^2)), data = df, engine = eng_et),
    NA
  )
})
