test_that("EL supports response-only with no auxiliaries (y ~ 1 | z)", {
  set.seed(1001)
  n <- 120
  Z <- rnorm(n)
  Y <- 1 + 0.3 * Z + rnorm(n, sd = 0.4)
  respond <- rbinom(n, 1, plogis(-0.2 + 0.5 * scale(Y)[, 1] + 0.4 * Z))
  df <- data.frame(Y_miss = Y, Z = Z)
  df$Y_miss[respond == 0] <- NA_real_
  eng <- el_engine(auxiliary_means = NULL, variance_method = "none", standardize = FALSE)
  fit <- nmar(Y_miss ~ 1 | Z, data = df, engine = eng)
  expect_s3_class(fit, "nmar_result_el")
  expect_true(isTRUE(fit$converged))
  obs <- which(respond == 1)
  w <- weights(fit)
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_equal(sum(w * df$Y_miss[obs]), fit$y_hat, tolerance = 1e-8)
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

test_that("Engines accept parentheses and transforms in RHS partitions", {
  set.seed(1004)
  n <- 80
  X1 <- rnorm(n); X2 <- rnorm(n); Z1 <- rnorm(n)
  Y <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
  p <- plogis(-0.5 + 0.5 * scale(Y)[, 1])
  R <- rbinom(n, 1, p)
  df <- data.frame(Y_miss = Y, X1 = X1, X2 = X2, Z1 = Z1)
  df$Y_miss[R == 0] <- NA_real_
# EL with grouped aux and transformed response predictor
  eng_el <- el_engine(auxiliary_means = c(X1 = 0, X2 = 0), variance_method = "none")
  fit_el <- nmar(Y_miss ~ (X1 + X2) | (Z1 + I(X1^2)), data = df, engine = eng_el)
  expect_true(isTRUE(fit_el$converged))
  resp_terms <- fit_el$model$coefficients
  expect_true(any(grepl("I\\(X1\\^2\\)", names(resp_terms))))
# ET with simple pattern (no overlap: response-only depends on Z1 and its transform)
  eng_et <- exptilt_engine(y_dens = "normal", variance_method = "delta", standardize = FALSE)
  fit_et <- nmar(Y_miss ~ (X1 + X2) | (Z1 + I(Z1^2)), data = df, engine = eng_et)
  expect_true(isTRUE(fit_et$converged))
  rn <- names(coef(fit_et))
  expect_true(any(grepl("Z1$", rn)))
  expect_true(any(grepl("I\\(Z1\\^2\\)", rn)))
  expect_false(any(grepl("X1|X2", rn)))
})
