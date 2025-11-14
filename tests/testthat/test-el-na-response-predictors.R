test_that("EL errors early when response-model predictors have NA among respondents (with transforms)", {
  set.seed(42)
  n <- 60
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 0.5 + 0.3 * X + 0.2 * Z + rnorm(n)
  p <- plogis(-0.2 + 0.5 * scale(Y)[, 1] + 0.1 * X)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_
# Inject NA in a response predictor among respondents
  idx <- which(!is.na(df$Y_miss))[1]
  df$Z[idx] <- NA_real_

  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
  expect_error(
    nmar(Y_miss ~ X | I(Z^2), data = df, engine = eng),
    regexp = "Missingness-model predictor",
    info = "NA in transformed response predictor should be trapped early"
  )
})

test_that("EL allows NA in response-model predictors only among nonrespondents", {
  set.seed(43)
  n <- 60
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 0.5 + 0.3 * X + 0.2 * Z + rnorm(n)
  p <- plogis(-0.2 + 0.5 * scale(Y)[, 1] + 0.1 * X)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_
# Inject NA in a response predictor at a nonrespondent row
  idx_nonresp <- which(is.na(df$Y_miss))[1]
  df$Z[idx_nonresp] <- NA_real_

  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
  fit <- nmar(Y_miss ~ X | Z, data = df, engine = eng)
  expect_s3_class(fit, "nmar_result_el")
  expect_true(isTRUE(fit$converged))
})
