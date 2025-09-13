test_that("EL engine runs with probit family (data.frame)", {
  set.seed(501)
  N <- 300
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- pnorm(-0.5 + 0.6 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_

  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = FALSE, family = "probit")
  fml <- list(outcome = ~Y_miss, covariates_outcome = ~X, covariates_missingness = ~NULL)
  fit <- nmar(formula = fml, data = df, engine = eng)
  expect_true(fit$converged)
  expect_true(is.finite(fit$y_hat))
  expect_true(is.finite(fit$se))
})

test_that("Estimating equations solved for probit family (max residual small)", {
  set.seed(502)
  N <- 300
  X <- rnorm(N)
  Y <- 2 + 0.3 * X + rnorm(N)
  p <- pnorm(-0.4 + 0.5 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_

  fit <- nmar(
    formula = list(outcome = ~Y_miss, covariates_outcome = ~X, covariates_missingness = ~NULL),
    data = df,
    engine = el_engine(auxiliary_means = c(X = 0), variance_method = "delta", family = "probit")
  )
  expect_true(fit$converged)
  expect_lt(fit$diagnostics$max_equation_residual, 1e-5)
})
