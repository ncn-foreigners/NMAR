test_that("EL engine runs with probit family (data.frame)", {
  set.seed(501)
  df <- make_iid_nmar(n = 300, alpha = 0.5, seed = 501)
  eng <- make_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = FALSE, family = "probit")
  fit <- nmar(formula = Y_miss ~ X, data = df, engine = eng)
  expect_true(fit$converged)
  expect_true(is.finite(fit[['estimate']]))
  expect_true(is.finite(fit[['std_error']]))
})

test_that("Estimating equations solved for probit family (max residual small)", {
  set.seed(502)
  df <- make_iid_nmar(n = 300, alpha = 0.5, seed = 502)
  fit <- nmar(formula = Y_miss ~ X, data = df,
              engine = make_engine(auxiliary_means = c(X = 0), variance_method = "delta", family = "probit"))
  expect_true(fit$converged)
  expect_lt(fit$diagnostics$max_equation_residual, 1e-5)
})
