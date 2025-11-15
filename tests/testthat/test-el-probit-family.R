test_that("EL engine runs with probit family (data.frame)", {
  set.seed(501)
  df <- make_iid_nmar(n = 300, alpha = 0.5, seed = 501)
  eng <- make_engine(auxiliary_means = c(X = 0), variance_method = "none", standardize = FALSE, family = "probit")
  fit <- nmar(formula = Y_miss ~ X, data = df, engine = eng)
  expect_true(fit$converged)
  expect_true(is.finite(fit[['y_hat']]))
  expect_true(is.na(fit[['se']]) || is.finite(fit[['se']]))
})

test_that("Estimating equations solved for probit family (max residual small)", {
  set.seed(502)
  df <- make_iid_nmar(n = 300, alpha = 0.5, seed = 502)
  fit <- nmar(formula = Y_miss ~ X, data = df,
              engine = make_engine(auxiliary_means = c(X = 0), variance_method = "none", family = "probit"))
  expect_true(fit$converged)
  expect_lt(fit$diagnostics$max_equation_residual, 1e-5)
})
