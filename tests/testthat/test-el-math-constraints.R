test_that("estimating equations solved to tight tolerance (max residual small)", {
  set.seed(3402)
  df <- make_iid_nmar(n = 300, alpha = 0.5, seed = 3402)
  fit <- NMAR:::el.data.frame(df, Y_miss ~ X,
    auxiliary_means = c(X = 0), standardize = TRUE,
    trim_cap = Inf, variance_method = "none"
  )
  expect_type(fit$converged, "logical")
  if (isTRUE(fit$converged)) {
    expect_lt(fit$diagnostics$max_equation_residual, 1e-6)
  }
  expect_true(is.finite(fit$diagnostics$min_denominator))
})
