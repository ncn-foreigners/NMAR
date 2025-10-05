test_that("response_predictors can include non-auxiliary variables", {
  set.seed(4242)
  df <- make_iid_nmar(n = 350, alpha = 0.6, include_z = TRUE, seed = 4242)

# Outcome RHS contains only X (auxiliary). Include Z as response-only predictor.
  res <- NMAR:::el.data.frame(
    data = df,
    formula = Y_miss ~ X,
    response_predictors = c("Z"),
    auxiliary_means = c(X = 0),
    standardize = FALSE,
    trim_cap = Inf,
    on_failure = "return",
    variance_method = "delta"
  )

  expect_s3_class(res, "nmar_result_el")
  expect_true(isTRUE(res$converged))
  expect_true(is.numeric(res[['std_error']]))
  expect_true(is.na(res[['std_error']]) || res[['std_error']] >= 0)
# Coefficient vector should include Z as a response predictor
  expect_true(any(grepl("(^|\\b)Z(\\b|$)", names(res$model$coefficients))))
})
