test_that("Response-model vcov names align with coefficients", {
  set.seed(7003)
  df <- make_iid_nmar(n = 300, alpha = 0.5, seed = 7003)
  fit <- nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = make_engine(auxiliary_means = c(X = 0), variance_method = "none", standardize = TRUE)
  )
  beta <- fit$model$coefficients
  V <- fit$model$vcov
  if (!is.null(beta) && !is.null(V)) {
    expect_true(all(colnames(V) == names(beta)))
    expect_true(all(rownames(V) == names(beta)))
  }
})
