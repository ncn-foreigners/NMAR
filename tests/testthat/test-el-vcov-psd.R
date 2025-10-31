test_that("vcov(beta) is symmetric PSD with aligned names (benign aux cell)", {
  skip("EL delta variance disabled; vcov not available")
  set.seed(7301)
  df <- make_iid_nmar(n = 600, alpha = 0.4, seed = 7301)
  fit <- nmar(Y_miss ~ X, data = df, engine = make_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = TRUE))
  expect_true(isTRUE(fit$converged))
  V <- fit$model$vcov
  beta <- fit$model$coefficients
  expect_true(is.matrix(V) && length(beta) == ncol(V))
  expect_true(all(rownames(V) == names(beta)))
  expect_true(all(colnames(V) == names(beta)))
# Symmetry and PSD up to numerical tolerance
  expect_true(isTRUE(all.equal(V, t(V), tolerance = 1e-10)))
  ev <- suppressWarnings(eigen(0.5 * (V + t(V)), symmetric = TRUE, only.values = TRUE)$values)
  expect_gte(min(ev), -1e-10)
})
