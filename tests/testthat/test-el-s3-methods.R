test_that("tidy/glance produce expected shapes", {
  set.seed(42)
  df <- make_iid_nmar(n = 200, alpha = 0.4, seed = 42)
  fit <- nmar(formula = Y_miss ~ X, data = df,
              engine = make_engine(auxiliary_means = c(X = 0), variance_method = "none"))

  td <- tidy(fit)
  gl <- glance(fit)
  expect_true(is.data.frame(td) && nrow(td) >= 1)
  expect_true(all(c("estimate", "std.error") %in% names(td)))
  expect_true(is.data.frame(gl) && nrow(gl) == 1)
  expect_true(all(c("y_hat", "std.error", "converged") %in% names(gl)))
})
