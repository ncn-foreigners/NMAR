set.seed(2025)

test_that("EL engine runs and returns expected structure (data.frame)", {
  dat <- make_iid_nmar(n = 400, alpha = 0.4, seed = 2025)
  eng <- make_engine(variance_method = "delta", auxiliary_means = c(X = 0))
  res <- nmar(formula = Y_miss ~ X, data = dat, engine = eng)

  expect_s3_class(res, "nmar_result_el")
  expect_true(isTRUE(res$converged))
  expect_true(is.numeric(res[['estimate']]))
  expect_true(is.numeric(res[['std_error']]))

  est <- res[['estimate']]
  expect_equal(as.numeric(est), res[['estimate']])

  ci <- confint(res)
  expect_true(is.matrix(ci) && nrow(ci) == 1)
  expect_equal(rownames(ci), "Y_miss")

  w <- weights(res)
  expect_true(is.numeric(w) || length(w) == 0)
})
