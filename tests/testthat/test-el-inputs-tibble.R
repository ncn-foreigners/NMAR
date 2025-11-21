test_that("EL input pipeline accepts tibble data", {
  skip_if_not_installed("tibble")
  df <- make_el_test_data(n = 120)
  tbl <- tibble::as_tibble(df)

  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")

  fit <- nmar(Y_miss ~ X | Z, data = tbl, engine = eng)
  expect_s3_class(fit, "nmar_result_el")
  expect_type(fit$converged, "logical")
})
