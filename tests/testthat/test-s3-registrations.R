skip_if_not_installed("generics")

make_test_result <- function() {
  new_nmar_result(
    estimate = 1.2,
    se = 0.05,
    converged = TRUE,
    inference = list(variance_method = "delta"),
    sample = list(n_total = 10L, n_respondents = 8L, is_survey = FALSE, design = NULL),
    diagnostics = list(trimmed_fraction = 0),
    class = "nmar_result_el"
  )
}

test_that("nmar_result S3 generics are registered", {
  res <- make_test_result()

  tidy_df <- generics::tidy(res)
  expect_s3_class(tidy_df, "data.frame")
  expect_true("estimand" %in% tidy_df$component)

  glance_df <- generics::glance(res)
  expect_s3_class(glance_df, "data.frame")
  expect_true(all(c("y_hat", "std.error") %in% names(glance_df)))

  sum_obj <- summary(res)
  expect_s3_class(sum_obj, "summary_nmar_result")

  capture.output(print(res))
  capture.output(print(sum_obj))

# expect_true(utils::isS3stdGeneric("autoplot"))
})
