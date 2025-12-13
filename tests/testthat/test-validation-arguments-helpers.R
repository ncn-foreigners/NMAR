test_that("validator helpers fail with clear errors", {
  expect_error(NMAR:::validator_assert_list(1, name = "x"), "must be a list", fixed = FALSE)
  expect_error(NMAR:::validator_assert_choice(5, choices = 0:3, name = "trace_level"), "should be one of", fixed = FALSE)

  expect_error(NMAR:::validator_assert_positive_number(NA_real_, name = "x"), "single numeric", fixed = FALSE)
  expect_error(NMAR:::validator_assert_positive_number(-1, name = "x"), "positive", fixed = FALSE)
  expect_error(NMAR:::validator_assert_positive_number(Inf, name = "x", allow_infinite = FALSE), "finite", fixed = FALSE)
  expect_silent(NMAR:::validator_assert_positive_number(Inf, name = "x", allow_infinite = TRUE))

  expect_silent(NMAR:::validator_assert_named_numeric(NULL, name = "x", allow_null = TRUE))
  expect_error(NMAR:::validator_assert_named_numeric(1:2, name = "x", allow_null = FALSE), "named numeric", fixed = FALSE)
})
