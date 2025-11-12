test_that("el_prepare_analysis_inputs validates mask and weights", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  mask <- c(TRUE, FALSE, TRUE, FALSE)

  expect_error(
    NMAR:::el_prepare_analysis_inputs(df, "Y_miss", mask[-1]),
    "mask must have the same length",
    fixed = FALSE
  )

  expect_error(
    NMAR:::el_prepare_analysis_inputs(df, "Y_miss", mask, weights_full = rep(1, 5)),
    "weights_full",
    fixed = FALSE
  )

  res <- NMAR:::el_prepare_analysis_inputs(df, "Y_miss", mask, variance_method = "none")
  expect_equal(res$N_pop, nrow(res$data_aug))
  expect_equal(length(res$respondent_weights), sum(mask))
  expect_equal(res$data_info$outcome_var, "Y_miss")
})
