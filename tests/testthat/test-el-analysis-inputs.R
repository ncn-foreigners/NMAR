test_that("el_prepare_analysis_context validates mask and weights", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  design <- NMAR:::el_prepare_design(Y_miss ~ X, df, require_na = FALSE)
  tampered <- design
  tampered$respondent_mask <- tampered$respondent_mask[-1]

  expect_error(
    NMAR:::el_prepare_analysis_context(df, tampered, variance_method = "none"),
    "mask must have the same length",
    fixed = FALSE
  )

  expect_error(
    NMAR:::el_prepare_analysis_context(df, design, weights_full = rep(1, 5), variance_method = "none"),
    "weights_full",
    fixed = FALSE
  )

  res <- NMAR:::el_prepare_analysis_context(df, design, variance_method = "none")
  expect_equal(res$N_pop, nrow(res$data_aug))
  expect_equal(length(res$respondent_weights), sum(design$respondent_mask))
  expect_equal(design$outcome_var, "Y_miss")
})
