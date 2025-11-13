test_that("el_prepare_analysis_context validates mask and weights", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  design <- NMAR:::el_prepare_design(Y_miss ~ X, df, require_na = FALSE)
  tampered <- design
  tampered$respondent_mask <- tampered$respondent_mask[-1]

  expect_error(
    NMAR:::el_prepare_analysis_context(df, tampered),
    "mask must have the same length",
    fixed = FALSE
  )

  expect_error(
    NMAR:::el_prepare_analysis_context(df, design, weights_full = rep(1, 5)),
    "weights_full",
    fixed = FALSE
  )

  res <- NMAR:::el_prepare_analysis_context(df, design)
  dat2 <- if (inherits(res$analysis_object, "survey.design")) res$analysis_object$variables else res$analysis_object
  expect_equal(res$N_pop, nrow(dat2))
  expect_equal(length(res$respondent_weights), sum(design$respondent_mask))
  expect_equal(design$outcome_var, "Y_miss")
})
