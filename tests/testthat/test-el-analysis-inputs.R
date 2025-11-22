test_that("el_prepare_inputs validates weights and respondent metadata", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )

  expect_error(
    el_prepare_inputs(
      formula = Y_miss ~ X,
      data = df,
      weights = rep(1, 5)
    ),
    "`weights` must align",
    fixed = FALSE
  )

  spec <- el_prepare_inputs(
    formula = Y_miss ~ X,
    data = df,
    weights = NULL
  )
  expect_equal(spec$N_pop, nrow(spec$analysis_data))
  expect_equal(length(spec$respondent_weights), sum(spec$respondent_mask))
  expect_equal(spec$outcome_expr, "Y_miss")
})

test_that("el_validate_design_spec detects tampered respondent mask", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  design <- el_prepare_inputs(Y_miss ~ X, df)
  tampered <- design
  tampered$respondent_mask <- tampered$respondent_mask[-1]
  expect_error(
    el_validate_design_spec(tampered, data_nrow = nrow(df)),
    "respondent mask length",
    fixed = FALSE
  )
})

test_that("el_prepare_inputs captures iid metadata", {
  set.seed(123)
  df <- data.frame(
    Y_miss = c(rnorm(5), NA),
    X = rnorm(6)
  )
  spec <- el_prepare_inputs(
    formula = Y_miss ~ X,
    data = df,
    weights = NULL
  )
  expect_equal(spec$outcome_expr, "Y_miss")
  expect_equal(spec$N_pop, nrow(df))
  expect_equal(length(spec$respondent_weights), sum(spec$respondent_mask))
  expect_true("..nmar_delta.." %in% names(spec$analysis_data))
  expect_identical(spec$respondent_mask, spec$analysis_data[["..nmar_delta.."]] == 1)
})
