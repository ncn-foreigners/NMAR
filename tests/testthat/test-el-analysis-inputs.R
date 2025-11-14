test_that("el_build_input_spec validates weights and respondent metadata", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )

  expect_error(
    NMAR:::el_build_input_spec(
      formula = Y_miss ~ X,
      data = df,
      weights_full = rep(1, 5),
      population_total = NULL,
      is_survey = FALSE,
      design_object = NULL,
      auxiliary_means = NULL
    ),
    "`weights_full` must align",
    fixed = FALSE
  )

  spec <- NMAR:::el_build_input_spec(
    formula = Y_miss ~ X,
    data = df,
    weights_full = NULL,
    population_total = NULL,
    is_survey = FALSE,
    design_object = NULL,
    auxiliary_means = NULL
  )
  dat2 <- if (inherits(spec$analysis_object, "survey.design")) spec$analysis_object$variables else spec$analysis_object
  expect_equal(spec$N_pop, nrow(dat2))
  expect_equal(length(spec$respondent_weights), sum(spec$respondent_mask))
  expect_equal(spec$outcome_var, "Y_miss")
})

test_that("el_validate_design_spec detects tampered respondent mask", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  design <- NMAR:::el_prepare_design(Y_miss ~ X, df, require_na = FALSE)
  tampered <- design
  tampered$respondent_mask <- tampered$respondent_mask[-1]
  expect_error(
    NMAR:::el_validate_design_spec(tampered, data_nrow = nrow(df), context_label = "data frame"),
    "respondent mask length",
    fixed = FALSE
  )
})

test_that("el_build_input_spec captures iid metadata", {
  set.seed(123)
  df <- data.frame(
    Y_miss = c(rnorm(5), NA),
    X = rnorm(6)
  )
  spec <- NMAR:::el_build_input_spec(
    formula = Y_miss ~ X,
    data = df,
    weights_full = NULL,
    population_total = NULL,
    is_survey = FALSE,
    design_object = NULL,
    auxiliary_means = c(X = 0)
  )
  expect_equal(spec$outcome_var, "Y_miss")
  expect_equal(spec$N_pop, nrow(df))
  expect_equal(length(spec$respondent_weights), length(spec$respondent_indices))
  expect_false(spec$is_survey)
  expect_true("..nmar_delta.." %in% names(spec$analysis_object))
  expect_identical(spec$respondent_mask, spec$analysis_object[["..nmar_delta.."]] == 1)
})
