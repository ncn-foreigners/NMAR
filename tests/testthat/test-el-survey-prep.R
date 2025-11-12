skip_if_not_installed("survey")

test_that("survey designs reuse EL prep workflow", {
  set.seed(100)
  df <- data.frame(
    Y_miss = c(1, 2, NA, 4),
    X = rnorm(4),
    w = c(1, 2, 1, 1)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)

  prep <- NMAR:::el_construct_design(Y_miss ~ X, df, require_na = FALSE)
  design_svy <- NMAR:::el_construct_design(Y_miss ~ X, design$variables, require_na = FALSE)

  expect_equal(prep$missingness_model_matrix, design_svy$missingness_model_matrix)
  expect_equal(prep$aux_mm_full[prep$respondent_mask, , drop = FALSE],
               design_svy$aux_mm_full[design_svy$respondent_mask, , drop = FALSE])
  expect_equal(prep$respondent_mask, design_svy$respondent_mask)

  expect_error(NMAR:::el_validate_respondents_only(Y_miss ~ X, df, NULL, context_label = "data frame"), NA)
  expect_error(NMAR:::el_validate_respondents_only(Y_miss ~ X, design$variables, NULL, context_label = "survey design"), NA)
})
