skip_if_not_installed("survey")

test_that("survey designs reuse EL prep workflow", {
  set.seed(100)
  df <- data.frame(
    Y_miss = c(1, 2, NA, 4),
    X = rnorm(4),
    w = c(1, 2, 1, 1)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)

  prep <- NMAR:::el_prepare_design(Y_miss ~ X, df, require_na = FALSE)
  design_svy <- NMAR:::el_prepare_design(Y_miss ~ X, design$variables, require_na = FALSE)

  expect_equal(prep$response, design_svy$response)
  expect_equal(prep$aux_resp, design_svy$aux_resp)
  expect_equal(prep$mask, design_svy$mask)

  expect_error(
    NMAR:::el_validate_respondents_only(Y_miss ~ X, df, NULL, context_label = "data frame"),
    NA
  )
  expect_error(
    NMAR:::el_validate_respondents_only(Y_miss ~ X, design$variables, NULL, context_label = "survey design"),
    NA
  )
})
