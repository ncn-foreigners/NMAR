skip_if_not_installed("survey")

test_that("survey designs reuse EL prep workflow", {
  set.seed(100)
  df <- data.frame(
    Y_miss = c(1, 2, NA, 4),
    X = rnorm(4),
    w = c(1, 2, 1, 1)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)

  prep <- NMAR:::el_parse_design(Y_miss ~ X, df, require_na = FALSE)
  design_svy <- NMAR:::el_parse_design(Y_miss ~ X, design$variables, require_na = FALSE)

  expect_equal(prep$response_matrix, design_svy$response_matrix)
  expect_equal(prep$aux_full[prep$respondent_mask, , drop = FALSE],
               design_svy$aux_full[design_svy$respondent_mask, , drop = FALSE])
  expect_equal(prep$respondent_mask, design_svy$respondent_mask)

  expect_error(NMAR:::el_validate_respondents_only(Y_miss ~ X, df, NULL, context_label = "data frame"), NA)
  expect_error(NMAR:::el_validate_respondents_only(Y_miss ~ X, design$variables, NULL, context_label = "survey design"), NA)
})

test_that("survey prep stores delta column and uses rescaled weights", {
  skip_if_not_installed("survey")
  set.seed(200)
  df <- data.frame(
    Y_miss = c(rnorm(5), NA),
    X = rnorm(6),
    w = c(1, 2, 3, 4, 5, 6)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  n_total <- 2 * sum(weights(des))
  res <- NMAR:::el.survey.design(
    data = des,
    formula = Y_miss ~ X,
    auxiliary_means = c(X = 0),
    n_total = n_total,
    variance_method = "none"
  )
  des_after <- res$sample$design
  expect_true("..nmar_delta.." %in% names(des_after$variables))
  expect_equal(sum(stats::weights(res, scale = "population")), n_total)
})
