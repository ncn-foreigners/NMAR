skip_if_not_installed("survey")

test_that("EL survey rejects negative design weights", {
  df <- data.frame(
    y_miss = c(1, NA, 2),
    w = c(1, -2, 3)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  expect_error(
    el.survey.design(
      data = des,
      formula = y_miss ~ 1,
      variance_method = "none"
    ),
    regexp = "weights must be nonnegative",
    fixed = FALSE
  )
})

test_that("EL survey strata extraction prefers design$strata over getCall", {
  set.seed(111)
  df <- data.frame(
    y_miss = c(1, NA, 2, NA, 3, NA),
    strata = factor(c("A", "A", "B", "B", "B", "B")),
    w = c(2, 2, 3, 3, 3, 3)
  )
  des <- survey::svydesign(ids = ~1, strata = ~strata, weights = ~w, data = df)
  des$call <- NULL

  eng <- el_engine(variance_method = "none", strata_augmentation = TRUE)
  fit <- nmar(y_miss ~ 1, data = des, engine = eng)

  aux_mat <- fit$diagnostics$auxiliary_matrix
  expect_true(any(grepl("^strata_", colnames(aux_mat))))
})
