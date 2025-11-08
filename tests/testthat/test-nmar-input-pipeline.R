test_that("parse_nmar_spec captures core pieces", {
  df <- data.frame(Y = c(1, NA, 3), X = 1:3, Z = 4:6)
  spec <- NMAR:::parse_nmar_spec(Y ~ X | Z + Z, df)
  expect_s3_class(spec, "nmar_input_spec")
  expect_equal(spec$outcome, "Y")
  expect_equal(spec$auxiliary_vars, "X")
  expect_equal(spec$response_predictors, "Z")
  expect_false(spec$is_survey)
  expect_true(is.data.frame(spec$data))
})

test_that("parse_nmar_spec works for survey designs", {
  skip_if_not_installed("survey")
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3), w = c(1, 1.5, 0.5))
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  spec <- NMAR:::parse_nmar_spec(Y ~ X, design)
  expect_true(spec$is_survey)
  expect_s3_class(spec$original_data, "survey.design")
  expect_equal(nrow(spec$data), nrow(df))
})

test_that("validate_nmar_args enforces defaults", {
  df <- data.frame(Y = c(1, NA, 3, 4), X = rnorm(4), Z = rnorm(4))
  spec <- NMAR:::parse_nmar_spec(Y ~ X | X, df)
  expect_error(NMAR:::validate_nmar_args(spec, list()), "mutually exclusive", fixed = FALSE)
})

test_that("validate_nmar_args can be relaxed for EL", {
  df <- data.frame(Y = c(1, NA, 3, 4), X = rnorm(4), Z = rnorm(4))
  spec <- NMAR:::parse_nmar_spec(Y ~ X | Y, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  expect_equal(spec$response_predictors_raw, "Y")
})

test_that("default traits block outcome on response model", {
  df <- data.frame(Y = c(1, NA, 2), X = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(Y ~ 1 | Y, df)
  expect_error(
    NMAR:::validate_nmar_args(spec, list()),
    "Outcome variable cannot appear",
    fixed = FALSE
  )
})

test_that("nonparam engine accepts multi-outcome formulas", {
  df <- data.frame(
    Voted_A = c(10, 12),
    Voted_B = c(5, 8),
    Other = c(2, 1),
    Gender = c(0, 1),
    Refusal = c(3, 4)
  )
  spec <- NMAR:::parse_nmar_spec(Voted_A + Voted_B + Other ~ Gender, df)
  traits_np <- NMAR:::engine_traits(exptilt_nonparam_engine(refusal_col = "Refusal"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits_np))
  traits_default <- NMAR:::engine_traits(exptilt_engine())
  expect_error(NMAR:::validate_nmar_args(spec, traits_default))
})

test_that("parse_nmar_spec resolves dots and factors", {
  df <- data.frame(
    Y = c(1, NA, 3),
    X = 1:3,
    Z = rnorm(3),
    F = factor(c("a", "b", "a"))
  )
  spec <- NMAR:::parse_nmar_spec(Y ~ . | F, df)
  expect_setequal(spec$auxiliary_vars, c("X", "Z", "F"))
  expect_equal(spec$response_predictors, "F")
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
})

test_that("partition rebuild preserves formula environments", {
  local_fn <- function(x) x^2
  env <- environment()
  df <- data.frame(
    Y = c(1, 2, NA, 4),
    X = rnorm(4)
  )
  formula <- Y ~ local_fn(X)
  environment(formula) <- env
  spec <- NMAR:::parse_nmar_spec(formula, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)
  rebuilt <- design$user_formula
  expect_identical(attr(rebuilt, ".Environment"), env)
  mm <- stats::model.matrix(rebuilt, df)
  expect_true(any(grepl("local_fn", colnames(mm))))
})

test_that("new_nmar_task preserves raw predictor metadata", {
  df <- data.frame(Y = c(1, NA, 3), X = 1:3, Z = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(Y ~ X | Z + Y, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  task <- NMAR:::new_nmar_task(spec, traits)
  expect_equal(task$response_predictors_raw, c("Z", "Y"))
  expect_equal(task$auxiliary_vars_raw, "X")
})

test_that("outcome transformations must be numeric vectors", {
  df <- data.frame(Y = c(1, 2, NA, 4), X = rnorm(4))
  spec <- NMAR:::parse_nmar_spec(I(as.character(Y)) ~ X, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  expect_error(
    NMAR:::prepare_nmar_design(task),
    "numeric vector",
    fixed = FALSE
  )
})

test_that("survey designs pick up derived outcome columns", {
  skip_if_not_installed("survey")
  df <- data.frame(
    Y = c(1, 2, NA, 4),
    X = rnorm(4),
    w = c(1, 1.5, 0.5, 2)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  spec <- NMAR:::parse_nmar_spec(log(Y + 10) ~ X, design)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  info <- NMAR:::prepare_nmar_design(task)
  expect_true(info$is_survey)
  expect_true(info$outcome_column %in% names(info$survey_design$variables))
  expect_equal(info$survey_design$variables[[info$outcome_column]], info$data[[info$outcome_column]])
})

test_that("multi-outcome validation enforces trait restrictions", {
  df <- data.frame(
    Y1 = c(1, 2, 3),
    Y2 = c(4, 5, 6),
    X = c(0.1, 0.2, 0.3),
    Z = c(1, 1, 0)
  )
  spec <- NMAR:::parse_nmar_spec(Y1 + Y2 ~ X | Z, df)
  traits_np <- NMAR:::engine_traits(exptilt_nonparam_engine(refusal_col = "Z"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits_np))
  traits_block <- NMAR:::engine_traits(exptilt_engine())
  expect_error(NMAR:::validate_nmar_args(spec, traits_block))
})

test_that("respondents-only inputs require allow_respondents_only trait", {
  df <- data.frame(Y = rnorm(5), X = rnorm(5))
  spec <- NMAR:::parse_nmar_spec(Y ~ X, df)
  traits_default <- NMAR:::engine_traits(exptilt_engine())
  expect_error(NMAR:::validate_nmar_args(spec, traits_default), "must contain NA", fixed = FALSE)
  traits_el <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0), n_total = 100))
  expect_silent(NMAR:::validate_nmar_args(spec, traits_el))
})
