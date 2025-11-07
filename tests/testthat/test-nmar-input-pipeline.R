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
