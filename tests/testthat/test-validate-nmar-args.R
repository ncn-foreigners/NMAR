test_that("default traits block outcome on RHS and overlap", {
  df <- data.frame(
    Y = c(1, NA, 3, 4),
    X = rnorm(4),
    Z = rnorm(4)
  )
  spec_outcome <- NMAR:::parse_nmar_spec(Y ~ 1 | Y, df)
  expect_error(
    NMAR:::validate_nmar_args(spec_outcome, NMAR:::NMAR_DEFAULT_TRAITS),
    "Outcome variable cannot appear",
    fixed = FALSE
  )

  spec_overlap <- NMAR:::parse_nmar_spec(Y ~ X | X + Z, df)
  expect_error(
    NMAR:::validate_nmar_args(spec_overlap, NMAR:::NMAR_DEFAULT_TRAITS),
    "Covariate sets must be mutually exclusive",
    fixed = FALSE
  )

  traits_el <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  expect_silent(NMAR:::validate_nmar_args(spec_outcome, traits_el))
  expect_silent(NMAR:::validate_nmar_args(spec_overlap, traits_el))
})

test_that("respondents-only allowance is trait gated", {
  df <- data.frame(Y = rnorm(5), X = rnorm(5))
  spec <- NMAR:::parse_nmar_spec(Y ~ X, df)
  expect_error(
    NMAR:::validate_nmar_args(spec, NMAR:::NMAR_DEFAULT_TRAITS),
    "must contain NA",
    fixed = FALSE
  )
  traits_el <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0), n_total = 10))
  expect_silent(NMAR:::validate_nmar_args(spec, traits_el))
})
