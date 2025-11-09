test_that("prepare_nmar_design builds response/aux matrices correctly (explicit RHS only)", {
  df <- data.frame(Y = c(1, NA, 3, 4), X = rnorm(4), Z = rnorm(4))
  spec <- NMAR:::parse_nmar_spec(Y ~ X | Z, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task, include_auxiliary = TRUE, include_response = TRUE)

# Response design matrix should contain only explicit RHS predictors (Z) and no intercept
  mm_resp <- design$design_matrices$response
  expect_true(is.matrix(mm_resp))
  expect_false("(Intercept)" %in% colnames(mm_resp))
  expect_false("Y" %in% colnames(mm_resp))
  expect_equal(colnames(mm_resp), "Z")

# Auxiliary matrix should contain only X (no intercept)
  mm_aux <- design$design_matrices$auxiliary
  expect_true(is.matrix(mm_aux))
  expect_false("(Intercept)" %in% colnames(mm_aux))
  expect_equal(colnames(mm_aux), "X")
})

test_that("prepare_nmar_design with no bar yields NULL response matrix and aux-only RHS", {
  df <- data.frame(Y = c(1, NA, 3, 4), X = rnorm(4))
  spec <- NMAR:::parse_nmar_spec(Y ~ X, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task, include_auxiliary = TRUE, include_response = TRUE)
  expect_null(design$design_matrices$response)
  expect_true(is.matrix(design$design_matrices$auxiliary))
  expect_equal(colnames(design$design_matrices$auxiliary), "X")
})

test_that("Non-finite values from LHS transform on observed rows are rejected with clear error", {
  df <- data.frame(Y = c(0, 1, NA, 2), X = rnorm(4))
# log(Y) produces -Inf at row 1 (observed), which should error
  spec <- NMAR:::parse_nmar_spec(log(Y) ~ X, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  expect_error(
    NMAR:::prepare_nmar_design(task),
    "produced non-finite values",
    fixed = FALSE
  )
})
