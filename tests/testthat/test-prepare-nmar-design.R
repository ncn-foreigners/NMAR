test_that("prepare_nmar_design materializes matrices with transforms", {
  set.seed(123)
  df <- data.frame(
    Y = c(1, 2, NA, 4, NA),
    X = rnorm(5),
    Z = factor(c("a", "b", "a", "b", "a"))
  )
  formula <- log(Y + 5) ~ I(X^2) | Z
  spec <- NMAR:::parse_nmar_spec(formula, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)

  expect_equal(design$outcome_label, "log(Y + 5)")
  expect_true(design$outcome_column %in% names(design$data))
  expect_true(all(is.na(df$Y) == is.na(design$data[[design$outcome_column]])))
  expect_identical(design$outcome_column, design$outcome_columns[[1]])
  expect_equal(design$user_formula, formula)
  expect_identical(
    design$engine_formula[[2]],
    as.name(design$outcome_column)
  )

  aux_matrix <- design$design_matrices$auxiliary
  response_matrix <- design$design_matrices$response
  expect_setequal(colnames(response_matrix), c("Za", "Zb"))
  expect_true(any(grepl("I\\(X\\^2\\)", colnames(aux_matrix))))
  expect_true(any(grepl("Zb", colnames(response_matrix), fixed = TRUE)))
})

test_that("prepare_nmar_design enforces weight length and reuses survey weights", {
  skip_if_not_installed("survey")
  df <- data.frame(
    Y = c(1, 2, NA, 4),
    X = rnorm(4),
    w = c(1, 2, 1, 2)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  spec <- NMAR:::parse_nmar_spec(Y ~ X, design)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design_info <- NMAR:::prepare_nmar_design(task)
  expect_true(design_info$is_survey)
  expect_equal(design_info$weights, df$w)
  expect_s3_class(design_info$survey_design, "survey.design")

  expect_error(
    NMAR:::prepare_nmar_design(task, design_weights = c(1, 2)),
    "must have the same length",
    fixed = FALSE
  )
})

test_that("nmar_materialize_design_matrices is method agnostic", {
  df <- data.frame(
    Y = rnorm(5),
    X = rnorm(5),
    Z = rnorm(5)
  )
  spec <- NMAR:::parse_nmar_spec(Y ~ X + Z, df)
  bp <- spec$blueprint
  mats <- NMAR:::nmar_materialize_design_matrices(bp, df)
  expect_true(all(colnames(mats$response) %in% c("(Intercept)", "X", "Z")))
  expect_true(all(colnames(mats$auxiliary) %in% c("X", "Z")))
})

test_that("prepare_nmar_design honors include flags", {
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3), Z = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(Y ~ X | Z, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  task <- NMAR:::new_nmar_task(spec, traits)
  design_aux_only <- NMAR:::prepare_nmar_design(task, include_response = FALSE)
  expect_equal(design_aux_only$response_predictors, character())
  expect_null(design_aux_only$design_matrices$response)
  design_resp_only <- NMAR:::prepare_nmar_design(task, include_auxiliary = FALSE)
  expect_equal(design_resp_only$auxiliary_vars, character())
  expect_null(design_resp_only$design_matrices$auxiliary)
  expect_null(design_aux_only$response_rhs_lang)
  expect_null(design_resp_only$aux_rhs_lang)
})

test_that("prepare_nmar_design supports multi-outcome expressions", {
  df <- data.frame(
    Y1 = c(1, 2, NA, 4),
    Y2 = c(NA, 3, 4, 5),
    X = rnorm(4)
  )
  formula <- cbind(Y1, Y2) ~ X
  traits <- utils::modifyList(NMAR:::NMAR_DEFAULT_TRAITS, list(requires_single_outcome = FALSE))
  spec <- NMAR:::parse_nmar_spec(formula, df, traits = traits)
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)
  expect_length(design$outcome_columns, 2)
  expect_true(all(design$outcome_columns %in% names(design$data)))
  expect_match(deparse(design$engine_formula[[2]]), "cbind", fixed = TRUE)
  expect_equal(design$user_formula, formula)
})
