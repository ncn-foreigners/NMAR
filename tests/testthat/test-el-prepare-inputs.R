test_that("el_prepare_inputs creates unique delta var name when colliding", {
  df <- data.frame(`..nmar_delta..` = 1:5, Y = c(1, 2, NA, 4, 5), X = rnorm(5))
  res <- prepare_el_inputs(Y ~ X, df, NULL)
  expect_true(res$delta_name != "..nmar_delta..")
  expect_true(res$delta_name %in% names(res$data))
})

test_that("el_prepare_inputs expands dot notation via Formula", {
  df <- data.frame(
    Y_miss = c(1, NA, 0, 1),
    X1 = rnorm(4),
    X2 = rnorm(4),
    Z = rnorm(4)
  )
  res <- prepare_el_inputs(Y_miss ~ . | ., df, require_na = FALSE)
  expect_setequal(colnames(res$auxiliary_matrix_full), c("X1", "X2", "Z"))
  expect_true(all(c("(Intercept)", "Y_miss", "X1", "X2", "Z") %in% colnames(res$response_matrix)))
})

test_that("response intercept is retained even when formula uses +0", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    Z = rnorm(4)
  )
  expect_warning(
    design <- NMAR:::el_construct_design(Y_miss ~ 1 | Z + 0, df, require_na = FALSE),
    "Missingness-model intercept is required",
    fixed = FALSE
  )
  expect_identical(colnames(design$missingness_model_matrix)[1], "(Intercept)")
  expect_true(all(design$missingness_model_matrix[, "(Intercept)"] == 1))
  expect_true("Z" %in% colnames(design$missingness_model_matrix))
  expect_equal(ncol(design$missingness_model_matrix), 3)
})

test_that("respondents-only data frame requires n_total", {
  df <- data.frame(Y_miss = 1:4)
  eng <- el_engine(variance_method = "none")
  expect_error(
    nmar(Y_miss ~ 1, data = df, engine = eng),
    "Respondents-only data detected",
    fixed = FALSE
  )
  eng_with_total <- el_engine(variance_method = "none", n_total = 100)
  expect_warning(
    expect_s3_class(
      nmar(Y_miss ~ 1, data = df, engine = eng_with_total),
      "nmar_result_el"
    ),
    "Auxiliary intercepts are ignored",
    fixed = FALSE
  )
})

test_that("dot expansion drops outcome-derived auxiliary terms", {
  set.seed(99)
  df <- data.frame(
    Y_miss = c(rnorm(3), NA),
    X = rnorm(4),
    F = factor(c("a", "b", "a", "b"))
  )
  res <- prepare_el_inputs(Y_miss ~ . | ., df, require_na = FALSE)
  expect_false("Y_miss" %in% colnames(res$auxiliary_matrix_full))
  expect_false(anyNA(res$auxiliary_matrix_full))
})

test_that("el_prepare_inputs forbids outcome in auxiliary constraints", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4),
    Z = rnorm(4)
  )
  expect_error(
    prepare_el_inputs(Y_miss ~ I(Y_miss^2) + X | Z, df),
    "outcome cannot appear",
    fixed = FALSE
  )
})

test_that("intercept-only auxiliaries are ignored without warnings", {
  df <- data.frame(
    Y_miss = c(1, NA, 3, NA),
    Z = rnorm(4)
  )
  expect_warning(
    res <- prepare_el_inputs(Y_miss ~ 1 | Z, df),
    "intercept",
    fixed = FALSE
  )
  expect_false(res$has_aux)
  expect_equal(ncol(res$auxiliary_matrix), 0)
})

test_that("explicit intercept plus auxiliaries triggers a warning", {
  df <- data.frame(
    Y_miss = c(1, NA, 3, NA),
    X = rnorm(4)
  )
  expect_warning(
    res <- prepare_el_inputs(Y_miss ~ 1 + X, df),
    "Auxiliary intercepts are ignored",
    fixed = FALSE
  )
  expect_equal(colnames(res$auxiliary_matrix), "X")
})
