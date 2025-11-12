test_that("el_prepare_inputs creates unique delta var name when colliding", {
  df <- data.frame(`..nmar_delta..` = 1:5, Y = c(1, 2, NA, 4, 5), X = rnorm(5))
  res <- NMAR:::el_prepare_inputs(Y ~ X, df, NULL)
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
  res <- NMAR:::el_prepare_inputs(Y_miss ~ . | ., df, require_na = FALSE)
  expect_setequal(colnames(res$auxiliary_matrix_full), c("X1", "X2", "Z"))
  expect_true(all(c("(Intercept)", "Y_miss", "X1", "X2", "Z") %in% colnames(res$response_matrix)))
})

test_that("dot expansion drops outcome-derived auxiliary terms", {
  set.seed(99)
  df <- data.frame(
    Y_miss = c(rnorm(3), NA),
    X = rnorm(4),
    F = factor(c("a", "b", "a", "b"))
  )
  res <- NMAR:::el_prepare_inputs(Y_miss ~ . | ., df, require_na = FALSE)
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
    NMAR:::el_prepare_inputs(Y_miss ~ I(Y_miss^2) + X | Z, df),
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
    res <- NMAR:::el_prepare_inputs(Y_miss ~ 1 | Z, df),
    NA
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
    res <- NMAR:::el_prepare_inputs(Y_miss ~ 1 + X, df),
    "dropping the requested intercept",
    fixed = FALSE
  )
  expect_equal(colnames(res$auxiliary_matrix), "X")
})
