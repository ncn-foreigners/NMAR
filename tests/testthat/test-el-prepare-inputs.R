test_that("el_prepare_inputs creates unique delta var name when colliding", {
  df <- data.frame(`..nmar_delta..` = 1:5, Y = c(1, 2, NA, 4, 5), X = rnorm(5))
  res <- NMAR:::el_prepare_inputs(Y ~ X, df, NULL)
  fml <- res$formula_list$response
  lhs <- all.vars(fml)[1]
  expect_true(lhs != "..nmar_delta..")
  expect_true(lhs %in% names(res$data))
})

test_that("el_prepare_inputs expands dot notation via Formula", {
  df <- data.frame(
    Y_miss = c(1, NA, 0, 1),
    X1 = rnorm(4),
    X2 = rnorm(4),
    Z = rnorm(4)
  )
  res <- NMAR:::el_prepare_inputs(Y_miss ~ . | ., df, require_na = FALSE)
  aux_terms <- res$formula_list$auxiliary
  expect_setequal(attr(stats::terms(aux_terms), "term.labels"), c("X1", "X2", "Z"))
  mm <- model.matrix(update(res$formula_list$response, NULL ~ .), data = res$data, na.action = stats::na.pass)
  expect_true(all(c("Y_miss", "X1", "X2", "Z") %in% colnames(mm)))
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
  expect_null(res$formula_list$auxiliary)
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
  aux_terms <- attr(stats::terms(res$formula_list$auxiliary), "term.labels")
  expect_equal(aux_terms, "X")
})
