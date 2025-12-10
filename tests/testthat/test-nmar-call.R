test_that("nmar result metadata records the user call", {
  skip_if_not_installed("nleqslv")

  set.seed(1)
  df <- data.frame(
    Y = c(1, NA, 2, NA, 3),
    X = rnorm(5)
  )
  eng <- el_engine(variance_method = "none", standardize = FALSE)

  fit <- nmar(Y ~ X, data = df, engine = eng)
  expect_s3_class(fit, "nmar_result")
  expect_true(is.call(fit$meta$call))
  expect_identical(fit$meta$call[[1]], as.name("nmar"))
  expect_equal(fit$meta$call$formula, Y ~ X)
  expect_true(identical(environment(fit$meta$call$formula), environment(Y ~ X)))
})

test_that("engine call is preserved in result metadata", {
  skip_if_not_installed("nleqslv")

  set.seed(2)
  df <- data.frame(
    Y = c(1, NA, 2, 4, NA),
    X = rnorm(5)
  )
  eng <- el_engine(variance_method = "none", standardize = FALSE)
  fit <- nmar(Y ~ X, data = df, engine = eng)

# meta$engine_call should capture the inner engine call (el(...))
  expect_true(is.call(fit$meta$engine_call))
# Method dispatch may record either 'el' or 'el.data.frame'
  expect_true(identical(fit$meta$engine_call[[1]], as.name("el")) ||
                identical(fit$meta$engine_call[[1]], as.name("el.data.frame")))
# The top-level call remains nmar(...)
  expect_identical(fit$meta$call[[1]], as.name("nmar"))
})

test_that("formula.nmar_result returns the estimation formula", {
  skip_if_not_installed("nleqslv")

  set.seed(3)
  df <- data.frame(
    Y = c(1, NA, 2, 5, NA),
    X = rnorm(5)
  )
  eng <- el_engine(variance_method = "none", standardize = FALSE)
  fit <- nmar(Y ~ X, data = df, engine = eng)
  expect_s3_class(fit, "nmar_result")
  expect_equal(formula(fit), Y ~ X)
})
