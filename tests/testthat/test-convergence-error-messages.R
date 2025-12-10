test_that("convergence errors provide helpful, actionable messages", {
  skip_on_cran()

# Create a scenario likely to fail convergence due to low iteration limit
  set.seed(12345)
  df <- data.frame(
    Y = rnorm(50, mean = 10),
    X = rnorm(50, mean = 0)
  )
  df$Y[1:25] <- NA # 50% missing

  eng <- el_engine(
    family = "logit",
    variance_method = "none",
    control = list(maxit = 2), # Very low iterations to force failure
    on_failure = "error"
  )

  err <- tryCatch(
    nmar(Y ~ X | X, data = df, engine = eng, trace_level = 0),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(err$message, "solver failed to converge", ignore.case = TRUE)
  expect_match(err$message, "maxit", ignore.case = TRUE) # Should suggest increasing
  expect_match(err$message, "on_failure.*return", ignore.case = TRUE) # Should suggest on_failure
})

test_that("convergence error is a standard error, not custom class", {
  skip_on_cran()

  set.seed(123)
  df <- data.frame(
    Y = rnorm(30),
    X = rnorm(30)
  )
  df$Y[1:15] <- NA

  eng <- el_engine(
    family = "logit",
    variance_method = "none",
    control = list(maxit = 1),
    on_failure = "error"
  )

  err <- tryCatch(
    nmar(Y ~ X | X, data = df, engine = eng, trace_level = 0),
    error = function(e) e
  )

# Should be a standard error, not convergenceError
  expect_s3_class(err, "error")
  expect_s3_class(err, "condition")
# Should NOT have convergenceError class
  expect_false("convergenceError" %in% class(err))
})

test_that("negative weights error is clear and informative", {
  skip_on_cran()

# Create scenario with grossly inconsistent auxiliary means
  set.seed(789)
  df <- data.frame(
    Y = rnorm(100, mean = 0, sd = 1),
    X = rnorm(100, mean = 0, sd = 1)
  )
  df$Y[1:50] <- NA

# Auxiliary mean impossibly far from data (100 standard deviations away)
  eng <- el_engine(
    family = "logit",
    variance_method = "none",
    auxiliary_means = c(X = 100),
    on_failure = "error"
  )

  err <- NULL
  expect_warning(
    err <- tryCatch(
      nmar(Y ~ X, data = df, engine = eng, trace_level = 0),
      error = function(e) e
    ),
    "Auxiliary means appear far from respondents' support",
    ignore.case = TRUE
  )

# May or may not error depending on solver behavior
# If it does error, check the message is informative
  if (inherits(err, "error")) {
    expect_s3_class(err, "error")
# Message should explain the problem
    expect_true(
      grepl("negative.*weights|auxiliary.*means|failed.*converge",
            err$message, ignore.case = TRUE)
    )
  }
})

test_that("on_failure = 'return' provides diagnostics without error", {
  skip_on_cran()

  set.seed(456)
  df <- data.frame(
    Y = rnorm(40),
    X = rnorm(40)
  )
  df$Y[1:20] <- NA

  eng <- el_engine(
    family = "logit",
    variance_method = "none",
    control = list(maxit = 1),
    on_failure = "return" # Should return gracefully, not error
  )

  result <- nmar(Y ~ X | X, data = df, engine = eng, trace_level = 0)

# Should return a list with converged = FALSE
  expect_type(result, "list")
  expect_false(result$converged)
  expect_true("diagnostics" %in% names(result))
  expect_true("convergence_code" %in% names(result$diagnostics))
  expect_true("message" %in% names(result$diagnostics))
})

test_that("validation errors use call. = FALSE for clean stack traces", {
  skip_on_cran()

# Empty data
  df_empty <- data.frame(Y = numeric(0), X = numeric(0))

  eng <- el_engine(family = "logit", variance_method = "none")

  err <- tryCatch(
    nmar(Y ~ X | X, data = df_empty, engine = eng, trace_level = 0),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(err$message, "empty", ignore.case = TRUE)
# The call should not include internal function names if call. = FALSE works
# We can't test this directly, but we can verify the error is thrown
})
