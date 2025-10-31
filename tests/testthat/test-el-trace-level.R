test_that("EL trace_level parameter works at all levels", {
  skip_on_cran()

# Create simple test data
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y_true <- 10 + 2 * x1 + 3 * x2 + rnorm(n)

# Create NMAR missingness
  eta_r <- -1 + 0.5 * x1 - 0.3 * y_true
  prob_response <- plogis(eta_r)
  R <- rbinom(n, 1, prob_response)
  y_obs <- ifelse(R == 1, y_true, NA)

  df <- data.frame(y = y_obs, x1 = x1, x2 = x2)

# Test trace_level = 0 (silent)
  engine <- el_engine(family = 'logit', variance_method = 'none')
  expect_silent({
    result0 <- nmar(y ~ x1 + x2 | x1, data = df, engine = engine, trace_level = 0)
  })
  expect_true(inherits(result0, "nmar_result"))
  expect_true(result0$converged)

# Test trace_level = 1 (basic output)
# We expect output but we don't test exact content to avoid brittle tests
  expect_output({
    result1 <- nmar(y ~ x1 + x2 | x1, data = df, engine = engine, trace_level = 1)
  }, "EMPIRICAL LIKELIHOOD ESTIMATION")
  expect_true(inherits(result1, "nmar_result"))
  expect_true(result1$converged)

# Test trace_level = 2 (detailed output)
  expect_output({
    result2 <- nmar(y ~ x1 + x2 | x1, data = df, engine = engine, trace_level = 2)
  }, "MODEL SPECIFICATION")
  expect_true(inherits(result2, "nmar_result"))
  expect_true(result2$converged)

# Test trace_level = 3 (full detail)
  expect_output({
    result3 <- nmar(y ~ x1 + x2 | x1, data = df, engine = engine, trace_level = 3)
  }, "DETAILED DIAGNOSTICS")
  expect_true(inherits(result3, "nmar_result"))
  expect_true(result3$converged)

# All results should be numerically identical
  expect_equal(result0$estimate, result1$estimate)
  expect_equal(result1$estimate, result2$estimate)
  expect_equal(result2$estimate, result3$estimate)
})

test_that("EL trace_level shows solver restart messages at level 3", {
  skip_on_cran()
  skip_if_not_installed("survey")

# This test uses a scenario that sometimes requires restarts
  set.seed(456)
  n <- 50
  x1 <- rnorm(n)
  y_true <- 5 + x1 + rnorm(n, sd = 0.5)

# Create difficult NMAR scenario
  eta_r <- -2 + 0.8 * x1 - 0.5 * y_true
  prob_response <- plogis(eta_r)
  R <- rbinom(n, 1, prob_response)
  y_obs <- ifelse(R == 1, y_true, NA)

  df <- data.frame(y = y_obs, x1 = x1)

# At level 3, if restarts occur, we should see restart messages
# We don't require restarts to happen, but if they do, level 3 should show them
  engine <- el_engine(family = 'logit', variance_method = 'none', on_failure = 'return')
  output <- capture.output({
    result <- nmar(y ~ x1 | x1, data = df, engine = engine, trace_level = 3)
  })

# Check that we got some output
  expect_true(length(output) > 0)

# If solver worked, result should exist
  if (!is.null(result)) {
    expect_true(inherits(result, "nmar_result"))
  }
})

test_that("EL trace_level parameter validates correctly", {
  skip_on_cran()

  set.seed(789)
  df <- data.frame(
    y = c(1, 2, NA, NA, 5),
    x = c(1, 2, 3, 4, 5)
  )

  engine <- el_engine(family = 'logit', variance_method = 'none')

# Valid trace_level values should work
  expect_no_error(nmar(y ~ x | x, data = df, engine = engine, trace_level = 0))
  expect_no_error(nmar(y ~ x | x, data = df, engine = engine, trace_level = 1))
  expect_no_error(nmar(y ~ x | x, data = df, engine = engine, trace_level = 2))
  expect_no_error(nmar(y ~ x | x, data = df, engine = engine, trace_level = 3))

# Invalid trace_level values should error
  expect_error(nmar(y ~ x | x, data = df, engine = engine, trace_level = -1))
  expect_error(nmar(y ~ x | x, data = df, engine = engine, trace_level = 4))
  expect_error(nmar(y ~ x | x, data = df, engine = engine, trace_level = "high"))
})
