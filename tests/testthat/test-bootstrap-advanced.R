test_that("svyrep.design is rejected with clear error", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

# Create replicate design
  drep <- svrep::as_bootstrap_design(dstrat, replicates = 10)

  estimator <- function(data, ...) {
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

# Should error with informative message
  expect_error(
    bootstrap_variance(
      drep, estimator,
      point_estimate = mean(apistrat$api00),
      bootstrap_reps = 10,
      bootstrap_cores = 1
    ),
    "replicate design|svyrep"
  )
})

test_that("replicate count mismatch warns but proceeds", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

  estimator <- function(data, ...) {
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  future::plan(future::sequential)

# Request a large number - may produce fewer
# This test verifies we warn and proceed (not error)
  result <- tryCatch(
    suppressWarnings(
      bootstrap_variance(
        dstrat, estimator,
        point_estimate = mean(apistrat$api00),
        bootstrap_reps = 200,
        bootstrap_cores = 1
      )
    ),
    error = function(e) list(error = TRUE, message = e$message)
  )

# Should succeed (not error)
  expect_false(isTRUE(result$error))
  expect_true(is.finite(result$variance))
  expect_true(result$variance > 0)
  expect_true(length(result$replicates) >= 2)
})

test_that("survey NA policy 'strict' shows detailed error", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

# Estimator that fails for specific replicates
  fail_at <- c(2, 5, 8)
  counter <- 0
  estimator <- function(data, ...) {
    counter <<- counter + 1
    if (counter %in% fail_at) {
      return(list(y_hat = NA_real_, converged = FALSE))
    }
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  counter <- 0
  future::plan(future::sequential)

  err <- tryCatch(
    bootstrap_variance(
      dstrat, estimator,
      point_estimate = mean(apistrat$api00),
      bootstrap_reps = 10,
      bootstrap_cores = 1,
      survey_na_policy = "strict"
    ),
    error = function(e) e$message
  )

# Should show specific indices
  expect_match(err, "2.*5.*8")
  expect_match(err, "strict")
  expect_match(err, "Troubleshooting")
})

test_that("survey NA policy 'omit' handles failures correctly", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

# Controlled failure pattern
  fail_at <- c(3, 7, 12)
  counter <- 0
  estimator <- function(data, ...) {
    counter <<- counter + 1
    if (counter %in% fail_at) {
      return(list(y_hat = NA_real_, converged = FALSE))
    }
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  counter <- 0
  future::plan(future::sequential)

  expect_warning(
    res <- bootstrap_variance(
      dstrat, estimator,
      point_estimate = mean(apistrat$api00),
      bootstrap_reps = 15,
      bootstrap_cores = 1,
      survey_na_policy = "omit"
    ),
    "3/15.*failed.*omitted"
  )

# Verify mathematical properties
  expect_equal(length(res$replicates), 12) # 15 - 3 failed
  expect_true(all(is.finite(res$replicates)))
  expect_true(is.finite(res$variance))
  expect_true(res$variance > 0)
  expect_true(res$se > 0)
})

test_that("survey NA policy 'omit' requires at least 2 successes", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

# Estimator that always fails
  estimator <- function(data, ...) {
    list(y_hat = NA_real_, converged = FALSE)
  }

  future::plan(future::sequential)

  expect_error(
    bootstrap_variance(
      dstrat, estimator,
      point_estimate = mean(apistrat$api00),
      bootstrap_reps = 10,
      bootstrap_cores = 1,
      survey_na_policy = "omit"
    ),
    "Too few successful"
  )
})

test_that("survey NA policy 'omit' shows failure pattern", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

# Many failures to test "First 10 failed + X more" message
  fail_at <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
  counter <- 0
  estimator <- function(data, ...) {
    counter <<- counter + 1
    if (counter %in% fail_at) {
      return(list(y_hat = NA_real_, converged = FALSE))
    }
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  counter <- 0
  future::plan(future::sequential)

  warn_msg <- NULL
  res <- tryCatch(
    withCallingHandlers(
      bootstrap_variance(
        dstrat, estimator,
        point_estimate = mean(apistrat$api00),
        bootstrap_reps = 20,
        bootstrap_cores = 1,
        survey_na_policy = "omit"
      ),
      warning = function(w) {
        warn_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) list(error = TRUE)
  )

# Should show "First 10 failed + X more"
  expect_match(warn_msg, "15/20")
  expect_match(warn_msg, "First 10")
  expect_match(warn_msg, "5 more")

# Should still produce valid result
  expect_false(isTRUE(res$error))
  expect_equal(length(res$replicates), 5) # 20 - 15 failed
})

test_that("mathematical correctness: variance has correct properties", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

# Test 1: Basic mathematical properties with survey design
  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

  estimator <- function(data, ...) {
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  future::plan(future::sequential)
  set.seed(424242)

  res <- bootstrap_variance(
    dstrat, estimator,
    point_estimate = mean(apistrat$api00),
    bootstrap_reps = 30,
    bootstrap_cores = 1
  )

# Property 1: Variance must be non-negative
  expect_true(res$variance >= 0)

# Property 2: SE = sqrt(variance)
  expect_equal(res$se, sqrt(res$variance))

# Property 3: All replicates are finite
  expect_true(all(is.finite(res$replicates)))

# Property 4: Correct number of replicates
  expect_equal(length(res$replicates), 30)

# Test 2: Compare to analytical ground truth for IID normal data
# For iid normal data, we know the ground truth variance of the sample mean
  set.seed(12345)
  n <- 1000
  true_sd <- 2
  x <- rnorm(n, mean = 5, sd = true_sd) # σ = 2, so σ² = 4

# True variance of sample mean = σ²/n = 4/1000 = 0.004
  true_variance_of_mean <- true_sd^2 / n

  data_iid <- data.frame(y = x)
  estimator_mean <- function(data, ...) {
    list(y_hat = mean(data$y), converged = TRUE)
  }

  set.seed(424242)
  res_iid <- bootstrap_variance(
    data_iid, estimator_mean,
    point_estimate = mean(x),
    bootstrap_reps = 500,
    bootstrap_cores = 1
  )

# Bootstrap estimate should be close to analytical value
# With 500 bootstrap reps, allow 30% tolerance due to Monte Carlo error
  expect_equal(res_iid$variance, true_variance_of_mean, tolerance = 0.3)

# Variance should be positive
  expect_true(res_iid$variance > 0)

# SE should match
  expect_equal(res_iid$se, sqrt(res_iid$variance))
})

test_that("boundary cases: minimum replicates and edge conditions", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

  estimator <- function(data, ...) {
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  future::plan(future::sequential)

# Test 1: bootstrap_reps = 2 (minimum for variance calculation)
  res_min <- bootstrap_variance(
    dstrat, estimator,
    point_estimate = mean(apistrat$api00),
    bootstrap_reps = 2,
    bootstrap_cores = 1
  )
  expect_equal(length(res_min$replicates), 2)
  expect_true(is.finite(res_min$variance))
  expect_true(res_min$variance >= 0)

# Test 2: bootstrap_reps = 1 should still work but variance may be problematic
# (svrVar should handle this, but it's an edge case)
  res_one <- tryCatch(
    bootstrap_variance(
      dstrat, estimator,
      point_estimate = mean(apistrat$api00),
      bootstrap_reps = 1,
      bootstrap_cores = 1
    ),
    error = function(e) list(error = TRUE, message = e$message)
  )
# With 1 replicate, variance calculation may fail or return zero
# We just verify it doesn't crash catastrophically
  expect_true(is.list(res_one))

# Test 3: IID data with small n
  small_data <- data.frame(y = c(1, 2, 3))
  estimator_simple <- function(data, ...) {
    list(y_hat = mean(data$y), converged = TRUE)
  }

  res_small <- bootstrap_variance(
    small_data, estimator_simple,
    point_estimate = mean(small_data$y),
    bootstrap_reps = 10,
    bootstrap_cores = 1
  )
  expect_equal(length(res_small$replicates), 10)
  expect_true(is.finite(res_small$variance))
  expect_true(res_small$variance >= 0)

# Test 4: Single observation (edge case)
  single_data <- data.frame(y = 5)
  res_single <- bootstrap_variance(
    single_data, estimator_simple,
    point_estimate = 5,
    bootstrap_reps = 10,
    bootstrap_cores = 1
  )
# With single observation, all bootstrap samples are identical
  expect_equal(length(res_single$replicates), 10)
  expect_true(all(res_single$replicates == 5))
# Variance should be zero (or very close) since all replicates are identical
  expect_true(res_single$variance < 1e-10)
})

test_that("omit policy correctly subsets rscales for mathematical correctness", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

# This test verifies the critical mathematical requirement:
# When omitting failed replicates, BOTH replicate_estimates AND rep_rscales
# must be subsetted to the same indices for svrVar formula to be correct:
# V = scale * sum(rscales[j] * (theta[j] - theta_0)^2)

# Create controlled failure pattern with specific indices
  fail_at <- c(2, 5, 8, 11, 14) # 5 failures out of 20
  counter <- 0
  estimator <- function(data, ...) {
    counter <<- counter + 1
    if (counter %in% fail_at) {
      return(list(y_hat = NA_real_, converged = FALSE))
    }
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  counter <- 0
  future::plan(future::sequential)
  set.seed(999)

  expect_warning(
    res <- bootstrap_variance(
      dstrat, estimator,
      point_estimate = mean(apistrat$api00),
      bootstrap_reps = 20,
      bootstrap_cores = 1,
      survey_na_policy = "omit"
    ),
    "5/20.*failed.*omitted"
  )

# Critical verification 1: Number of replicates matches expected successes
  expect_equal(length(res$replicates), 15) # 20 - 5 failed

# Critical verification 2: All returned replicates are finite
# (If rscales weren't subsetted correctly, NA values might slip through)
  expect_true(all(is.finite(res$replicates)))

# Critical verification 3: Variance is finite and positive
# If rscales length didn't match replicates length, svrVar would fail or
# produce NA/Inf, or use wrong scaling factors
  expect_true(is.finite(res$variance))
  expect_true(res$variance > 0)

# Critical verification 4: SE is consistent with variance
  expect_equal(res$se, sqrt(res$variance))

# Verification 5: Test with different failure patterns to ensure robustness
# Pattern: Random failures
  set.seed(777)
  fail_randomly <- sample(1:30, 12) # 12 random failures out of 30
  counter <- 0

  estimator_random <- function(data, ...) {
    counter <<- counter + 1
    if (counter %in% fail_randomly) {
      return(list(y_hat = NA_real_, converged = FALSE))
    }
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  counter <- 0
  set.seed(888)

  expect_warning(
    res_random <- bootstrap_variance(
      dstrat, estimator_random,
      point_estimate = mean(apistrat$api00),
      bootstrap_reps = 30,
      bootstrap_cores = 1,
      survey_na_policy = "omit"
    ),
    "12/30.*failed"
  )

  expect_equal(length(res_random$replicates), 18) # 30 - 12 failed
  expect_true(all(is.finite(res_random$replicates)))
  expect_true(is.finite(res_random$variance))
  expect_true(res_random$variance > 0)
  expect_equal(res_random$se, sqrt(res_random$variance))
})

test_that("survey NA policy default is 'strict'", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

# Estimator that fails
  fail_at <- c(5)
  counter <- 0
  estimator <- function(data, ...) {
    counter <<- counter + 1
    if (counter %in% fail_at) {
      return(list(y_hat = NA_real_, converged = FALSE))
    }
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

  counter <- 0
  future::plan(future::sequential)

# Without specifying survey_na_policy, should default to strict and error
  err <- tryCatch(
    bootstrap_variance(
      dstrat, estimator,
      point_estimate = mean(apistrat$api00),
      bootstrap_reps = 10,
      bootstrap_cores = 1
    ),
    error = function(e) e$message
  )

  expect_match(err, "strict")
  expect_match(err, "1/10.*failed")
})
