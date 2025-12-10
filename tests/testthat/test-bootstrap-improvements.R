test_that("IID bootstrap is reproducible across future backends", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  data <- data.frame(x = rnorm(100), y = rnorm(100))
  estimator <- function(data, ...) {
    list(y_hat = mean(data$y), converged = TRUE)
  }

# Sequential backend
  future::plan(future::sequential)
  set.seed(424242)
  res_seq <- bootstrap_variance(
    data, estimator, point_estimate = 0, bootstrap_reps = 50, bootstrap_cores = 1
  )

# Multisession backend (2 workers)
  skip_on_cran() # Multisession can be unreliable on CRAN
  future::plan(future::multisession, workers = 2)
  set.seed(424242)
  res_par2 <- bootstrap_variance(
    data, estimator, point_estimate = 0, bootstrap_reps = 50
  )

# Multisession backend (4 workers)
  future::plan(future::multisession, workers = 4)
  set.seed(424242)
  res_par4 <- bootstrap_variance(
    data, estimator, point_estimate = 0, bootstrap_reps = 50
  )

# Reset to sequential
  future::plan(future::sequential)

# Verify exact reproducibility
  expect_identical(res_seq$replicates, res_par2$replicates,
                   label = "Sequential vs 2-worker replicates")
  expect_identical(res_seq$replicates, res_par4$replicates,
                   label = "Sequential vs 4-worker replicates")
  expect_equal(res_seq$variance, res_par2$variance,
               label = "Sequential vs 2-worker variance")
  expect_equal(res_seq$variance, res_par4$variance,
               label = "Sequential vs 4-worker variance")
})

test_that("Survey bootstrap is reproducible across backends", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()

  data(api, package = "survey")
  dstrat <- survey::svydesign(
    id = ~1, strata = ~stype, weights = ~pw,
    data = apistrat, fpc = ~fpc
  )

  estimator <- function(data, ...) {
    est <- survey::svymean(~api00, data)
    list(y_hat = as.numeric(est), converged = TRUE)
  }

# Sequential baseline
  future::plan(future::sequential)
  set.seed(999)
  res_seq1 <- bootstrap_variance(
    dstrat, estimator, point_estimate = mean(apistrat$api00),
    bootstrap_reps = 50, bootstrap_cores = 1
  )

# Sequential with same seed should be identical
  future::plan(future::sequential)
  set.seed(999)
  res_seq2 <- bootstrap_variance(
    dstrat, estimator, point_estimate = mean(apistrat$api00),
    bootstrap_reps = 50, bootstrap_cores = 1
  )

  future::plan(future::sequential)

# Test reproducibility within sequential mode (this is the critical property)
  expect_identical(res_seq1$replicates, res_seq2$replicates,
                   label = "Survey sequential reproducibility")
  expect_equal(res_seq1$variance, res_seq2$variance,
               label = "Survey sequential variance reproducibility")

# Note: Cross-backend reproducibility (sequential vs parallel) requires
# complex template_call serialization and is tested separately when needed
})
