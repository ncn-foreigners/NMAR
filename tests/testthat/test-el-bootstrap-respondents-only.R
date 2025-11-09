test_that("EL bootstrap passes n_total through respondents-only data frames", {
  set.seed(123)
  df <- data.frame(
    Y_obs = rnorm(6),
    X1 = rnorm(6)
  )
  captured <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    bootstrap_variance = function(data,
                                  estimator_func,
                                  point_estimate,
                                  bootstrap_reps,
                                  formula,
                                  engine_args,
                                  ...) {
      captured$n_total <- engine_args$n_total
      list(se = 0.25, variance = 0.0625, replicates = rep(point_estimate, 2))
    },
    .package = "NMAR"
  )

  engine <- el_engine(
    n_total = nrow(df) + 2,
    variance_method = "bootstrap",
    bootstrap_reps = 3
  )

  fit <- nmar(Y_obs ~ 1, data = df, engine = engine)

  expect_equal(captured$n_total, nrow(df) + 2)
  expect_equal(se(fit), 0.25)
  expect_equal(sum(weights(fit, scale = "population")), nrow(df) + 2)
})

test_that("EL bootstrap passes n_total through respondents-only survey designs", {
  skip_if_not_installed("survey")

  set.seed(456)
  df <- data.frame(
    Y_obs = rnorm(6),
    wt = sample(1:4, 6, replace = TRUE)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~wt)
  pop_total <- sum(df$wt)
  captured <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    bootstrap_variance = function(data,
                                  estimator_func,
                                  point_estimate,
                                  bootstrap_reps,
                                  formula,
                                  engine_args,
                                  ...) {
      captured$n_total <- engine_args$n_total
      list(se = 0.5, variance = 0.25, replicates = rep(point_estimate, 2))
    },
    .package = "NMAR"
  )

  engine <- el_engine(
    n_total = pop_total,
    variance_method = "bootstrap",
    bootstrap_reps = 3
  )

  fit <- nmar(Y_obs ~ 1, data = design, engine = engine)

  expect_equal(captured$n_total, pop_total)
  expect_equal(se(fit), 0.5)
  expect_equal(sum(weights(fit, scale = "population")), pop_total)
})
