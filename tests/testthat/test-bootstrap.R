bootstrap_dummy_estimator <- function(data, on_failure = "return") {
  if (inherits(data, "survey.design")) {
    y <- data$variables$y
    w <- stats::weights(data)
    est <- sum(y * w) / sum(w)
  } else {
    est <- mean(data$y)
  }
  structure(
    list(
      y_hat = est,
      estimate_name = "y",
      se = NA_real_,
      converged = TRUE,
      sample = list(is_survey = inherits(data, "survey.design")),
      inference = list(),
      model = list(),
      weights_info = list(values = NULL, trimmed_fraction = NA_real_)
    ),
    class = c("bootstrap_dummy_result", "nmar_result")
  )
}
test_that("bootstrap_variance handles iid data", {
  set.seed(123)
  df <- data.frame(y = c(1, 2, 3, 4))
  point_est <- mean(df$y)
  res <- bootstrap_variance(df,
    estimator_func = bootstrap_dummy_estimator,
    point_estimate = point_est,
    bootstrap_reps = 10
  )
  expect_length(res$replicates, 10)
  expect_true(is.finite(res$se))
})

test_that("bootstrap_variance handles survey designs", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")
  set.seed(123)
  df <- data.frame(y = c(1, 2, 3, 4), w = c(1, 2, 1, 3))
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  point_est <- sum(df$y * df$w) / sum(df$w)
  res <- bootstrap_variance(design,
    estimator_func = bootstrap_dummy_estimator,
    point_estimate = point_est,
    bootstrap_reps = 5
  )
  expect_length(res$replicates, 5)
  expect_true(is.finite(res$se))
})

test_that("bootstrap_variance forwards settings to svrep", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")
  df <- data.frame(y = c(1, 2, 3, 4), w = c(1, 2, 1, 3))
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  point_est <- sum(df$y * df$w) / sum(df$w)
  original <- svrep::as_bootstrap_design
  observed <- new.env(parent = emptyenv())
  observed$type <- NULL
  observed$mse <- NULL

  testthat::local_mocked_bindings(
    as_bootstrap_design = function(design, replicates, ...) {
      args <- list(...)
      observed$type <- args$type
      observed$mse <- args$mse
      do.call(original, c(list(design = design, replicates = replicates), args))
    },
    .package = "svrep"
  )

  res <- suppressWarnings(bootstrap_variance(design,
    estimator_func = bootstrap_dummy_estimator,
    point_estimate = point_est,
    bootstrap_reps = 4,
    bootstrap_type = "Preston",
    bootstrap_mse = TRUE
  ))
  expect_length(res$replicates, 4)
  expect_equal(observed$type, "Preston")
  expect_true(isTRUE(observed$mse))
})

test_that("bootstrap_variance supports subset() and update() survey designs", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")
  df <- data.frame(y = c(1, 2, 3, 4), w = c(1, 2, 1, 3))
  base_design <- survey::svydesign(ids = ~1, data = df, weights = ~w)

  subset_design <- subset(base_design, y > 1)
  subset_point <- sum(subset_design$variables$y * stats::weights(subset_design)) / sum(stats::weights(subset_design))
  res_subset <- bootstrap_variance(subset_design,
    estimator_func = bootstrap_dummy_estimator,
    point_estimate = subset_point,
    bootstrap_reps = 3
  )
  expect_length(res_subset$replicates, 3)
  expect_true(is.finite(res_subset$se))

  update_design <- update(base_design, y2 = 2 * y)
  update_point <- sum(update_design$variables$y * stats::weights(update_design)) / sum(stats::weights(update_design))
  res_update <- bootstrap_variance(update_design,
    estimator_func = bootstrap_dummy_estimator,
    point_estimate = update_point,
    bootstrap_reps = 3
  )
  expect_length(res_update$replicates, 3)
  expect_true(is.finite(res_update$se))
})

test_that("bootstrap_variance errors on calibrated survey designs", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")
  df <- data.frame(y = c(1, 2, 3, 4), w = c(1, 2, 1, 3))
  base_design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  calibrated <- survey::calibrate(base_design, formula = ~1, population = c("(Intercept)" = 10))

  expect_error(
    bootstrap_variance(calibrated,
      estimator_func = bootstrap_dummy_estimator,
      point_estimate = 1,
      bootstrap_reps = 3
    ),
    regexp = "calibrated|post-stratified",
    fixed = FALSE
  )
})
