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
  res <- NMAR:::bootstrap_variance(df,
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
  res <- NMAR:::bootstrap_variance(design,
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

  res <- suppressWarnings(NMAR:::bootstrap_variance(design,
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

test_that("bootstrap_variance errors when design call not svydesign", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")
  df <- data.frame(y = c(1, 2, 3, 4), w = c(1, 2, 1, 3))
  base_design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  subset_design <- subset(base_design, y > 1)
  expect_error(
    NMAR:::bootstrap_variance(subset_design,
      estimator_func = bootstrap_dummy_estimator,
      point_estimate = 1,
      bootstrap_reps = 3
    ),
    "survey::svydesign",
    fixed = FALSE
  )
})

test_that("bootstrap_variance errors for svydesigns that use probs/pps", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")
  df <- data.frame(y = c(1, 2, 3, 4), pi = c(0.2, 0.3, 0.25, 0.25))
  design <- survey::svydesign(ids = ~1, probs = ~pi, data = df)
  expect_error(
    NMAR:::bootstrap_variance(design,
      estimator_func = bootstrap_dummy_estimator,
      point_estimate = 1,
      bootstrap_reps = 3
    ),
    "probabilities/PPS arguments",
    fixed = FALSE
  )
})

test_that("nmar_adjust_bootstrap_scale rescales omit-policy replicates", {
  scalar_scale <- 1 / 50
  adj_scalar <- NMAR:::nmar_adjust_bootstrap_scale(scalar_scale, total_reps = 50, keep_idx = 1:10)
  expect_equal(adj_scalar, scalar_scale * 5)

  vector_scale <- rep(0.1, 4)
  adj_vector <- NMAR:::nmar_adjust_bootstrap_scale(vector_scale, total_reps = 4, keep_idx = c(1, 3))
  expect_equal(adj_vector, vector_scale[c(1, 3)])

  expect_error(
    NMAR:::nmar_adjust_bootstrap_scale(rep(0.1, 3), total_reps = 4, keep_idx = 1:2),
    "Replicate scale length",
    fixed = FALSE
  )
})
