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

suppress_bootstrap_assumption_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("injects replicate analysis weights", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
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

test_that("bootstrap_settings cannot override design/replicates", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  df <- data.frame(y = c(1, 2, 3, 4), w = c(1, 2, 1, 3))
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)

  expect_error(
    bootstrap_variance(design,
      estimator_func = bootstrap_dummy_estimator,
      point_estimate = 1,
      bootstrap_reps = 3,
      bootstrap_settings = list(replicates = 999)
    ),
    "must not include `design` or `replicates`",
    fixed = TRUE
  )
})

test_that("bootstrap_mse requires a finite point_estimate", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  df <- data.frame(y = c(1, 2, 3, 4), w = c(1, 2, 1, 3))
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)

  expect_error(
    bootstrap_variance(design,
      estimator_func = bootstrap_dummy_estimator,
      point_estimate = NA_real_,
      bootstrap_reps = 3,
      bootstrap_mse = TRUE
    ),
    "point_estimate",
    fixed = FALSE
  )
})

test_that("survey-only bootstrap options are ignored for iid data", {
  skip_if_not_installed("future.apply")

  set.seed(123)
  df <- data.frame(y = c(1, 2, 3, 4))
  res <- bootstrap_variance(df,
    estimator_func = bootstrap_dummy_estimator,
    point_estimate = mean(df$y),
    bootstrap_reps = 5,
    bootstrap_type = "Preston",
    bootstrap_mse = TRUE,
    bootstrap_settings = list(type = "Rao-Wu-Yue-Beaumont"),
    bootstrap_options = list(mse = FALSE),
    survey_na_policy = "omit"
  )
  expect_length(res$replicates, 5)
  expect_true(is.finite(res$variance))
})

test_that("iid-only resample_guard is ignored for survey designs", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  df <- data.frame(y = c(1, 2, 3, 4), w = c(1, 2, 1, 3))
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  point_est <- sum(df$y * df$w) / sum(df$w)

  set.seed(123)
  res <- bootstrap_variance(design,
    estimator_func = bootstrap_dummy_estimator,
    point_estimate = point_est,
    bootstrap_reps = 5,
    resample_guard = function(indices, data) FALSE
  )
  expect_length(res$replicates, 5)
  expect_true(is.finite(res$variance))
})

test_that("bootstrap falls back to sequential lapply when future.apply is unavailable", {
  testthat::local_mocked_bindings(
    nmar_has_future_apply = function() FALSE,
    .package = "NMAR"
  )
  testthat::local_mocked_bindings(
    nmar_future_workers = function() 1L,
    .package = "NMAR"
  )
  opt <- "NMAR.bootstrap.warned_no_future_apply"
  old_opt <- getOption(opt)
  options(setNames(list(FALSE), opt))
  on.exit(options(setNames(list(old_opt), opt)), add = TRUE)

  set.seed(123)
  df <- data.frame(y = c(1, 2, 3, 4))
  point_est <- mean(df$y)

  res1 <- bootstrap_variance(df,
    estimator_func = bootstrap_dummy_estimator,
    point_estimate = point_est,
    bootstrap_reps = 5
  )
  expect_length(res1$replicates, 5)

  expect_silent(
    res2 <- bootstrap_variance(df,
      estimator_func = bootstrap_dummy_estimator,
      point_estimate = point_est,
      bootstrap_reps = 5
    )
  )
  expect_length(res2$replicates, 5)

  testthat::local_mocked_bindings(
    nmar_future_workers = function() 2L,
    .package = "NMAR"
  )
  expect_warning(
    bootstrap_variance(df,
      estimator_func = bootstrap_dummy_estimator,
      point_estimate = point_est,
      bootstrap_reps = 5
    ),
    "future\\.apply.*not installed|future\\.apply.*Install",
    fixed = FALSE
  )
})

test_that("survey bootstrap falls back to sequential lapply when future.apply is unavailable", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  testthat::local_mocked_bindings(
    nmar_has_future_apply = function() FALSE,
    .package = "NMAR"
  )
  testthat::local_mocked_bindings(
    nmar_future_workers = function() 1L,
    .package = "NMAR"
  )
  opt <- "NMAR.bootstrap.warned_no_future_apply"
  old_opt <- getOption(opt)
  options(setNames(list(FALSE), opt))
  on.exit(options(setNames(list(old_opt), opt)), add = TRUE)

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

  testthat::local_mocked_bindings(
    nmar_future_workers = function() 2L,
    .package = "NMAR"
  )
  expect_warning(
    bootstrap_variance(design,
      estimator_func = bootstrap_dummy_estimator,
      point_estimate = point_est,
      bootstrap_reps = 5
    ),
    "future\\.apply.*not installed|future\\.apply.*Install",
    fixed = FALSE
  )
})

test_that("IID bootstrap is reproducible across future backends", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  data <- data.frame(x = rnorm(100), y = rnorm(100))
  estimator <- function(data, ...) {
    list(y_hat = mean(data$y), converged = TRUE)
  }

# Force future.apply for all backends so the RNG stream matches.
  old_opt <- getOption("nmar.bootstrap_apply", NULL)
  on.exit(options(nmar.bootstrap_apply = old_opt), add = TRUE)
  options(nmar.bootstrap_apply = "future")

# Sequential backend
  future::plan(future::sequential)
  set.seed(424242)
  res_seq <- bootstrap_variance(
    data, estimator, point_estimate = 0, bootstrap_reps = 50
  )

# Multisession backend (2 workers)
  skip_on_cran() # Multisession can be unreliable on CRAN
  tryCatch(
    future::plan(future::multisession, workers = 2),
    error = function(e) skip(paste("future::multisession not available:", conditionMessage(e)))
  )
  set.seed(424242)
  res_par2 <- bootstrap_variance(
    data, estimator, point_estimate = 0, bootstrap_reps = 50
  )

# Multisession backend (4 workers)
  tryCatch(
    future::plan(future::multisession, workers = 4),
    error = function(e) skip(paste("future::multisession not available:", conditionMessage(e)))
  )
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
  res_seq1 <- suppress_bootstrap_assumption_warning(
    bootstrap_variance(
      dstrat, estimator, point_estimate = mean(apistrat$api00),
      bootstrap_reps = 50
    )
  )

# Sequential with same seed should be identical
  future::plan(future::sequential)
  set.seed(999)
  res_seq2 <- suppress_bootstrap_assumption_warning(
    bootstrap_variance(
      dstrat, estimator, point_estimate = mean(apistrat$api00),
      bootstrap_reps = 50
    )
  )

  future::plan(future::sequential)

# Test reproducibility within sequential mode
  expect_identical(res_seq1$replicates, res_seq2$replicates,
                   label = "Survey sequential reproducibility")
  expect_equal(res_seq1$variance, res_seq2$variance,
               label = "Survey sequential variance reproducibility")
})
