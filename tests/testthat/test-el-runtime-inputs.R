test_that("el runtime inputs create unique delta names", {
  df <- data.frame(`..nmar_delta..` = 1:5, Y = c(1, 2, NA, 4, 5), X = rnorm(5))
  engine <- el_engine(variance_method = "none")
  spec <- NMAR:::parse_nmar_spec(Y ~ X, df)
  traits <- NMAR:::engine_traits(engine)
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task, standardize = engine$standardize, auxiliary_means = engine$auxiliary_means)
  runtime <- NMAR:::el_build_runtime_inputs(
    data = design$data,
    design_info = design,
    auxiliary_means = engine$auxiliary_means,
    n_total = engine$n_total,
    require_na = TRUE,
    context = "data.frame"
  )
  expect_true(any(grepl("^\\.\\.nmar_delta", names(runtime$data))))
  expect_true(runtime$response_var %in% names(runtime$data))
  expect_true(all(runtime$data[[runtime$response_var]] %in% c(0, 1)))
})

test_that("el runtime inputs enforce auxiliary mean guard on respondents-only data", {
  df <- data.frame(Y = c(1, 2, 3), X = rnorm(3))
  engine <- el_engine(variance_method = "none", n_total = 50)
  spec <- NMAR:::parse_nmar_spec(Y ~ X, df)
  traits <- NMAR:::engine_traits(engine)
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task, standardize = engine$standardize, auxiliary_means = engine$auxiliary_means)
  expect_error(
    NMAR:::el_build_runtime_inputs(
      data = design$data,
      design_info = design,
      auxiliary_means = NULL,
      n_total = 50,
      require_na = FALSE,
      context = "data.frame"
    ),
    "auxiliary_means",
    fixed = FALSE
  )
})

test_that("el runtime inputs require n_total for respondents-only data", {
  df <- data.frame(Y = c(1, 2, 3), X = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(Y ~ X, df)
  traits <- utils::modifyList(NMAR:::NMAR_DEFAULT_TRAITS, list(allow_respondents_only = TRUE))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task, standardize = TRUE, auxiliary_means = c(X = 0))
  expect_error(
    NMAR:::el_build_runtime_inputs(
      data = design$data,
      design_info = design,
      auxiliary_means = c(X = 0),
      n_total = NULL,
      require_na = FALSE,
      context = "data.frame"
    ),
    "n_total",
    fixed = FALSE
  )
})
