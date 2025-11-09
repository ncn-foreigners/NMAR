make_el_test_data <- function(n = 100) {
  set.seed(123)
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.5 * X + 0.3 * Z + rnorm(n)
  p <- plogis(-0.3 + 0.4 * scale(Y)[, 1] + 0.2 * Z)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_
  df
}

test_that("el engine accepts valid inputs", {
  df <- make_el_test_data()
  el_base <- nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = el_engine(auxiliary_means = c(X = 0), variance_method = 'none')
  )
  expect_s3_class(el_base, "nmar_result_el")
  expect_type(el_base$converged, "logical")
  expect_equal(formula(el_base), Y_miss ~ X)

  el_resp_extra <- nmar(
    formula = Y_miss ~ X | Z,
    data = df,
    engine = el_engine(auxiliary_means = c(X = 0), variance_method = 'none')
  )
  expect_s3_class(el_resp_extra, "nmar_result_el")
  expect_type(el_resp_extra$converged, "logical")
  expect_equal(formula(el_resp_extra), Y_miss ~ X | Z)

  skip_if_not_installed("survey")
  design <- survey::svydesign(ids = ~1, data = df, weights = ~1)
  el_svy <- nmar(
    formula = Y_miss ~ X | Z,
    data = design,
    engine = el_engine(auxiliary_means = c(X = 0), variance_method = 'none')
  )
  expect_s3_class(el_svy, "nmar_result_el")
  expect_true(is.numeric(stats::weights(el_svy)))
})

test_that("el engine rejects faulty inputs", {
  base_df <- make_el_test_data()
  good_engine <- el_engine(auxiliary_means = c(X = 0), variance_method = 'none')

  bad_cases <- list(
    list(
      name = "missing outcome column",
      data = base_df[, c("X", "Z")],
      formula = Y_miss ~ X,
      response_predictors = NULL,
      message = "Outcome variable"
    ),
    list(
      name = "outcome without NA",
      data = transform(base_df, Y_miss = ifelse(is.na(Y_miss), 0, Y_miss)),
      formula = Y_miss ~ X,
      response_predictors = NULL,
      message = "must contain NA"
    ),
    list(
      name = "missing covariate",
      data = base_df,
      formula = Y_miss ~ W,
      response_predictors = NULL,
      message = "Variables not found"
    ),
    list(
      name = "missing response predictor",
      data = base_df,
      formula = Y_miss ~ X | W,
      message = "Variables not found"
    ),
    list(
      name = "non numeric outcome",
      data = transform(base_df, Y_miss = ifelse(is.na(Y_miss), NA, as.character(Y_miss))),
      formula = Y_miss ~ X,
      response_predictors = NULL,
      message = "must be numeric"
    ),
    list(
      name = "NA in covariate",
      data = transform(base_df, X = replace(X, 1, NA_real_)),
      formula = Y_miss ~ X,
      response_predictors = NULL,
      message = "contains NA values"
    )
  )

  for (case in bad_cases) {
    expect_error(
      nmar(
        formula = case$formula,
        data = case$data,
        engine = good_engine
      ),
      regexp = case$message,
      info = case$name
    )
  }
})

test_that("auxiliary means align with model-matrix columns", {
  df <- make_el_test_data(40)
  spec <- NMAR:::parse_nmar_spec(Y_miss ~ X, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)

  runtime <- NMAR:::el_build_runtime_inputs(
    data = design$data,
    design_info = design,
    auxiliary_means = NULL,
    n_total = NULL,
    require_na = TRUE,
    context = "data.frame"
  )
  precomputed <- NMAR:::el_build_precomputed_design(
    design_matrices = design$design_matrices,
    estimation_data = design$data,
    outcome_var = runtime$outcome_name,
    respondent_indices = runtime$observed_indices
  )

  aux_res <- NMAR:::el_resolve_auxiliaries(
    full_data = design$data,
    respondent_data = design$data[runtime$observed_indices, , drop = FALSE],
    aux_formula = runtime$internal_formula$auxiliary,
    auxiliary_means = c(X = 0, bogus = 5),
    precomputed_resp = precomputed$auxiliary_resp,
    precomputed_full = precomputed$auxiliary_full
  )

  expect_equal(colnames(aux_res$matrix), "X")
  expect_equal(names(aux_res$means), "X")
  expect_true(aux_res$has_aux)
})
