test_that("EL print method shows convergence message on failure objects", {
  res <- NMAR:::new_nmar_result(
    estimate = NA_real_,
    estimate_name = "Y_miss",
    se = NA_real_,
    converged = FALSE,
    model = list(coefficients = NULL, vcov = NULL),
    weights_info = list(values = numeric(0), trimmed_fraction = NA_real_),
    sample = list(n_total = 10, n_respondents = 5, is_survey = FALSE, design = NULL),
    inference = list(variance_method = "none", df = NA_real_, message = "failed"),
    diagnostics = list(message = "solver failed"),
    meta = list(engine_name = "empirical_likelihood", formula = Y_miss ~ 1),
    extra = list(),
    class = "nmar_result_el"
  )

  out <- capture.output(print(res))
  expect_true(any(grepl("Method:", out)))
  expect_true(any(grepl("Convergence message:", out)))
})

test_that("EL print and summary printing cover vcov and df branches", {
  res <- NMAR:::new_nmar_result(
    estimate = 1,
    estimate_name = "Y_miss",
    se = 1,
    converged = TRUE,
    model = list(
      coefficients = c("(Intercept)" = 0.1, "Y_miss" = 0.2),
      vcov = diag(c(0.04, 0.09))
    ),
    weights_info = list(values = rep(1, 5), trimmed_fraction = 0),
    sample = list(n_total = 10, n_respondents = 5, is_survey = FALSE, design = NULL),
    inference = list(variance_method = "none", df = 10, message = NA_character_),
    diagnostics = list(max_equation_residual = 1e-8, constraint_sum_W = 1e-8, constraint_sum_aux = c(X = 0)),
    meta = list(engine_name = "empirical_likelihood", formula = Y_miss ~ 1),
    extra = list(),
    class = "nmar_result_el"
  )

  out_print <- capture.output(print(res))
  expect_true(any(grepl("Max equation residual:", out_print)))
  expect_true(any(grepl("Constraint sum \\(W\\):", out_print)))
  expect_true(any(grepl("Constraint sums \\(aux\\):", out_print)))

  old <- options(nmar.show_call = TRUE)
  on.exit(options(old), add = TRUE)
  out_sum <- capture.output(print(summary(res)))
  expect_true(any(grepl("Missingness-model coefficients", out_sum)))
  expect_true(any(grepl("t value", out_sum)))
  expect_true(any(grepl("Pr\\(>\\|t\\|\\)", out_sum)))

  res2 <- res
  res2$inference$df <- NA_real_
  out_sum2 <- capture.output(print(summary(res2)))
  expect_true(any(grepl("z value", out_sum2)))
  expect_true(any(grepl("Pr\\(>\\|z\\|\\)", out_sum2)))
})
