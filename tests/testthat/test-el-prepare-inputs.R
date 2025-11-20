test_that("delta column name avoids collisions and is present", {
  df <- data.frame(`..nmar_delta..` = 1:5, Y = c(1, 2, NA, 4, 5), X = rnorm(5))
  out <- el_make_delta_column_name(df, outcome_var = "Y", respondent_mask = !is.na(df$Y))
  expect_true(out$delta_column_name != "..nmar_delta..")
  expect_true(out$delta_column_name %in% names(out$data))
})

test_that("el_process_design expands dot notation via Formula", {
  df <- data.frame(
    Y_miss = c(1, NA, 0, 1),
    X1 = rnorm(4),
    X2 = rnorm(4),
    Z = rnorm(4)
  )
  design <- el_process_design(Y_miss ~ . | ., df)
  expect_setequal(colnames(design$auxiliary_design_full), c("X1", "X2", "Z"))
  expect_true(all(c("(Intercept)", "Y_miss", "X1", "X2", "Z") %in% colnames(design$missingness_design)))
})

test_that("el_process_design accepts transformed outcomes and tracks source column", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  design <- el_process_design(log(Y_miss) ~ X, df)
  expect_identical(design$outcome_var, "log(Y_miss)")
  expect_identical(design$outcome_source, "Y_miss")
  expect_true("log(Y_miss)" %in% colnames(design$missingness_design))
})

test_that("LHS expressions referencing multiple outcome variables are rejected", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    Y2 = rnorm(4),
    X = rnorm(4)
  )
  expect_error(
    el_process_design(I(Y_miss + Y2) ~ X, df),
    regexp = "Left-hand side may reference only one outcome variable",
    fixed = FALSE
  )
})

test_that("transforms that introduce NA for respondents are rejected", {
  df <- data.frame(
    Y_miss = c(1, NA, -2, NA),
    X = rnorm(4)
  )
  expect_error(
    el_process_design(log(Y_miss) ~ X, df),
    regexp = "produced NA/NaN",
    fixed = FALSE
  )
})

test_that("response intercept is retained even when formula uses +0", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    Z = rnorm(4)
  )
  expect_warning(
    design <- el_process_design(Y_miss ~ 1 | Z + 0, df),
    "Missingness-model intercept is required",
    fixed = FALSE
  )
  expect_identical(colnames(design$missingness_design)[1], "(Intercept)")
  expect_true(all(design$missingness_design[, "(Intercept)"] == 1))
  expect_true("Z" %in% colnames(design$missingness_design))
  expect_equal(ncol(design$missingness_design), 3)
})

test_that("respondents-only data frame requires n_total", {
  df <- data.frame(Y_miss = 1:4)
  eng <- el_engine(variance_method = "none")
  expect_error(
    nmar(Y_miss ~ 1, data = df, engine = eng),
    "Respondents-only data detected",
    fixed = FALSE
  )
  eng_with_total <- el_engine(variance_method = "none", n_total = 100)
  expect_s3_class(
    nmar(Y_miss ~ 1, data = df, engine = eng_with_total),
    "nmar_result_el"
  )
})

test_that("dot expansion drops outcome-derived auxiliary terms", {
  set.seed(99)
  df <- data.frame(
    Y_miss = c(rnorm(3), NA),
    X = rnorm(4),
    F = factor(c("a", "b", "a", "b"))
  )
  design <- el_process_design(Y_miss ~ . | ., df)
  expect_false("Y_miss" %in% colnames(design$auxiliary_design_full))
  expect_false(anyNA(design$auxiliary_design_full))
})

test_that("el_process_design forbids outcome in auxiliary constraints", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4),
    Z = rnorm(4)
  )
  expect_error(
    el_process_design(Y_miss ~ I(Y_miss^2) + X | Z, df),
    "outcome cannot appear",
    fixed = FALSE
  )
})

test_that("intercept-only auxiliaries are ignored silently", {
  df <- data.frame(
    Y_miss = c(1, NA, 3, NA),
    Z = rnorm(4)
  )
  design <- el_process_design(Y_miss ~ 1 | Z, df)
  expect_equal(ncol(design$auxiliary_design_full), 0)
})

test_that("explicit intercept plus auxiliaries triggers a warning", {
  df <- data.frame(
    Y_miss = c(1, NA, 3, NA),
    X = rnorm(4)
  )
  expect_warning(
    design <- el_process_design(Y_miss ~ 1 + X, df),
    "Auxiliary intercepts are ignored",
    fixed = FALSE
  )
  expect_equal(colnames(design$auxiliary_design_full), "X")
})

test_that("auxiliary '-1' or '+0' are no-ops (no warning, no intercept)", {
  set.seed(123)
  df <- data.frame(
    Y_miss = c(1, NA, 3, NA),
    X = rnorm(4),
    Z = rnorm(4)
  )
# '-1' should remove any implicit intercept silently
  design_m1 <- el_process_design(Y_miss ~ X - 1 | Z, df)
  expect_false("(Intercept)" %in% colnames(design_m1$auxiliary_design_full))
  expect_equal(colnames(design_m1$auxiliary_design_full), "X")

# '+0' should behave identically and silently
  design_p0 <- el_process_design(Y_miss ~ X + 0 | Z, df)
  expect_false("(Intercept)" %in% colnames(design_p0$auxiliary_design_full))
  expect_equal(colnames(design_p0$auxiliary_design_full), "X")
})

test_that("el_process_design rejects formulas with more than two RHS partitions", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X1 = rnorm(4),
    X2 = rnorm(4),
    X3 = rnorm(4)
  )
  expect_error(
    el_process_design(Formula::Formula(Y_miss ~ X1 | X2 | X3), df),
    "at most two RHS partitions",
    fixed = FALSE
  )
})

test_that("auxiliary constraints cannot reference the outcome", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  expect_error(
    el_process_design(Y_miss ~ Y_miss + X, df),
    regexp = "outcome cannot appear",
    fixed = FALSE
  )
})

test_that("constant auxiliary columns created via transforms are rejected", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  expect_error(
    el_process_design(Y_miss ~ I(0 * X + 1) + X, df),
    regexp = "Auxiliary covariate .* zero variance",
    fixed = FALSE
  )
})

test_that("el_process_design handles language objects coerced to formula", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  lang_formula <- as.call(list(as.name("~"), quote(Y_miss), quote(X)))
  parsed <- el_process_design(lang_formula, df)
  expect_identical(parsed$outcome_var, "Y_miss")
  expect_identical(parsed$outcome_source, "Y_miss")
  expect_true("(Intercept)" %in% colnames(parsed$missingness_design))
})

test_that("response-only formulas drop RHS2 entirely", {
  set.seed(55)
  df <- data.frame(
    Y_miss = c(rnorm(4), NA),
    X = rnorm(5)
  )
  parsed <- el_process_design(Y_miss ~ X, df)
  expect_identical(colnames(parsed$missingness_design), c("(Intercept)", "Y_miss"))
})

test_that("offset() usage is rejected on auxiliary and response partitions", {
  set.seed(123)
  df <- data.frame(
    Y_miss = c(rnorm(4), NA),
    X = rnorm(5),
    Z = rnorm(5)
  )
  expect_error(
    el_process_design(Y_miss ~ offset(X) + Z, df),
    "Offsets .* auxiliary predictors",
    fixed = FALSE
  )
  expect_error(
    el_process_design(Y_miss ~ Z | offset(X), df),
    "Offsets .* response predictors",
    fixed = FALSE
  )
})

test_that("auxiliary predictors with zero variance among respondents trigger an error", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    A = c(5, 2, 5, 3)
  )
  expect_error(
    el_process_design(Y_miss ~ A, df),
    regexp = "Auxiliary covariate .* zero variance",
    fixed = FALSE
  )
})

test_that("missingness predictors with zero variance among respondents only warn", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4),
    Z = c(4, 1, 4, 2)
  )
  expect_warning(
    design <- el_process_design(Y_miss ~ X | Z, df),
    regexp = "Missingness-model predictor .* zero variance",
    fixed = FALSE
  )
  expect_s3_class(design, "el_design_spec")
})
