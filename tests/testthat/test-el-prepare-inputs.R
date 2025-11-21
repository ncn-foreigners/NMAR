test_that("delta column name avoids collisions and is present", {
  df <- data.frame(`..nmar_delta..` = 1:5, Y = c(1, 2, NA, 4, 5), X = rnorm(5))
  spec <- el_prepare_inputs(Y ~ X, df)
  cols <- names(spec$analysis_data)
# Original column should remain
  expect_true("..nmar_delta.." %in% cols)
# A new delta column should have been added with a different name
  delta_like <- grep("^\\.\\.nmar_delta\\.\\.", cols, value = TRUE)
  expect_true(any(delta_like != "..nmar_delta.."))
})

test_that("el_prepare_inputs expands dot notation via Formula", {
  df <- data.frame(
    Y_miss = c(1, NA, 0, 1),
    X1 = rnorm(4),
    X2 = rnorm(4),
    Z = rnorm(4)
  )
  design <- el_prepare_inputs(Y_miss ~ . | ., df)
  expect_setequal(colnames(design$aux_design_full), c("X1", "X2", "Z"))
  expect_true(all(c("(Intercept)", "Y_miss", "X1", "X2", "Z") %in% colnames(design$missingness_design)))
})

test_that("el_prepare_inputs accepts transformed outcomes and tracks source column", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  design <- el_prepare_inputs(log(Y_miss) ~ X, df)
  expect_identical(design$outcome_expr, "log(Y_miss)")
  expect_true("log(Y_miss)" %in% colnames(design$missingness_design))
})

test_that("LHS expressions referencing multiple outcome variables are rejected", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    Y2 = rnorm(4),
    X = rnorm(4)
  )
  expect_error(
    el_prepare_inputs(I(Y_miss + Y2) ~ X, df),
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
    el_prepare_inputs(log(Y_miss) ~ X, df),
    regexp = "produced NA/NaN",
    fixed = FALSE
  )
})

test_that("respondent mask is based on raw outcome, not transformed LHS", {
  df <- data.frame(
    Y_miss = c(NA, 2, NA, 4),
    X = rnorm(4)
  )
# transform maps NA to 0 but mask should still treat NA rows as nonrespondents
  design <- el_prepare_inputs(ifelse(is.na(Y_miss), 0, Y_miss) ~ X, df)
  expect_identical(design$respondent_mask, c(FALSE, TRUE, FALSE, TRUE))
})

test_that("outcome terms on RHS2 are rejected", {
  df <- data.frame(
    Y_miss = c(1, 2, NA, 4),
    X = rnorm(4),
    Z = rnorm(4)
  )
  expect_error(
    el_prepare_inputs(log(Y_miss) ~ X | Y_miss + Z, df),
    regexp = "Outcome cannot appear in missingness predictors",
    fixed = FALSE
  )
})

test_that("logical and two-level factor outcomes are coerced to numeric with a warning", {
  set.seed(42)
  df <- data.frame(
    Y_logical = c(TRUE, FALSE, TRUE, NA),
    Y_factor = factor(c("a", "b", "a", NA)),
    X = rnorm(4)
  )
  expect_warning(des_log <- el_prepare_inputs(Y_logical ~ X, df), "Coercing logical outcome")
  expect_true(is.list(des_log))
  expect_type(des_log$missingness_design[, "Y_logical"], "double")

  expect_warning(des_fac <- el_prepare_inputs(Y_factor ~ X, df), "Coercing two-level factor outcome")
  expect_true(is.list(des_fac))
  expect_true("Y_factor" %in% colnames(des_fac$missingness_design))
})

test_that("response intercept is retained even when formula uses +0", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    Z = rnorm(4)
  )
  design <- el_prepare_inputs(Y_miss ~ 1 | Z + 0, df)
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
  design <- el_prepare_inputs(Y_miss ~ . | ., df)
  expect_false("Y_miss" %in% colnames(design$aux_design_full))
  expect_false(anyNA(design$aux_design_full))
})

test_that("el_prepare_inputs forbids outcome in auxiliary constraints", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4),
    Z = rnorm(4)
  )
  expect_error(
    el_prepare_inputs(Y_miss ~ I(Y_miss^2) + X | Z, df),
    "outcome cannot appear",
    fixed = FALSE
  )
})

test_that("intercept-only auxiliaries are ignored silently", {
  df <- data.frame(
    Y_miss = c(1, NA, 3, NA),
    Z = rnorm(4)
  )
  design <- el_prepare_inputs(Y_miss ~ 1 | Z, df)
  expect_equal(ncol(design$aux_design_full), 0)
})

test_that("explicit auxiliary intercept is dropped silently", {
  df <- data.frame(
    Y_miss = c(1, NA, 3, NA),
    X = rnorm(4)
  )
  design <- el_prepare_inputs(Y_miss ~ 1 + X, df)
  expect_equal(colnames(design$aux_design_full), "X")
})

test_that("auxiliary '-1' or '+0' are overridden to k-1 (warning)", {
  set.seed(123)
  df <- data.frame(
    Y_miss = c(1, NA, 3, NA),
    X = rnorm(4),
    Z = rnorm(4)
  )
# '-1' should still drop intercept
  design_m1 <- el_prepare_inputs(Y_miss ~ X - 1 | Z, df)
  expect_false("(Intercept)" %in% colnames(design_m1$aux_design_full))
  expect_equal(colnames(design_m1$aux_design_full), "X")

# '+0' should behave identically
  design_p0 <- el_prepare_inputs(Y_miss ~ X + 0 | Z, df)
  expect_false("(Intercept)" %in% colnames(design_p0$aux_design_full))
  expect_equal(colnames(design_p0$aux_design_full), "X")
})

test_that("el_prepare_inputs rejects formulas with more than two RHS partitions", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X1 = rnorm(4),
    X2 = rnorm(4),
    X3 = rnorm(4)
  )
  expect_error(
    el_prepare_inputs(Formula::Formula(Y_miss ~ X1 | X2 | X3), df),
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
    el_prepare_inputs(Y_miss ~ Y_miss + X, df),
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
    el_prepare_inputs(Y_miss ~ I(0 * X + 1) + X, df),
    regexp = "Auxiliary covariate .* zero variance",
    fixed = FALSE
  )
})

test_that("el_prepare_inputs handles language objects coerced to formula", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    X = rnorm(4)
  )
  lang_formula <- as.call(list(as.name("~"), quote(Y_miss), quote(X)))
  parsed <- el_prepare_inputs(lang_formula, df)
  expect_identical(parsed$outcome_expr, "Y_miss")
  expect_true("(Intercept)" %in% colnames(parsed$missingness_design))
})

test_that("response-only formulas drop RHS2 entirely", {
  set.seed(55)
  df <- data.frame(
    Y_miss = c(rnorm(4), NA),
    X = rnorm(5)
  )
  parsed <- el_prepare_inputs(Y_miss ~ X, df)
  expect_identical(colnames(parsed$missingness_design), c("(Intercept)", "Y_miss"))
  expect_equal(colnames(parsed$aux_design_full), "X")
})

test_that("offset() usage is rejected on auxiliary and missingness partitions", {
  set.seed(123)
  df <- data.frame(
    Y_miss = c(rnorm(4), NA),
    X = rnorm(5),
    Z = rnorm(5)
  )
  expect_error(
    el_prepare_inputs(Y_miss ~ offset(X) + Z, df),
    "Offsets .* auxiliary predictors",
    fixed = FALSE
  )
  expect_error(
    el_prepare_inputs(Y_miss ~ Z | offset(X), df),
    "Offsets .* missingness predictors",
    fixed = FALSE
  )
})

test_that("auxiliary predictors with zero variance among respondents trigger an error", {
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA),
    A = c(5, 2, 5, 3)
  )
  expect_error(
    el_prepare_inputs(Y_miss ~ A, df),
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
    design <- el_prepare_inputs(Y_miss ~ X | Z, df),
    regexp = "Missingness-model predictor .* zero variance",
    fixed = FALSE
  )
  expect_true(is.list(design))
})

test_that("auxiliary factors use k-1 coding even with explicit -1/0", {
  df <- data.frame(
    Y_miss = c(1, 2, 3, NA),
    F = factor(c("a", "b", "c", "a"))
  )
  des <- el_prepare_inputs(Y_miss ~ 0 + F, df)
  expect_identical(colnames(des$aux_design_full), c("Fb", "Fc"))
  des2 <- el_prepare_inputs(Y_miss ~ F - 1, df)
  expect_identical(colnames(des2$aux_design_full), c("Fb", "Fc"))
})

test_that("el_prepare_inputs handles dot formulas with mixed numeric and factor auxiliaries", {
  set.seed(2025)
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA, 3, 4),
    X1 = rnorm(6),
    X2 = rnorm(6),
    F = factor(c("a", "b", "c", "a", "b", "c"))
  )
  design <- el_prepare_inputs(Y_miss ~ . | ., df)
  aux_cols <- colnames(design$aux_design_full)

# Outcome and intercept should not appear in auxiliaries
  expect_false("Y_miss" %in% aux_cols)
  expect_false("(Intercept)" %in% aux_cols)

# Factor should contribute k-1 columns
  fac_cols <- grep("^F", aux_cols, value = TRUE)
  expect_equal(length(fac_cols), nlevels(df$F) - 1)
})
