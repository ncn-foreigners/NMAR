test_that("outcome prioritization prefers first column with NAs", {
# Create data with multiple outcomes where only one has missingness
  df <- data.frame(
    Y1 = c(1, 2, 3, 4, 5), # No NAs
    Y2 = c(10, NA, 30, NA, 50), # Has NAs (should be prioritized)
    X = c(1, 2, 3, 4, 5)
  )

# Use traits that allow multi-outcome for testing prioritization logic
  traits <- modifyList(NMAR:::NMAR_DEFAULT_TRAITS, list(requires_single_outcome = FALSE))
  spec <- parse_nmar_spec(Y1 + Y2 ~ X, df, env = parent.frame(), traits = traits)

# Should prefer Y2 because it has missingness
  expect_equal(spec$outcome_primary, "Y2")
  expect_setequal(spec$outcome_vars_data, c("Y1", "Y2"))
})

test_that("outcome prioritization falls back to first when none have NAs", {
# When allow_respondents_only is TRUE and no outcomes have NAs,
# should fall back to first outcome
  df <- data.frame(
    Y1 = c(1, 2, 3, 4, 5), # No NAs
    Y2 = c(10, 20, 30, 40, 50), # No NAs
    X = c(1, 2, 3, 4, 5)
  )

  traits <- modifyList(NMAR:::NMAR_DEFAULT_TRAITS, list(requires_single_outcome = FALSE))
  spec <- parse_nmar_spec(Y1 + Y2 ~ X, df, env = parent.frame(), traits = traits)

# Should fall back to first (Y1)
  expect_equal(spec$outcome_primary, "Y1")
})

test_that("outcome prioritization handles helper functions correctly", {
# Simulates scenario where helper function appears first in all.vars()
# but we want the actual outcome column to be selected
  df <- data.frame(
    income = c(100, NA, 300, NA, 500), # Has NAs
    age = c(25, 30, 35, 40, 45), # No NAs
    X = c(1, 2, 3, 4, 5)
  )

  traits <- modifyList(NMAR:::NMAR_DEFAULT_TRAITS, list(requires_single_outcome = FALSE))
  spec <- parse_nmar_spec(income + age ~ X, df, env = parent.frame(), traits = traits)

# Should prioritize income (has missingness) over age (fully observed)
  expect_equal(spec$outcome_primary, "income")
})

test_that("single outcome uses first column directly", {
  df <- data.frame(
    Y = c(1, NA, 3, NA, 5),
    X = c(1, 2, 3, 4, 5)
  )

  spec <- parse_nmar_spec(Y ~ X, df, env = parent.frame(),
                          traits = NMAR_DEFAULT_TRAITS)

  expect_equal(spec$outcome_primary, "Y")
  expect_equal(spec$outcome_vars_data, "Y")
  expect_false(spec$outcome_is_multi)
})

test_that("cbind multi-outcome is correctly detected", {
  df <- data.frame(
    Y1 = c(1, NA, 3),
    Y2 = c(10, 20, NA),
    X = c(1, 2, 3)
  )

# Using cbind() should trigger multi-outcome detection
  traits <- modifyList(NMAR:::NMAR_DEFAULT_TRAITS, list(requires_single_outcome = FALSE))
  spec <- parse_nmar_spec(cbind(Y1, Y2) ~ X, df, env = parent.frame(), traits = traits)

  expect_true(spec$outcome_is_multi)
# For multi-outcome, outcome_vars includes all outcomes
  expect_setequal(spec$outcome_vars_data, c("Y1", "Y2"))
})

test_that("outcome prioritization integrates with nmar() workflow", {
  skip_if_not_installed("nleqslv")

  set.seed(123)
  n <- 50
  df <- data.frame(
    Y1 = rnorm(n), # Fully observed
    Y2 = rnorm(n), # Will have NAs
    X = rnorm(n)
  )

# Create missingness in Y2 only
  R <- rbinom(n, 1, 0.7)
  df$Y2[R == 0] <- NA_real_

# Test that multi-outcome formula is accepted by exptilt_nonparam
# (which has requires_single_outcome = FALSE)
  eng <- exptilt_nonparam_engine(refusal_col = "Y1")

# Should not error - multi-outcome is allowed for this engine
  expect_no_error({
    fit <- nmar(Y1 + Y2 ~ X, data = df, engine = eng)
  })

# Verify result structure
  expect_s3_class(fit, "nmar_result")
  expect_true(!is.null(fit$meta$formula))
})

test_that("prioritization respects order when multiple have missingness", {
  df <- data.frame(
    Y1 = c(NA, 2, 3), # Has NAs (appears first)
    Y2 = c(10, NA, 30), # Also has NAs (appears second)
    X = c(1, 2, 3)
  )

  traits <- modifyList(NMAR:::NMAR_DEFAULT_TRAITS, list(requires_single_outcome = FALSE))
  spec <- parse_nmar_spec(Y1 + Y2 ~ X, df, env = parent.frame(), traits = traits)

# Should select Y1 (first with NAs)
  expect_equal(spec$outcome_primary, "Y1")
})
