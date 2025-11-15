test_that("auxiliary_means names match transformed columns (I(X^2))", {
  set.seed(123)
  n <- 60
  X <- rnorm(n)
  Y_true <- 1 + 0.2 * X + rnorm(n, sd = 0.3)
  R <- rbinom(n, 1, plogis(-0.3 + 0.5 * scale(Y_true)[, 1]))
  df <- data.frame(Y_miss = Y_true, X = X)
  df$Y_miss[R == 0] <- NA_real_

# Correctly named auxiliary means for transformed column
  eng <- el_engine(auxiliary_means = c(`I(X^2)` = 0), variance_method = "none")
  expect_silent({
    fit <- nmar(Y_miss ~ I(X^2), data = df, engine = eng)
  })
  expect_s3_class(fit, "nmar_result_el")
})

test_that("misnamed auxiliary_means trigger an error", {
  set.seed(124)
  n <- 50
  X <- rnorm(n)
  Y_true <- 0.5 + 0.4 * X + rnorm(n, sd = 0.4)
  R <- rbinom(n, 1, plogis(-0.2 + 0.5 * scale(Y_true)[, 1]))
  df <- data.frame(Y_miss = Y_true, X = X)
  df$Y_miss[R == 0] <- NA_real_

# Build aux-only formula for EL internal utilities (no input pipeline)
  aux_fml <- stats::as.formula("~ 0 + I(X^2)")
  resp_df <- df[!is.na(df$Y_miss), , drop = FALSE]

# Misnamed means (no match) -> informative error
  auxiliary_design_full <- model.matrix(aux_fml, data = df)
  respondent_mask <- !is.na(df$Y_miss)
  expect_error(
    NMAR:::el_resolve_auxiliaries(
      auxiliary_design_full,
      respondent_mask = respondent_mask,
      auxiliary_means = c(X2 = 0) # misnamed
    ),
    "auxiliary_means must supply entries",
    fixed = FALSE
  )
})

test_that("partially mismatched auxiliary_means warn and use only matched names", {
  set.seed(42)
  n <- 120
  X <- rnorm(n)
  Y <- 1 + 0.5 * X + rnorm(n, sd = 0.3)
  R <- rbinom(n, 1, plogis(-0.3 + 0.6 * scale(Y)[, 1]))
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[R == 0] <- NA_real_
# Provide one matching and one unmatched name; expect a warning about dropping unmatched
  eng <- el_engine(auxiliary_means = c(X = 0, foo = 10), variance_method = "none")
  expect_warning(
    fit <- nmar(Y_miss ~ X, data = df, engine = eng),
    "Ignoring unused names in 'auxiliary_means'",
    fixed = FALSE
  )
  expect_s3_class(fit, "nmar_result_el")
  expect_true(isTRUE(fit$converged))
})

test_that("nmar errors when auxiliary_means omit required columns", {
  set.seed(101)
  n <- 80
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- 1 + 0.5 * X1 - 0.2 * X2 + rnorm(n)
  R <- rbinom(n, 1, plogis(-0.3 + 0.4 * scale(Y)[, 1]))
  df <- data.frame(Y_miss = Y, X1 = X1, X2 = X2)
  df$Y_miss[R == 0] <- NA_real_
  eng <- el_engine(auxiliary_means = c(X1 = 0), variance_method = "none")
  expect_error(
    nmar(Y_miss ~ X1 + X2, data = df, engine = eng),
    "auxiliary_means must supply entries",
    fixed = FALSE
  )
})
