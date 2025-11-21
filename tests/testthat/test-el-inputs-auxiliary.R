test_that("el_resolve_auxiliaries works for data.frame with level drops", {
  set.seed(1)
  n <- 20
  full <- data.frame(
    Y = rnorm(n),
    f = factor(sample(c("A", "B"), n, replace = TRUE))
  )
  aux_formula <- ~ f - 1
  auxiliary_design_full <- model.matrix(aux_formula, data = full)
  respondent_mask <- full$f == "A"
  out <- el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = NULL)
  expect_true(is.matrix(out$auxiliary_design))
  expect_true(length(out$means) >= 1)
  expect_setequal(colnames(out$auxiliary_design), names(out$means))
  expect_equal(nrow(out$auxiliary_design), sum(respondent_mask))
})

test_that("el_resolve_auxiliaries computes design-weighted means for survey.design", {
  skip_if_not_installed("survey")
  set.seed(2)
  n <- 30
  full <- data.frame(
    Y = rnorm(n),
    f = factor(sample(c("A", "B"), n, replace = TRUE)),
    w = runif(n, 0.5, 2)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = full)
  aux_formula <- ~ f - 1
  auxiliary_design_full <- model.matrix(aux_formula, data = des$variables)
  respondent_mask <- des$variables$f == "A"
  out <- el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = NULL, weights_full = weights(des))
  mm_full <- model.matrix(aux_formula, data = full)
  mu_expected <- as.numeric(colSums(mm_full * full$w) / sum(full$w))
  names(mu_expected) <- colnames(mm_full)
  mu_expected <- mu_expected[colnames(out$auxiliary_design)]
  expect_equal(unname(out$means), unname(mu_expected), tolerance = 1e-12)
})

test_that("el_resolve_auxiliaries warns on extra names in auxiliary_means and ignores them", {
  set.seed(3)
  n <- 25
  full <- data.frame(
    Y = rnorm(n),
    X = rnorm(n)
  )
  auxiliary_design_full <- model.matrix(~ X - 1, data = full)
  respondent_mask <- rep(TRUE, n)
  aux_means_supplied <- c(X = 0, EXTRA = 123)
  out <- expect_warning(
    el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = aux_means_supplied),
    regexp = "Ignoring unused names in 'auxiliary_means'",
    fixed = FALSE
  )
  expect_equal(names(out$means), colnames(out$auxiliary_design))
  expect_false("EXTRA" %in% names(out$means))
})

test_that("el_resolve_auxiliaries rejects NA in full auxiliary data when means missing", {
  auxiliary_design_full <- rbind(
    c(X = 1),
    c(X = NA_real_),
    c(X = 0)
  )
  respondent_mask <- c(TRUE, FALSE, TRUE)
  expect_error(
    el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = NULL),
    regexp = "Auxiliary variables contain NA values",
    fixed = TRUE
  )
})

test_that("el_resolve_auxiliaries allows NA among nonrespondents when auxiliary_means supplied", {
  auxiliary_design_full <- rbind(
    c(X = 1),
    c(X = NA_real_),
    c(X = 0)
  )
  respondent_mask <- c(TRUE, FALSE, TRUE)
  out <- el_resolve_auxiliaries(
    auxiliary_design_full,
    respondent_mask = respondent_mask,
    auxiliary_means = c(X = 0.5)
  )
  expect_equal(nrow(out$auxiliary_design), sum(respondent_mask))
  expect_equal(colnames(out$auxiliary_design), "X")
  expect_equal(out$means[["X"]], 0.5)
})

test_that("auxiliary_means names match transformed columns (I(X^2))", {
  set.seed(123)
  n <- 60
  X <- rnorm(n)
  Y_true <- 1 + 0.2 * X + rnorm(n, sd = 0.3)
  R <- rbinom(n, 1, plogis(-0.3 + 0.5 * scale(Y_true)[, 1]))
  df <- data.frame(Y_miss = Y_true, X = X)
  df$Y_miss[R == 0] <- NA_real_
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

  aux_fml <- stats::as.formula("~ 0 + I(X^2)")
  auxiliary_design_full <- model.matrix(aux_fml, data = df)
  respondent_mask <- !is.na(df$Y_miss)
  expect_error(
    el_resolve_auxiliaries(
      auxiliary_design_full,
      respondent_mask = respondent_mask,
      auxiliary_means = c(X2 = 0)
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

test_that("EL errors early when auxiliaries contain NA (before solver)", {
  set.seed(123)
  n <- 50
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.5 * X + 0.3 * Z + rnorm(n)
  p <- plogis(-0.3 + 0.4 * scale(Y)[, 1] + 0.2 * Z)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_
  idx <- which(!is.na(df$Y_miss))[1]
  df$X[idx] <- NA_real_

  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
  expect_error(
    nmar(Y_miss ~ X, data = df, engine = eng),
    regexp = "contains NA values",
    info = "Auxiliary NA should be caught before solver setup"
  )
})

test_that("Infeasible/inconsistent auxiliaries trigger a warning and diagnostics", {
  set.seed(7002)
  N <- 800
  X <- rnorm(N)
  Y <- 1 + 0.5 * X + rnorm(N)
  p <- plogis(-0.4 + 0.3 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df$Y_miss[!R] <- NA_real_
  eng <- el_engine(auxiliary_means = c(X = 10), variance_method = "none", on_failure = "return")
  expect_warning(
    fit <- nmar(Y_miss ~ X, data = df, engine = eng),
    regexp = "Auxiliary means appear far from respondents' support"
  )
  expect_true("auxiliary_inconsistency_max_z" %in% names(fit$diagnostics))
  expect_true("auxiliary_inconsistency_cols" %in% names(fit$diagnostics))
})

test_that("EL removes auxiliary intercept (+1) silently", {
  set.seed(101)
  n <- 80
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.5 * X + rnorm(n)
  p <- plogis(-0.3 + 0.4 * scale(Y)[, 1])
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_

  eng <- el_engine(variance_method = "none")

  res1 <- nmar(Y_miss ~ X + 1, data = df, engine = eng, trace_level = 0)
  expect_s3_class(res1, "nmar_result_el")

  res2 <- nmar(Y_miss ~ X + 1 | Z, data = df, engine = eng, trace_level = 0)
  expect_s3_class(res2, "nmar_result_el")
})
