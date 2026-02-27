test_that("induced-logit engine runs end-to-end via nmar() (IID, variance_method='none')", {
  set.seed(9091)
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 2 + 0.5 * x1 - 0.3 * x2 + rnorm(n)

  p <- plogis(-0.6 + 0.3 * scale(y)[, 1] + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  eng <- induced_logit_engine(variance_method = "none", on_failure = "error")
  res <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng)

  expect_s3_class(res, "nmar_result_induced_logit")
  expect_true(res$converged)
  expect_true(is.finite(res$y_hat))
  expect_true(is.na(res$se))

  expect_equal(res$meta$engine_name, "induced_logistic")
  expect_equal(res$inference$variance_method, "none")
  expect_equal(res$diagnostics$survey_design_policy, "strict")
  expect_true(is.finite(res$diagnostics$response_model_rank))
  expect_true(is.finite(res$diagnostics$response_model_ncol))
  expect_true(is.finite(res$diagnostics$response_model_condition_number))
  expect_null(res$diagnostics$engine)
  expect_null(res$extra$raw)

# Coefficients and fitted response probabilities are exposed through shared S3.
  b <- coef(res)
  expect_true(is.numeric(b) && length(b) >= 2)
  fv <- fitted(res)
  expect_true(is.numeric(fv))
  expect_equal(length(fv), nrow(df))

  expect_error(utils::capture.output(print(res)), NA)
  s <- summary(res)
  expect_s3_class(s, "summary_nmar_result_induced_logit")
  expect_error(utils::capture.output(print(s)), NA)
})

test_that("induced-logit on_failure='return' returns informative non-converged results", {
  set.seed(9092)
  n <- 30
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x1 + x2 + rnorm(n)
  r <- rbinom(n, 1, 0.7)
  if (all(r == 0L)) r[1] <- 1L
  if (all(r == 1L)) r[1] <- 0L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  res <- induced_logit.data.frame(
    data = df,
    formula = y_miss ~ .,
    variance_method = "none",
    on_failure = "return"
  )

  expect_s3_class(res, "nmar_result_induced_logit")
  expect_false(res$converged)
  expect_equal(res$estimate_name, "y_miss")
  expect_equal(res$sample$n_total, nrow(df))
  expect_equal(res$sample$n_respondents, sum(!is.na(df$y_miss)))
  expect_true(is.character(res$diagnostics$message) && nzchar(res$diagnostics$message))
})

test_that("induced-logit control can force nonconvergence and is surfaced", {
  set.seed(90921)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 2 + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  p <- plogis(-0.6 + 0.3 * scale(y)[, 1] + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

# maxit=1 should fail the missingness GLM convergence check.
  eng <- induced_logit_engine(
    variance_method = "none",
    on_failure = "return",
    control = list(maxit = 1)
  )

  res <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng)
  expect_s3_class(res, "nmar_result_induced_logit")
  expect_false(res$converged)
  expect_true(grepl("did not converge", res$diagnostics$message))
})

test_that("induced-logit invalid control list is surfaced as a user-facing failure", {
  set.seed(90922)
  n <- 180
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.5 + 0.5 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.1 * x1)
  r <- il_force_mixed_01(rbinom(n, 1, p))
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  eng <- induced_logit_engine(
    variance_method = "none",
    on_failure = "return",
    control = list(maxit = "not-an-integer")
  )
  res <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng)
  expect_false(res$converged)
  expect_match(res$diagnostics$message, "Invalid `control` for glm\\.control\\(\\)")
})

test_that("induced-logit bootstrap variance works end-to-end (IID)", {
  set.seed(9093)
  n <- 250
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)

  p <- plogis(-0.4 + 0.25 * y + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  eng <- induced_logit_engine(variance_method = "bootstrap", bootstrap_reps = 30, on_failure = "error")
  res <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng)

  expect_true(res$converged)
  expect_true(is.finite(res$y_hat))
  expect_true(is.finite(res$se))
  expect_true(res$se >= 0)
  expect_equal(res$inference$variance_method, "bootstrap")
})

test_that("induced-logit standardize=TRUE is a pure reparameterization (IID)", {
  set.seed(9094)
  n <- 600
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.7 + 0.9 * x1 - 0.6 * x2 + rnorm(n)

# Mildly imbalanced missingness to exercise the induced GLM.
  p <- plogis(-0.8 + 0.25 * y + 0.2 * x1 - 0.1 * x2)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  eng0 <- induced_logit_engine(variance_method = "none", standardize = FALSE, on_failure = "error")
  eng1 <- induced_logit_engine(variance_method = "none", standardize = TRUE, on_failure = "error")

# Keep x1 smaller than the mu-model so mu_hat is not collinear with x1
  res0 <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng0)
  res1 <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng1)

  expect_true(res0$converged)
  expect_true(res1$converged)

# Primary estimate and paper gamma should be invariant up to tolerance.
  expect_equal(res1$y_hat, res0$y_hat, tolerance = 1e-10)
  expect_equal(res1$diagnostics$gamma_hat_paper, res0$diagnostics$gamma_hat_paper, tolerance = 1e-10)

# Fitted response probabilities are invariant
  expect_equal(fitted(res1), fitted(res0), tolerance = 1e-10)

# Missingness-model coefficients/vcov are reported on the original scale
  expect_equal(coef(res1), coef(res0), tolerance = 1e-8)
  if (is.matrix(res0$model$vcov) && is.matrix(res1$model$vcov)) {
    expect_equal(res1$model$vcov, res0$model$vcov, tolerance = 1e-6)
  }

  expect_true(is.null(res0$extra$scaling))
  expect_true(is.list(res1$extra$scaling))
  expect_null(res0$extra$raw)
  expect_null(res1$extra$raw)
})

test_that("induced-logit keep_fits=TRUE stores raw fit objects", {
  set.seed(9095)
  n <- 220
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 0.8 + 0.5 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.25 * y + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  eng <- induced_logit_engine(variance_method = "none", keep_fits = TRUE, on_failure = "error")
  res <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng)

  expect_true(res$converged)
  expect_true(is.list(res$extra$raw))
  expect_true(all(c("spec", "mu_fit", "induced_glm") %in% names(res$extra$raw)))
})

test_that("induced-logit on_failure='error' propagates computation errors", {
  set.seed(9096)
  n <- 140
  x1 <- rnorm(n)
  x_dup <- x1
  x2 <- rnorm(n)
  y <- 1 + 0.6 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x_dup = x_dup, x2 = x2)

  eng <- induced_logit_engine(variance_method = "none", on_failure = "error")
  expect_error(
    nmar(y_miss ~ x1 + x2 | x1 + x_dup, data = df, engine = eng),
    "rank deficient"
  )
})
