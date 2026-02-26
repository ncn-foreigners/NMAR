test_that("induced-logit runs end-to-end via nmar() on survey.design (variance_method='none')", {
  skip_if_not_installed("survey")

  set.seed(9101)
  n <- 400
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 2 + 0.6 * x1 - 0.2 * x2 + rnorm(n)

# Create nonresponse depending on y and x1; keep some missingness.
  p <- plogis(-0.5 + 0.25 * scale(y)[, 1] + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  w <- runif(n, 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  eng <- induced_logit_engine(variance_method = "none", on_failure = "error")
  res <- nmar(y_miss ~ x1 + x2 | x1, data = des, engine = eng)

  expect_s3_class(res, "nmar_result_induced_logit")
  expect_true(res$converged)
  expect_true(is.finite(res$y_hat))
  expect_true(is.na(res$se))
  expect_true(isTRUE(res$sample$is_survey))

  expect_true(is.finite(res$inference$df))
  expect_equal(res$inference$df, survey::degf(des))

  expect_true(is.finite(res$diagnostics$alpha0_hat_paper))
  expect_equal(res$diagnostics$survey_design_policy, "strict")
  expect_true(is.list(res$diagnostics$survey_assumptions))
  expect_true(is.finite(res$diagnostics$response_model_rank))
  expect_true(is.finite(res$diagnostics$response_model_ncol))
  expect_true(is.finite(res$diagnostics$response_model_condition_number))
  expect_null(res$diagnostics$engine)
  expect_null(res$extra$raw)

# Coefficients and fitted probabilities are exposed through shared S3
  b <- coef(res)
  expect_true(is.numeric(b) && length(b) >= 2)
  fv <- fitted(res)
  expect_true(is.numeric(fv))
  expect_equal(length(fv), nrow(df))
})

test_that("induced-logit bootstrap variance works end-to-end on survey.design", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  set.seed(9102)
  n <- 250
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.5 + 0.7 * x1 - 0.4 * x2 + rnorm(n)

  p <- plogis(-0.4 + 0.25 * scale(y)[, 1] + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  w <- runif(n, 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  eng <- induced_logit_engine(variance_method = "bootstrap", bootstrap_reps = 20, on_failure = "error")
  res <- nmar(y_miss ~ x1 + x2 | x1, data = des, engine = eng)

  expect_true(res$converged)
  expect_true(is.finite(res$y_hat))
  expect_true(is.finite(res$se))
  expect_true(res$se >= 0)
  expect_equal(res$inference$variance_method, "bootstrap")
})

test_that("induced-logit keep_fits=TRUE stores raw fit objects on survey.design", {
  skip_if_not_installed("survey")

  set.seed(9104)
  n <- 220
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 0.9 + 0.6 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.4 + 0.2 * scale(y)[, 1] + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  des <- survey::svydesign(ids = ~1, weights = ~ runif(n, 0.5, 2), data = df)
  eng <- induced_logit_engine(variance_method = "none", keep_fits = TRUE, on_failure = "error")
  res <- nmar(y_miss ~ x1 + x2 | x1, data = des, engine = eng)

  expect_true(res$converged)
  expect_true(is.list(res$extra$raw))
  expect_true(all(c("spec", "mu_fit", "induced_glm") %in% names(res$extra$raw)))
})
