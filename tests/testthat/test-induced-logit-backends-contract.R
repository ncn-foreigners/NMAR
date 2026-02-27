expect_il_backend_fit_contract <- function(fit, n, standardized) {
  expect_true(is.list(fit))

  required_names <- c(
    "point", "models", "warnings", "diagnostics", "sample"
  )
  expect_true(all(required_names %in% names(fit)))

  scalar_num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  expect_true(is.list(fit$point))
  expect_true(scalar_num(fit$point$tau_hat))
  expect_true(scalar_num(fit$point$eta_hat))
  expect_true(scalar_num(fit$point$mu_bar))
  expect_true(scalar_num(fit$point$gamma_hat_paper))
  expect_true(scalar_num(fit$point$m2_over_m1))
  expect_true(scalar_num(fit$point$log_m1_hat))
  expect_true(is.numeric(fit$point$alpha_hat_paper) && length(fit$point$alpha_hat_paper) == 1L)
  expect_true(is.numeric(fit$point$alpha0_hat_paper) && length(fit$point$alpha0_hat_paper) == 1L)

  expect_true(is.list(fit$models))
  expect_false(is.null(fit$models$mu_fit))
  expect_false(is.null(fit$models$induced_glm))

  expect_true(is.numeric(fit$models$induced_glm_coef))
  expect_true(!is.null(names(fit$models$induced_glm_coef)))
  expect_true(length(fit$models$induced_glm_coef) >= 2L)
  expect_true(all(is.finite(fit$models$induced_glm_coef)))

  expect_true(is.matrix(fit$models$induced_glm_vcov))
  expect_equal(nrow(fit$models$induced_glm_vcov), length(fit$models$induced_glm_coef))
  expect_equal(ncol(fit$models$induced_glm_vcov), length(fit$models$induced_glm_coef))
  expect_equal(rownames(fit$models$induced_glm_vcov), names(fit$models$induced_glm_coef))
  expect_equal(colnames(fit$models$induced_glm_vcov), names(fit$models$induced_glm_coef))

  expect_true(is.numeric(fit$models$fitted_values))
  expect_equal(length(fit$models$fitted_values), n)

  expect_true(is.list(fit$warnings))
  expect_true(all(c("mu", "induced_glm") %in% names(fit$warnings)))
  expect_true(is.character(fit$warnings$mu))
  expect_true(is.character(fit$warnings$induced_glm))
  expect_true(is.list(fit$diagnostics))
  rm <- fit$diagnostics$response_model %||% list()
  expect_true(scalar_num(rm$rank))
  expect_true(scalar_num(rm$ncol))
  expect_true(scalar_num(rm$condition_number))
  expect_true(rm$rank <= rm$ncol)

  expect_true(is.list(fit$sample))
  expect_true(is.numeric(fit$sample$n_total) && length(fit$sample$n_total) == 1L && fit$sample$n_total == n)
  expect_true(is.numeric(fit$sample$n_respondents) && length(fit$sample$n_respondents) == 1L)
  expect_true(fit$sample$n_respondents >= 1L && fit$sample$n_respondents <= n)

  if (isTRUE(standardized)) {
    expect_true(is.list(fit$diagnostics$scaling))
    expect_true(all(c("mu", "x1") %in% names(fit$diagnostics$scaling)))
  } else {
    expect_true(is.null(fit$diagnostics$scaling))
  }
}

test_that("induced-logit backends return a stable fit contract (IID)", {
  set.seed(9301)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.3 + 0.7 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.6 + 0.25 * y + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  f <- y_miss ~ x1 + x2 | x1
  spec <- induced_logit_prepare_inputs(formula = f, data = df)
  backend <- il_backend_iid(df)

  fit0 <- il_fit_from_backend(spec, backend = backend, standardize = FALSE)
  expect_il_backend_fit_contract(fit0, n = n, standardized = FALSE)

  fit1 <- il_fit_from_backend(spec, backend = backend, standardize = TRUE)
  expect_il_backend_fit_contract(fit1, n = n, standardized = TRUE)
})

test_that("induced-logit backends return a stable fit contract (survey)", {
  skip_if_not_installed("survey")

  set.seed(9302)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.3 + 0.7 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.6 + 0.25 * scale(y)[, 1] + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  w <- runif(n, 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  f <- y_miss ~ x1 + x2 | x1
  spec <- induced_logit_prepare_inputs(formula = f, data = des$variables)
  backend <- il_backend_survey(des)

  fit0 <- il_fit_from_backend(spec, backend = backend, standardize = FALSE)
  expect_il_backend_fit_contract(fit0, n = n, standardized = FALSE)

  fit1 <- il_fit_from_backend(spec, backend = backend, standardize = TRUE)
  expect_il_backend_fit_contract(fit1, n = n, standardized = TRUE)
})

test_that("induced-logit identifiability diagnostics enforce rank/conditioning thresholds", {
  expect_error(
    il_enforce_response_model_identifiability(list(rank = 2, ncol = 3, condition_number = 10)),
    "rank deficient"
  )

  expect_error(
    il_enforce_response_model_identifiability(list(rank = 3, ncol = 3, condition_number = 1e13)),
    "nearly non-identifiable"
  )

  weak <- il_enforce_response_model_identifiability(list(rank = 3, ncol = 3, condition_number = 1e9))
  expect_true(is.list(weak))
  expect_true(is.character(weak$warning) && nzchar(weak$warning))
  expect_equal(weak$diag$rank, 3)
  expect_equal(weak$diag$ncol, 3)
})

test_that("induced-logit prepare_inputs caches model matrices consistently", {
  set.seed(9303)
  n <- 120
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.4 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.1 * x1)
  r <- il_force_mixed_01(rbinom(n, 1, p))
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  spec <- induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df)
  expect_true(is.matrix(spec$mu_mat_full))
  expect_true(is.matrix(spec$x1_mat_full))
  expect_equal(nrow(spec$mu_mat_full), n)
  expect_equal(nrow(spec$x1_mat_full), n)

  mm_mu <- stats::model.matrix(spec$mu_rhs, data = spec$model_frame)
  mm_x1 <- stats::model.matrix(spec$x1_rhs, data = spec$model_frame)
  expect_equal(colnames(spec$mu_mat_full), colnames(mm_mu))
  expect_equal(colnames(spec$x1_mat_full), colnames(mm_x1))
  expect_equal(spec$mu_mat_full, mm_mu, tolerance = 0)
  expect_equal(spec$x1_mat_full, mm_x1, tolerance = 0)
})
