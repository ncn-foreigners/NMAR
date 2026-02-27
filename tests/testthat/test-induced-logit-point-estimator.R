test_that("induced-logit point estimator runs and is numerically stable (IID)", {
  set.seed(24601)
  n <- 800

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eps <- rnorm(n)

  b0 <- 1.2
  b1 <- 0.8
  b2 <- -0.5
  mu <- b0 + b1 * x1 + b2 * x2
  y <- mu + eps

  alpha0 <- -1.3
  beta1 <- 0.4
  gamma <- 0.7

  p_r <- 1 / (1 + exp(alpha0 + beta1 * x1 + gamma * y))
  r <- rbinom(n, 1, p_r)
  if (all(r == 1L)) r[sample.int(n, 1)] <- 0L
  if (all(r == 0L)) r[sample.int(n, 1)] <- 1L

  df <- data.frame(
    y_miss = ifelse(r == 1L, y, NA_real_),
    x1 = x1,
    x2 = x2
  )

  spec <- induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df)
  fit <- il_fit_from_backend(spec = spec, backend = il_backend_iid(df), standardize = FALSE)

  expect_true(is.finite(fit$point$tau_hat))
  expect_true(is.finite(fit$point$gamma_hat_paper))
  expect_true(is.finite(fit$point$alpha0_hat_paper))
  expect_true(fit$sample$n_respondents > 0)
  expect_true(fit$sample$n_total == nrow(df))

# True E(Y) under the DGP is b0 (E[x1]=E[x2]=0, E[eps]=0).
  expect_lt(abs(fit$point$tau_hat - b0), 0.15)

# Stress-test the stabilized ratio helper directly (should not overflow).
  expect_equal(il_m2_over_m1_ratio(c(1, 2, 3), gamma = 1000), 3, tolerance = 1e-12)
  expect_equal(il_m2_over_m1_ratio(c(1, 2, 3), gamma = -1000), 1, tolerance = 1e-12)

# log M1 helper should be stable and match direct computation on moderate inputs
  eps <- c(-1, 0, 2)
  g <- 0.5
  expect_equal(il_log_m1_hat(eps, gamma = g), log(mean(exp(g * eps))), tolerance = 1e-12)

# Weighted helpers should match direct weighted expressions
  w <- c(2, 1, 3)
  num <- sum(w * eps * exp(g * eps))
  den <- sum(w * exp(g * eps))
  expect_equal(il_m2_over_m1_ratio_weighted(eps, gamma = g, weights = w), num / den, tolerance = 1e-12)
  expect_equal(il_log_m1_hat_weighted(eps, gamma = g, weights = w), log(sum(w * exp(g * eps)) / sum(w)), tolerance = 1e-12)

# Equal weights should reduce to the unweighted helpers
  w1 <- rep(1, length(eps))
  expect_equal(il_m2_over_m1_ratio_weighted(eps, gamma = g, weights = w1), il_m2_over_m1_ratio(eps, gamma = g), tolerance = 1e-12)
  expect_equal(il_log_m1_hat_weighted(eps, gamma = g, weights = w1), il_log_m1_hat(eps, gamma = g), tolerance = 1e-12)

# Zero weights should be safe: they must not affect the log-sum-exp shift
  eps2 <- c(1000, -1000)
  w2 <- c(0, 1)
  expect_equal(il_m2_over_m1_ratio_weighted(eps2, gamma = 1, weights = w2), -1000, tolerance = 1e-12)
  expect_equal(il_log_m1_hat_weighted(eps2, gamma = 1, weights = w2), -1000, tolerance = 1e-12)
})

test_that("induced-logit point estimator enforces fully observed covariates", {
  set.seed(24602)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x1 + x2 + rnorm(n)
  r <- rep(1L, n)
  r[sample.int(n, n %/% 3)] <- 0L

  df <- data.frame(
    y_miss = ifelse(r == 1L, y, NA_real_),
    x1 = x1,
    x2 = x2
  )

  df_bad <- df
  df_bad$x2[3] <- NA_real_

  expect_error(induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df_bad), "fully observed")

  df_bad2 <- df
  df_bad2$x1[7] <- NA_real_

  expect_error(induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df_bad2), "fully observed")
})

test_that("induced-logit IID estimate equals direct Eq. (13) reconstruction", {
  set.seed(24603)
  n <- 350
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.1 + 0.7 * x1 - 0.3 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.15 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  eng <- induced_logit_engine(variance_method = "none", on_failure = "error", keep_fits = TRUE)
  res <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng)
  expect_true(res$converged)

  spec <- res$extra$raw$spec
  mu_fit <- res$extra$raw$mu_fit
  glm_fit <- res$extra$raw$induced_glm

  mu_hat <- as.numeric(stats::predict(mu_fit, newdata = df))
  eps_hat <- as.numeric(df$y_miss[spec$respondent_mask] - mu_hat[spec$respondent_mask])
  gamma_hat <- as.numeric(-stats::coef(glm_fit)[["..nmar_mu_hat.."]])
  eta_hat <- mean(as.integer(spec$respondent_mask))
  m2_over_m1 <- sum(eps_hat * exp(gamma_hat * eps_hat)) / sum(exp(gamma_hat * eps_hat))
  tau_manual <- mean(mu_hat) + (1 - eta_hat) * m2_over_m1

  expect_equal(res$y_hat, tau_manual, tolerance = 1e-10)
  expect_true(is.finite(res$diagnostics$response_model_rank))
  expect_true(is.finite(res$diagnostics$response_model_ncol))
  expect_true(is.finite(res$diagnostics$response_model_condition_number))
  expect_equal(res$diagnostics$response_model_rank, res$diagnostics$response_model_ncol)
})

test_that("induced-logit IID response model fails on rank deficiency", {
  set.seed(24604)
  n <- 220
  x1 <- rnorm(n)
  x_dup <- x1
  x2 <- rnorm(n)
  y <- 1.2 + 0.7 * x1 - 0.4 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(
    y_miss = ifelse(r == 1L, y, NA_real_),
    x1 = x1,
    x_dup = x_dup,
    x2 = x2
  )
  f <- y_miss ~ x1 + x2 | x1 + x_dup

  expect_error(
    nmar(f, data = df, engine = induced_logit_engine(variance_method = "none", on_failure = "error")),
    "rank deficient"
  )

  res_fail <- induced_logit.data.frame(
    data = df,
    formula = f,
    variance_method = "none",
    on_failure = "return"
  )
  expect_false(res_fail$converged)
  expect_match(res_fail$diagnostics$message, "rank deficient")
})

test_that("induced-logit IID glm.fit path matches formula-based glm reference", {
  set.seed(24605)
  n <- 260
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 0.9 + 0.6 * x1 - 0.3 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.15 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  spec <- induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df)
  mu_res <- il_fit_mu_iid_core(spec = spec, standardize = FALSE)
  mu_hat <- as.numeric(stats::predict(mu_res$mu_fit, newdata = df))

  fit_fast <- il_fit_resp_iid_core(
    spec = spec,
    mu_hat = mu_hat,
    standardize = FALSE,
    glm_control = stats::glm.control()
  )

  d_ref <- df
  d_ref$..nmar_r.. <- as.integer(spec$respondent_mask)
  d_ref$..nmar_mu_hat.. <- mu_hat
  x1_terms <- spec$x1_term_labels
  if (is.null(x1_terms)) x1_terms <- character(0)
  x1_intercept <- spec$x1_intercept
  if (is.null(x1_intercept)) x1_intercept <- TRUE
  f_ref <- stats::reformulate(
    termlabels = c(x1_terms, "..nmar_mu_hat.."),
    response = "..nmar_r..",
    intercept = isTRUE(x1_intercept)
  )
  fit_ref <- stats::glm(f_ref, family = stats::binomial(), data = d_ref)

  expect_equal(
    unname(fit_fast$beta_glm[names(stats::coef(fit_ref))]),
    unname(stats::coef(fit_ref)),
    tolerance = 1e-10
  )
  expect_equal(
    fit_fast$fitted,
    as.numeric(stats::fitted(fit_ref)),
    tolerance = 1e-10
  )
  expect_equal(
    fit_fast$gamma_hat_paper,
    as.numeric(-stats::coef(fit_ref)[["..nmar_mu_hat.."]]),
    tolerance = 1e-12
  )
})

test_that("induced-logit captures separation-style warnings from glm.fit", {
  set.seed(24606)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)

# Outcome is unrelated to x1, response is near-separated by x1, which should
# trigger the classic "fitted probabilities numerically 0 or 1 occurred" warning
# while still converging
  y <- 0.2 * x2 + rnorm(n)
  p <- plogis(20 * x1)
  r <- il_force_mixed_01(rbinom(n, 1, p))

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  spec <- induced_logit_prepare_inputs(y_miss ~ x2 | x1, data = df)
  mu_res <- il_fit_mu_iid_core(spec = spec, standardize = FALSE)
  mu_hat <- as.numeric(stats::predict(mu_res$mu_fit, newdata = df))

  fit_fast <- il_fit_resp_iid_core(
    spec = spec,
    mu_hat = mu_hat,
    standardize = FALSE,
    glm_control = stats::glm.control(maxit = 100)
  )

  expect_true(is.character(fit_fast$warnings))
  expect_true(length(fit_fast$warnings) > 0L)
  expect_true(any(grepl("fitted probabilities numerically 0 or 1", fit_fast$warnings)))
})
