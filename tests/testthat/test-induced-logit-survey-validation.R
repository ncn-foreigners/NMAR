test_that("induced-logit survey.design matches IID when ids=~1 and weights=1", {
  skip_if_not_installed("survey")

  set.seed(9201)
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.7 * x1 - 0.4 * x2 + rnorm(n)

  p <- plogis(-0.6 + 0.15 * scale(y)[, 1] + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  des <- survey::svydesign(ids = ~1, weights = ~1, data = df)

  for (standardize in c(FALSE, TRUE)) {
    eng <- induced_logit_engine(variance_method = "none", standardize = standardize, on_failure = "error")
    res_df <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng)
    res_des <- nmar(y_miss ~ x1 + x2 | x1, data = des, engine = eng)

    expect_true(res_df$converged)
    expect_true(res_des$converged)

    expect_equal(res_des$y_hat, res_df$y_hat, tolerance = 1e-8)
    expect_equal(res_des$diagnostics$gamma_hat_paper, res_df$diagnostics$gamma_hat_paper, tolerance = 1e-8)
    expect_equal(fitted(res_des), fitted(res_df), tolerance = 1e-8)
    expect_equal(coef(res_des), coef(res_df), tolerance = 1e-6)
  }
})

test_that("induced-logit survey estimate is invariant to global weight scaling", {
  skip_if_not_installed("survey")

  set.seed(92015)
  n <- 420
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.3 + 0.6 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.1 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  w <- runif(n, 0.5, 2)
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2, w = w, w2 = 10 * w)
  des1 <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  des2 <- survey::svydesign(ids = ~1, weights = ~w2, data = df)
  eng <- induced_logit_engine(variance_method = "none", on_failure = "error", keep_fits = TRUE)

  res1 <- nmar(y_miss ~ x1 + x2 | x1, data = des1, engine = eng)
  res2 <- nmar(y_miss ~ x1 + x2 | x1, data = des2, engine = eng)

  expect_true(res1$converged)
  expect_true(res2$converged)
  expect_equal(res2$y_hat, res1$y_hat, tolerance = 1e-10)
  expect_equal(res2$diagnostics$gamma_hat_paper, res1$diagnostics$gamma_hat_paper, tolerance = 1e-10)
  expect_equal(coef(res2), coef(res1), tolerance = 1e-8)
  expect_equal(fitted(res2), fitted(res1), tolerance = 1e-8)
})

test_that("induced-logit survey standardize preserves no-intercept mu model", {
  skip_if_not_installed("survey")

  set.seed(92016)
  n <- 320
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 0.9 + 0.4 * x1 - 0.3 * x2 + rnorm(n)
  p <- plogis(-0.6 + 0.2 * y + 0.1 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  w <- runif(n, 0.5, 2)
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2, w = w)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  f <- y_miss ~ x1 + x2 - 1 | x1
  eng0 <- induced_logit_engine(variance_method = "none", standardize = FALSE, on_failure = "error")
  eng1 <- induced_logit_engine(variance_method = "none", standardize = TRUE, on_failure = "error")

  res0 <- nmar(f, data = des, engine = eng0)
  res1 <- nmar(f, data = des, engine = eng1)

  expect_true(res0$converged)
  expect_true(res1$converged)
  expect_equal(res1$y_hat, res0$y_hat, tolerance = 1e-8)
  expect_equal(res1$diagnostics$gamma_hat_paper, res0$diagnostics$gamma_hat_paper, tolerance = 1e-8)
})

il_manual_m2_over_m1_weighted <- function(eps_hat, gamma, weights) {
  eps_hat <- as.numeric(eps_hat)
  weights <- as.numeric(weights)
  if (length(eps_hat) != length(weights)) stop("length mismatch")
  if (length(eps_hat) == 0L) return(NA_real_)
  if (any(!is.finite(eps_hat)) || any(!is.finite(weights)) || any(weights < 0)) return(NA_real_)
  if (!is.finite(gamma)) return(NA_real_)
  lw <- log(weights) + gamma * eps_hat
  a <- max(lw, na.rm = TRUE)
  w_tilt <- exp(lw - a)
  den <- sum(w_tilt)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num <- sum(eps_hat * w_tilt)
  as.numeric(num / den)
}

il_manual_induced_logit_survey <- function(design, formula) {
  stopifnot(inherits(design, "survey.design"))
  vars <- design$variables

  spec <- induced_logit_prepare_inputs(formula = formula, data = vars)
  respondent_mask <- spec$respondent_mask
  r_vec <- as.integer(respondent_mask)

  w_full <- as.numeric(stats::weights(design))
  sum_w <- sum(w_full)

  mu_terms <- spec$mu_term_labels
  if (is.null(mu_terms)) mu_terms <- character(0)
  mu_intercept <- spec$mu_intercept
  if (is.null(mu_intercept)) mu_intercept <- TRUE
  mu_formula <- stats::reformulate(
    termlabels = mu_terms,
    response = spec$outcome,
    intercept = isTRUE(mu_intercept)
  )
  design_resp <- subset(design, respondent_mask)
  mu_fit <- survey::svyglm(mu_formula, design = design_resp)

  mu_mat_full <- spec$mu_mat_full
  mu_coef <- stats::coef(mu_fit)
  mu_hat <- as.numeric(mu_mat_full %*% mu_coef[colnames(mu_mat_full)])

  design_aug <- stats::update(design,
    ..nmar_r.. = r_vec,
    ..nmar_mu_hat.. = mu_hat
  )

  x1_terms <- spec$x1_term_labels
  if (is.null(x1_terms)) x1_terms <- character(0)
  x1_intercept <- spec$x1_intercept
  if (is.null(x1_intercept)) x1_intercept <- TRUE
  resp_formula <- stats::reformulate(
    termlabels = c(x1_terms, "..nmar_mu_hat.."),
    response = "..nmar_r..",
    intercept = isTRUE(x1_intercept)
  )
  resp_fit <- survey::svyglm(resp_formula, design = design_aug, family = stats::quasibinomial())
  beta <- stats::coef(resp_fit)
  gamma_hat_paper <- as.numeric(-beta[["..nmar_mu_hat.."]])

  eta_hat <- sum(w_full * r_vec) / sum_w
  mu_bar <- sum(w_full * mu_hat) / sum_w

  y_vec <- vars[[spec$outcome]]
  eps_hat <- as.numeric(y_vec[respondent_mask] - mu_hat[respondent_mask])
  w_resp <- w_full[respondent_mask]
  m2_over_m1 <- il_manual_m2_over_m1_weighted(eps_hat = eps_hat, gamma = gamma_hat_paper, weights = w_resp)

  list(
    tau_hat = as.numeric(mu_bar + (1 - eta_hat) * m2_over_m1),
    mu_hat = mu_hat,
    gamma_hat_paper = gamma_hat_paper,
    coef = beta,
    fitted = as.numeric(stats::fitted(resp_fit))
  )
}

test_that("induced-logit survey.design matches a manual reference calculation", {
  skip_if_not_installed("survey")

  set.seed(9202)
  n <- 450
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 2 + 0.6 * x1 - 0.2 * x2 + rnorm(n)

  p <- plogis(-0.4 + 0.2 * scale(y)[, 1] + 0.15 * x1 - 0.05 * x2)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  w <- runif(n, 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  f <- y_miss ~ x1 + x2 | x1

  ref <- il_manual_induced_logit_survey(design = des, formula = f)

  eng0 <- induced_logit_engine(variance_method = "none", standardize = FALSE, on_failure = "error")
  eng1 <- induced_logit_engine(variance_method = "none", standardize = TRUE, on_failure = "error")

  res0 <- nmar(f, data = des, engine = eng0)
  res1 <- nmar(f, data = des, engine = eng1)

  expect_equal(res0$y_hat, ref$tau_hat, tolerance = 1e-8)
  expect_equal(res0$diagnostics$gamma_hat_paper, ref$gamma_hat_paper, tolerance = 1e-8)
  expect_equal(fitted(res0), ref$fitted, tolerance = 1e-8)
  expect_equal(coef(res0), ref$coef, tolerance = 1e-6)

# standardize=TRUE should not change point estimate nor fitted probabilities.
  expect_equal(res1$y_hat, ref$tau_hat, tolerance = 1e-8)
  expect_equal(res1$diagnostics$gamma_hat_paper, ref$gamma_hat_paper, tolerance = 1e-8)
  expect_equal(fitted(res1), ref$fitted, tolerance = 1e-8)
  expect_equal(coef(res1), ref$coef, tolerance = 1e-6)
})

test_that("induced-logit survey.design works for clustered + stratified designs", {
  skip_if_not_installed("survey")

  set.seed(9203)
  n_strata <- 4
  psu_per_stratum <- 6
  n_per_psu <- 20
  stratum <- rep(seq_len(n_strata), each = psu_per_stratum * n_per_psu)
  psu <- rep(seq_len(n_strata * psu_per_stratum), each = n_per_psu)
  n <- length(psu)

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.2 + 0.8 * x1 - 0.3 * x2 + rnorm(n)

  p <- plogis(-0.5 + 0.2 * scale(y)[, 1] + 0.15 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(
    y_miss = ifelse(r == 1L, y, NA_real_),
    x1 = x1,
    x2 = x2,
    stratum = factor(stratum),
    psu = factor(psu)
  )
  w <- runif(n, 0.5, 2)
  df$w <- w

  des <- survey::svydesign(ids = ~psu, strata = ~stratum, weights = ~w, data = df)
  f <- y_miss ~ x1 + x2 | x1

  eng0 <- induced_logit_engine(variance_method = "none", standardize = FALSE, on_failure = "error")
  eng1 <- induced_logit_engine(variance_method = "none", standardize = TRUE, on_failure = "error")
  res0 <- nmar(f, data = des, engine = eng0)
  res1 <- nmar(f, data = des, engine = eng1)

  expect_true(res0$converged)
  expect_true(res1$converged)
  expect_true(is.finite(res0$y_hat))
  expect_true(is.finite(res1$y_hat))
  expect_equal(res1$y_hat, res0$y_hat, tolerance = 1e-8)
  expect_equal(fitted(res1), fitted(res0), tolerance = 1e-8)
  expect_equal(coef(res1), coef(res0), tolerance = 1e-6)
  expect_equal(res0$inference$df, survey::degf(des))
})

test_that("induced-logit survey estimate equals direct weighted Eq. (13) reconstruction", {
  skip_if_not_installed("survey")

  set.seed(92031)
  n <- 360
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  p <- plogis(-0.6 + 0.25 * y + 0.15 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  w <- runif(n, 0.5, 2)
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2, w = w)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  eng <- induced_logit_engine(variance_method = "none", on_failure = "error", keep_fits = TRUE)
  res <- nmar(y_miss ~ x1 + x2 | x1, data = des, engine = eng)
  expect_true(res$converged)

  spec <- res$extra$raw$spec
  mu_fit <- res$extra$raw$mu_fit
  glm_fit <- res$extra$raw$induced_glm

  vars <- des$variables
  w_full <- as.numeric(stats::weights(des))
  mu_hat <- as.numeric(stats::predict(mu_fit, newdata = vars))
  eps_hat <- as.numeric(vars$y_miss[spec$respondent_mask] - mu_hat[spec$respondent_mask])
  gamma_hat <- as.numeric(-stats::coef(glm_fit)[["..nmar_mu_hat.."]])
  eta_hat <- sum(w_full * as.integer(spec$respondent_mask)) / sum(w_full)
  w_resp <- w_full[spec$respondent_mask]
  m2_over_m1 <- sum(w_resp * eps_hat * exp(gamma_hat * eps_hat)) / sum(w_resp * exp(gamma_hat * eps_hat))
  tau_manual <- sum(w_full * mu_hat) / sum(w_full) + (1 - eta_hat) * m2_over_m1

  expect_equal(res$y_hat, tau_manual, tolerance = 1e-10)
  expect_true(is.finite(res$diagnostics$response_model_rank))
  expect_true(is.finite(res$diagnostics$response_model_ncol))
  expect_true(is.finite(res$diagnostics$response_model_condition_number))
  expect_equal(res$diagnostics$response_model_rank, res$diagnostics$response_model_ncol)
})

test_that("induced-logit survey design policy enforces supported assumptions", {
  skip_if_not_installed("survey")

  set.seed(92032)
  n <- 260
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.5 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.1 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  w <- runif(n, 0.5, 2)
  df <- data.frame(
    y_miss = ifelse(r == 1L, y, NA_real_),
    x1 = x1,
    x2 = x2,
    w = w,
    fpc = rep(n + 100L, n)
  )
  des_fpc <- survey::svydesign(ids = ~1, weights = ~w, fpc = ~fpc, data = df)
  f <- y_miss ~ x1 + x2 | x1

  eng_strict <- induced_logit_engine(variance_method = "none", on_failure = "error", survey_design_policy = "strict")
  expect_error(
    nmar(f, data = des_fpc, engine = eng_strict),
    "weight-based estimating equations"
  )

  eng_warn <- induced_logit_engine(variance_method = "none", on_failure = "error", survey_design_policy = "warn")
  res_warn <- NULL
  expect_warning(
    {
      res_warn <- nmar(f, data = des_fpc, engine = eng_warn)
    },
    "weight-based estimating equations"
  )
  expect_true(res_warn$converged)
  expect_equal(res_warn$diagnostics$survey_design_policy, "warn")
  expect_true(is.list(res_warn$diagnostics$survey_assumptions))

  des_base <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  pop <- c(`(Intercept)` = sum(w), x1 = sum(w * x1))
  des_cal <- survey::calibrate(des_base, formula = ~x1, population = pop)
  expect_error(
    nmar(f, data = des_cal, engine = eng_warn),
    "does not support calibrated/post-stratified designs"
  )
})

test_that("induced-logit survey path validates extracted design weights", {
  skip_if_not_installed("survey")

  set.seed(92033)
  n <- 180
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.1 + 0.5 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.5 + 0.2 * y + 0.1 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2, w = runif(n, 0.5, 2))
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  f <- y_miss ~ x1 + x2 | x1
  eng <- induced_logit_engine(variance_method = "none", on_failure = "error")

  des_neg <- des
  des_neg$prob <- rep(-1, n)
  expect_error(
    nmar(f, data = des_neg, engine = eng),
    "weights must be nonnegative"
  )

  des_na <- des
  des_na$prob <- c(NA_real_, rep(1, n - 1))
  expect_error(
    nmar(f, data = des_na, engine = eng),
    "weights must be finite and aligned"
  )

  des_zero <- des
  des_zero$prob <- rep(Inf, n)
  expect_error(
    nmar(f, data = des_zero, engine = eng),
    "weights must sum to a positive number"
  )
})

test_that("induced-logit survey on_failure='return' returns non-converged result on failure", {
  skip_if_not_installed("survey")

  set.seed(92034)
  n <- 200
  x1 <- rnorm(n)
  x_dup <- x1
  x2 <- rnorm(n)
  y <- 1 + 0.6 * x1 - 0.2 * x2 + rnorm(n)
  p <- plogis(-0.4 + 0.2 * y + 0.1 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x_dup = x_dup, x2 = x2, w = runif(n, 0.5, 2))
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  res <- induced_logit.survey.design(
    data = des,
    formula = y_miss ~ x1 + x2 | x1 + x_dup,
    variance_method = "none",
    on_failure = "return"
  )

  expect_s3_class(res, "nmar_result_induced_logit")
  expect_false(res$converged)
  expect_true(isTRUE(res$sample$is_survey))
  expect_true(is.character(res$diagnostics$message) && nzchar(res$diagnostics$message))
})

test_that("induced-logit survey.design supports factors and interactions under standardize", {
  skip_if_not_installed("survey")

  set.seed(9204)
  n <- 500
  x1f <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  x2 <- rnorm(n)
  x1n <- as.numeric(x1f) - 2
  y <- 0.5 + 0.3 * x1n + 0.6 * x2 + 0.4 * x1n * x2 + rnorm(n)

  p <- plogis(-0.7 + 0.2 * scale(y)[, 1] + 0.3 * x1n)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L
# Ensure outcome-model factor levels are present among respondents so the
# outcome regression is identifiable.
  for (lvl in levels(x1f)) {
    idx <- which(x1f == lvl)
    if (length(idx) >= 3) r[idx[1:3]] <- 1L
  }

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1f = x1f, x2 = x2)
  w <- runif(n, 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  f <- y_miss ~ x1f * x2 | x1f
  eng0 <- induced_logit_engine(variance_method = "none", standardize = FALSE, on_failure = "error")
  eng1 <- induced_logit_engine(variance_method = "none", standardize = TRUE, on_failure = "error")
  res0 <- nmar(f, data = des, engine = eng0)
  res1 <- nmar(f, data = des, engine = eng1)

  expect_true(res0$converged)
  expect_true(res1$converged)
  expect_equal(res1$y_hat, res0$y_hat, tolerance = 1e-8)
  expect_equal(res1$diagnostics$gamma_hat_paper, res0$diagnostics$gamma_hat_paper, tolerance = 1e-8)
  expect_equal(fitted(res1), fitted(res0), tolerance = 1e-8)
  expect_equal(coef(res1), coef(res0), tolerance = 1e-6)
})

test_that("induced-logit survey bootstrap SE is invariant to standardize", {
  skip_if_not_installed("survey")
  skip_if_not_installed("svrep")

  set.seed(9205)
  n <- 350
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1.4 + 0.7 * x1 - 0.3 * x2 + rnorm(n)

  p <- plogis(-0.6 + 0.2 * scale(y)[, 1] + 0.2 * x1)
  r <- rbinom(n, 1, p)
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L

  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
  w <- runif(n, 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  f <- y_miss ~ x1 + x2 | x1

  reps <- 20
  set.seed(9206)
  eng0 <- induced_logit_engine(variance_method = "bootstrap", bootstrap_reps = reps, standardize = FALSE, on_failure = "error")
  res0 <- nmar(f, data = des, engine = eng0)

  set.seed(9206)
  eng1 <- induced_logit_engine(variance_method = "bootstrap", bootstrap_reps = reps, standardize = TRUE, on_failure = "error")
  res1 <- nmar(f, data = des, engine = eng1)

  expect_true(res0$converged)
  expect_true(res1$converged)
  expect_true(is.finite(res0$se))
  expect_true(is.finite(res1$se))
  expect_equal(res1$y_hat, res0$y_hat, tolerance = 1e-8)
  expect_equal(res1$se, res0$se, tolerance = 1e-6)
})
