test_that("induced-logit input preprocessing parses partitions and defaults RHS2 to | 1", {
  set.seed(331)
  n <- 120
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x1 + x2 + rnorm(n)
  r <- rbinom(n, 1, plogis(-0.3 + 0.2 * y))
  if (all(r == 1L)) r[1] <- 0L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  inp1 <- induced_logit_prepare_inputs(y_miss ~ x1 + x2, data = df)
  expect_equal(inp1$outcome, "y_miss")
  expect_equal(paste(deparse(inp1$mu_rhs[[2L]]), collapse = ""), "x1 + x2")
  expect_equal(paste(deparse(inp1$x1_rhs[[2L]]), collapse = ""), "1")

  inp2 <- induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df)
  expect_equal(paste(deparse(inp2$x1_rhs[[2L]]), collapse = ""), "x1")

  spec_no_bar <- induced_logit_prepare_inputs(y_miss ~ x1 + x2, data = df)
  fit_no_bar <- il_fit_from_backend(spec = spec_no_bar, backend = il_backend_iid(df), standardize = FALSE)

  spec_bar_1 <- induced_logit_prepare_inputs(y_miss ~ x1 + x2 | 1, data = df)
  fit_bar_1 <- il_fit_from_backend(spec = spec_bar_1, backend = il_backend_iid(df), standardize = FALSE)

  expect_equal(fit_no_bar$point$tau_hat, fit_bar_1$point$tau_hat, tolerance = 1e-10)
})

test_that("induced-logit v1 rejects transformed LHS and dot expansion", {
  set.seed(332)
  n <- 80
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x1 + rnorm(n)
  r <- rbinom(n, 1, plogis(-0.2 + 0.2 * y))
  if (all(r == 1L)) r[1] <- 0L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  expect_error(induced_logit_prepare_inputs(scale(y_miss) ~ x1, data = df), "left-hand side")
  expect_error(induced_logit_prepare_inputs(y_miss ~ ., data = df), "`.` expansion")
  expect_error(induced_logit_prepare_inputs(y_miss ~ x1 | ., data = df), "`.` expansion")
})

test_that("induced-logit rejects all-respondent and all-nonrespondent samples", {
  set.seed(3321)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.3 * x1 - 0.2 * x2 + rnorm(n)

  df_all_resp <- data.frame(y_miss = y, x1 = x1, x2 = x2)
  expect_error(
    induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df_all_resp),
    "No nonrespondents detected"
  )

  df_all_nonresp <- data.frame(y_miss = rep(NA_real_, n), x1 = x1, x2 = x2)
  expect_error(
    induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df_all_nonresp),
    "No respondents detected"
  )
})

test_that("induced-logit rejects non-numeric outcomes", {
  set.seed(3322)
  n <- 40
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  y_fac <- factor(sample(c("a", "b"), n, replace = TRUE))
  y_log <- sample(c(TRUE, FALSE), n, replace = TRUE)

  df_fac <- data.frame(y_miss = y_fac, x1 = x1, x2 = x2)
  df_log <- data.frame(y_miss = y_log, x1 = x1, x2 = x2)

  expect_error(induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df_fac), "Outcome must be numeric")
  expect_error(induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df_log), "Outcome must be numeric")
})

test_that("induced-logit missingness model must include intercept", {
  set.seed(3323)
  n <- 90
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x1 + x2 + rnorm(n)
  r <- rbinom(n, 1, plogis(-0.3 + 0.2 * y))
  if (all(r == 1L)) r[1] <- 0L
  if (all(r == 0L)) r[1] <- 1L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  expect_error(
    induced_logit_prepare_inputs(y_miss ~ x1 + x2 | 0 + x1, data = df),
    "must include an intercept"
  )
})

test_that("induced-logit standardize must be a TRUE/FALSE scalar", {
  expect_error(induced_logit_engine(standardize = NA), "standardize")
  expect_error(induced_logit_engine(standardize = 1), "standardize")
})

test_that("induced-logit survey_design_policy is validated and defaults to strict", {
  eng_default <- induced_logit_engine()
  expect_equal(eng_default$survey_design_policy, "strict")
  expect_error(induced_logit_engine(survey_design_policy = "bad"), "one of")
})

test_that("induced-logit keep_fits is validated and defaults to FALSE", {
  eng_default <- induced_logit_engine()
  expect_false(eng_default$keep_fits)
  expect_error(induced_logit_engine(keep_fits = NA), "keep_fits")
  expect_error(induced_logit_engine(keep_fits = 1), "keep_fits")
})

test_that("induced-logit input validation enforces fully observed and finite covariates", {
  set.seed(333)
  n <- 90
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + x1 + x2 + rnorm(n)
  r <- rbinom(n, 1, plogis(-0.2 + 0.2 * y))
  if (all(r == 1L)) r[1] <- 0L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)

  df_na <- df
  df_na$x2[5] <- NA_real_
  expect_error(induced_logit_prepare_inputs(y_miss ~ x1 + x2 | x1, data = df_na), "fully observed")

  df_inf <- df
  df_inf$x1[1] <- 0
# log(x1^2) = -Inf when x1==0 (but no NaN), should be rejected as non-finite
  expect_error(induced_logit_prepare_inputs(y_miss ~ log(x1^2) + x2 | x2, data = df_inf), "finite")
})

test_that("induced-logit rejects outcome on RHS", {
  set.seed(334)
  n <- 60
  x1 <- rnorm(n)
  y <- 1 + x1 + rnorm(n)
  r <- rbinom(n, 1, plogis(-0.1 + 0.1 * y))
  if (all(r == 1L)) r[1] <- 0L
  df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1)

  expect_error(induced_logit_prepare_inputs(y_miss ~ x1 + y_miss, data = df), "Outcome cannot appear")
  expect_error(induced_logit_prepare_inputs(y_miss ~ x1 | y_miss, data = df), "Outcome cannot appear")
})

test_that("induced-logit rejects mu-model factor levels that occur only among nonrespondents", {
  set.seed(335)
  n <- 80
  x1 <- rnorm(n)
  grp <- rep(c("A", "B"), each = n / 2)
  y <- 1 + 0.3 * x1 + rnorm(n)

# Make all "B" units nonrespondents
  r <- ifelse(grp == "A", 1L, 0L)
  df <- data.frame(
    y_miss = ifelse(r == 1L, y, NA_real_),
    x1 = x1,
    grp = factor(grp)
  )

  expect_error(
    induced_logit_prepare_inputs(y_miss ~ grp + x1 | x1, data = df),
    "factor levels must be observed among respondents"
  )
})
