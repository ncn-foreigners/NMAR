test_that("delta uses analytic ∇g when untrimmed and numeric when trimmed (indirect sensitivity test)", {
# Idea: when trim_cap = Inf (smooth), analytic gradient path is used, so
# tuning numeric gradient options should not affect the SE. When trimming is
# active, the numeric gradient path is used and SE becomes sensitive to those
# options. We compare SEs under two very different numeric gradient settings.
  on.exit({options(nmar.grad_eps = NULL, nmar.grad_d = NULL)}, add = TRUE)
  df <- make_iid_nmar(n = 200, alpha = 0.5, seed = 123)
  aux_means <- c(X = 0)

# Untrimmed: should be insensitive to grad_numeric options (analytic ∇g)
  options(nmar.grad_eps = 1e-6, nmar.grad_d = 1e-3)
  eng_inf_1 <- make_engine(variance_method = "delta", family = "logit",
                           auxiliary_means = aux_means, trim_cap = Inf,
                           standardize = TRUE)
  fit_inf_1 <- nmar(Y_miss ~ X, data = df, engine = eng_inf_1)
  se_inf_1 <- fit_inf_1$std_error

  options(nmar.grad_eps = 1e-2, nmar.grad_d = 1e-1)
  eng_inf_2 <- make_engine(variance_method = "delta", family = "logit",
                           auxiliary_means = aux_means, trim_cap = Inf,
                           standardize = TRUE)
  fit_inf_2 <- nmar(Y_miss ~ X, data = df, engine = eng_inf_2)
  se_inf_2 <- fit_inf_2$std_error

# Expect near equality under smooth case
  expect_equal(se_inf_1, se_inf_2, tolerance = 1e-6)

# Trimmed: numeric gradient path; sensitivity to grad options likely
  options(nmar.grad_eps = 1e-6, nmar.grad_d = 1e-3)
  eng_trim_1 <- make_engine(variance_method = "delta", family = "logit",
                            auxiliary_means = aux_means, trim_cap = 1.1,
                            standardize = TRUE)
  fit_trim_1 <- nmar(Y_miss ~ X, data = df, engine = eng_trim_1)
  se_trim_1 <- fit_trim_1$std_error
# Ensure trimming engaged
  expect_gt(fit_trim_1$diagnostics$trimmed_fraction, 0)

  options(nmar.grad_eps = 1e-2, nmar.grad_d = 1e-1)
  eng_trim_2 <- make_engine(variance_method = "delta", family = "logit",
                            auxiliary_means = aux_means, trim_cap = 1.1,
                            standardize = TRUE)
  fit_trim_2 <- nmar(Y_miss ~ X, data = df, engine = eng_trim_2)
  se_trim_2 <- fit_trim_2$std_error

# Expect larger difference under trimmed (non-smooth) case
  expect_gt(abs(se_trim_1 - se_trim_2), 1e-5)
})
