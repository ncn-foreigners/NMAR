test_that("EL weights normalize correctly (probability and population)", {
  set.seed(9101)
  df <- make_iid_nmar(n = 300, alpha = 0.5, seed = 9101)

  eng <- make_engine(auxiliary_means = c(X = 0), variance_method = "none", standardize = TRUE, trim_cap = Inf)
  fit <- nmar(Y_miss ~ X, data = df, engine = eng)
  expect_true(fit$converged)

# Probability weights sum to 1
  wp <- weights(fit, scale = "probability")
  expect_true(abs(sum(wp) - 1) < 1e-12)

# Population weights sum to N_pop
  wN <- weights(fit, scale = "population")
  N_pop <- fit$sample$n_total
  expect_true(abs(sum(wN) - N_pop) < 1e-8)

# Relationship holds numerically (ignore attributes)
  expect_equal(unname(as.numeric(wN)), as.numeric(N_pop) * unname(as.numeric(wp)), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("Trimming caps EL masses and normalization still holds", {
  set.seed(9102)
  df <- make_iid_nmar(n = 350, alpha = 0.8, seed = 9102)

# Pilot fit without trimming to estimate a cap that activates trimming
  eng0 <- make_engine(auxiliary_means = c(X = 0), variance_method = "none", standardize = TRUE, trim_cap = Inf)
  fit0 <- nmar(Y_miss ~ X, data = df, engine = eng0)
  info0 <- NMAR:::nmar_result_get_weights_info(fit0)
  base_masses <- as.numeric(info0$values)
  expect_true(length(base_masses) > 0)

# Choose a cap at the 90th percentile of untrimmed masses to induce non-zero trimming
  cap <- as.numeric(stats::quantile(base_masses, probs = 0.90, na.rm = TRUE))

  eng <- make_engine(auxiliary_means = c(X = 0), variance_method = "none", standardize = TRUE, trim_cap = cap)
  fit <- nmar(Y_miss ~ X, data = df, engine = eng)
  expect_true(fit$converged)

# Extract trimmed (unnormalized) masses from result internals
  info <- NMAR:::nmar_result_get_weights_info(fit)
  trimmed_masses <- as.numeric(info$values)
# At least some trimming occurred
  expect_true((info$trimmed_fraction %||% 0) > 0)
# All trimmed masses are <= cap (within tiny tolerance)
  expect_true(max(trimmed_masses) <= cap + 1e-8)

# Normalization invariants still hold
  wp <- weights(fit, scale = "probability")
  expect_true(abs(sum(wp) - 1) < 1e-12)
  wN <- weights(fit, scale = "population")
  N_pop <- fit$sample$n_total
  expect_true(abs(sum(wN) - N_pop) < 1e-8)
  expect_equal(unname(as.numeric(wN)), as.numeric(N_pop) * unname(as.numeric(wp)), tolerance = 1e-12, ignore_attr = TRUE)
})
