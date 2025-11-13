test_that("Survey weights normalization matches user-supplied n_total (mismatch triggers warning)", {
  skip_if_not_installed("survey")
  set.seed(9201)
  N <- 600
  x <- rnorm(N)
  y <- 1 + x + rnorm(N)
# NMAR response
  pr <- plogis(0.5 - 0.3 * y)
  R <- rbinom(N, 1, pr)
  df <- data.frame(y_miss = ifelse(R == 1, y, NA_real_), x = x)

# Create a survey design with rescaled weights (analysis-scale)
  base_w <- rep(1, N)
  d2 <- base_w / mean(base_w) # still 1s here, but keep structure
  des <- survey::svydesign(ids = ~1, weights = ~d2, data = df)

# Supply n_total that differs from sum(weights(design)) to trigger rescaling
  n_total <- N * 2 # mismatch by 100%
  eng <- el_engine(auxiliary_means = c(x = mean(df$x)), variance_method = "none", standardize = TRUE, n_total = n_total)
  expect_warning(
    fit <- nmar(y_miss ~ x, data = des, engine = eng, trace_level = 0),
    regexp = "Scale mismatch detected|Large scale mismatch",
    fixed = FALSE
  )

# Population weights from the result must sum to the supplied n_total
  wN <- weights(fit, scale = "population")
  expect_true(abs(sum(wN) - n_total) < 1e-8)
# Probability weights should always sum to 1
  wp <- weights(fit, scale = "probability")
  expect_true(abs(sum(wp) - 1) < 1e-12)
})
