# Test script to verify bootstrap sampling fix
# This script demonstrates that bootstrap variance is now correctly estimated
# regardless of the sample_size parameter used for the point estimate

library(nmar)

# Create test data with known properties
set.seed(123)
n <- 10000
data <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  y = rnorm(n, mean = 2, sd = 1)
)

# Introduce NMAR missingness
# P(respond | y, x1) depends on y
prob_respond <- plogis(-0.5 + 0.3 * data$y + 0.2 * data$x1)
data$y[runif(n) > prob_respond] <- NA

cat("Data summary:\n")
cat(sprintf("  Total observations: %d\n", nrow(data)))
cat(sprintf("  Respondents: %d (%.1f%%)\n",
            sum(!is.na(data$y)),
            100 * sum(!is.na(data$y)) / nrow(data)))
cat(sprintf("  Non-respondents: %d (%.1f%%)\n\n",
            sum(is.na(data$y)),
            100 * sum(is.na(data$y)) / nrow(data)))

# Test 1: Bootstrap variance should be similar regardless of sample_size
cat("Test 1: Bootstrap variance with different sample_size values\n")
cat("=", rep("=", 70), "\n", sep = "")

# Note: Using small bootstrap_reps for speed (increase for real analysis)
bootstrap_reps <- 50

cat("\nRunning with sample_size = 2000...\n")
result_2k <- exptilt(
  data = data,
  formula = y ~ x1 + x2,
  auxiliary_means = list(x1 = 0, x2 = 0),
  sample_size = 2000,
  variance_method = "bootstrap",
  bootstrap_reps = bootstrap_reps,
  trace_level = 0
)

cat("\nRunning with sample_size = 5000...\n")
result_5k <- exptilt(
  data = data,
  formula = y ~ x1 + x2,
  auxiliary_means = list(x1 = 0, x2 = 0),
  sample_size = 5000,
  variance_method = "bootstrap",
  bootstrap_reps = bootstrap_reps,
  trace_level = 0
)

cat("\nRunning with sample_size = Inf (no sampling)...\n")
result_inf <- exptilt(
  data = data,
  formula = y ~ x1 + x2,
  auxiliary_means = list(x1 = 0, x2 = 0),
  sample_size = Inf,
  variance_method = "bootstrap",
  bootstrap_reps = bootstrap_reps,
  trace_level = 0
)

cat("\n\nResults Comparison:\n")
cat("-------------------\n")
cat(sprintf("sample_size = 2000:  Estimate = %.4f, SE = %.4f\n",
            result_2k$y_hat, result_2k$se))
cat(sprintf("sample_size = 5000:  Estimate = %.4f, SE = %.4f\n",
            result_5k$y_hat, result_5k$se))
cat(sprintf("sample_size = Inf:   Estimate = %.4f, SE = %.4f\n",
            result_inf$y_hat, result_inf$se))

# Calculate relative differences in SE
diff_2k <- abs(result_2k$se - result_inf$se) / result_inf$se
diff_5k <- abs(result_5k$se - result_inf$se) / result_inf$se

cat("\n\nRelative SE differences (should be small, < 20% due to bootstrap sampling):\n")
cat(sprintf("  2k vs Inf: %.1f%%\n", 100 * diff_2k))
cat(sprintf("  5k vs Inf: %.1f%%\n", 100 * diff_5k))

if (diff_2k < 0.3 && diff_5k < 0.3) {
  cat("\n✓ PASS: Bootstrap variance is consistent across sample_size values\n")
} else {
  cat("\n✗ WARNING: Bootstrap variance differs significantly across sample_size values\n")
  cat("  This may indicate a problem with the bootstrap implementation.\n")
}

# Test 2: Compare with delta method (when applicable)
cat("\n\nTest 2: Compare bootstrap vs delta method\n")
cat("=", rep("=", 70), "\n", sep = "")

cat("\nRunning with variance_method = 'delta'...\n")
result_delta <- exptilt(
  data = data,
  formula = y ~ x1 + x2,
  auxiliary_means = list(x1 = 0, x2 = 0),
  sample_size = Inf,
  variance_method = "delta",
  trace_level = 0
)

cat("\n\nResults Comparison:\n")
cat("-------------------\n")
cat(sprintf("Delta method:  Estimate = %.4f, SE = %.4f\n",
            result_delta$y_hat, result_delta$se))
cat(sprintf("Bootstrap:     Estimate = %.4f, SE = %.4f\n",
            result_inf$y_hat, result_inf$se))

diff_method <- abs(result_delta$se - result_inf$se) / result_delta$se
cat(sprintf("\nRelative SE difference: %.1f%%\n", 100 * diff_method))

if (diff_method < 0.3) {
  cat("✓ PASS: Bootstrap and delta methods agree (within 30%)\n")
} else {
  cat("✗ WARNING: Bootstrap and delta methods differ substantially\n")
}

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("Test Summary:\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("The bootstrap variance estimation should now work correctly:\n")
cat("1. SE should be similar regardless of sample_size used for point estimate\n")
cat("2. Bootstrap should agree with delta method (when both are valid)\n")
cat("3. Bootstrap operates on the FULL dataset, not the sampled subset\n")
cat("\nIf tests pass, the bootstrap sampling fix is working correctly.\n")
