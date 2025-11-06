test_that("el_denominator computes guarded pack correctly", {
  lambda_W <- 2
  W <- 0.3
  p <- c(0.2, 0.4)
  Xc_lambda <- c(0.1, -0.1)
  floor <- 1.0

  dp <- NMAR:::el_denominator(lambda_W, W, Xc_lambda, p, floor)

  raw <- 1 + lambda_W * (p - W) + Xc_lambda
  expect_equal(dp$denom, pmax(raw, floor))
  expect_equal(dp$active, as.numeric(raw > floor))
  expect_equal(dp$inv, 1 / pmax(raw, floor))
  expect_equal(dp$inv_sq, dp$inv * dp$inv)
})

test_that("el_masses returns masses and probabilities (no trimming)", {
  weights <- c(2, 3)
  denom <- c(1, 2)
  out <- NMAR:::el_masses(weights, denom, floor = 1e-8, trim_cap = Inf)
  expect_equal(out$mass_untrim, weights / denom)
  expect_equal(sum(out$prob_mass), 1, tolerance = 1e-12)
  expect_equal(out$prob_mass, out$mass_trimmed / sum(out$mass_trimmed))
  expect_equal(out$trimmed_fraction, 0)
})

test_that("el_masses respects trimming cap and normalizes probabilities", {
# Feasible trimming that preserves total
  weights <- c(4, 1)
  denom <- c(1, 1)
  cap <- 3
  out <- NMAR:::el_masses(weights, denom, floor = 1e-8, trim_cap = cap)
  expect_true(all(out$mass_trimmed <= cap + 1e-12))
  expect_equal(sum(out$prob_mass), 1, tolerance = 1e-12)
  expect_true(out$trimmed_fraction >= 0)
})
