test_that("trim_weights caps and preserves mass when possible", {
  w <- c(0.1, 0.2, 0.3, 10)
# Choose cap high enough so total capacity (n*cap) >= sum(w)
  cap <- 4.0
  out <- NMAR:::trim_weights(w, cap)
  wt <- out$weights
# All weights <= cap
  expect_true(all(wt <= cap + 1e-12))
# Mass preserved because cap is not binding overall
  expect_equal(sum(wt), sum(w), tolerance = 1e-10)
# Fraction at cap should be 1/length(w) in this simple case (only the big one)
  expect_equal(out$trimmed_fraction, 1 / length(w))
})

test_that("trim_weights warns when mass cannot be preserved", {
  w <- c(5, 5, 5)
  cap <- 1.0
  expect_warning({
    out <- NMAR:::trim_weights(w, cap)
  })
  wt <- out$weights
  expect_true(all(wt <= cap + 1e-12))
# In this extreme case no eligible mass remains to redistribute to, so total decreases
  expect_lt(sum(wt), sum(w))
# All are at cap
  expect_equal(out$trimmed_fraction, 1)
})
