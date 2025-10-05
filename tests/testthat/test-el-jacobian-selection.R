test_that("kappa ratio gate flips to numeric A when analytic is ill-conditioned", {
# Toy system: linear root at (1,2,0)
  f <- function(x) c(x[1] - 1, x[2] - 2, x[3])
# Analytic Jacobian intentionally near-singular (two identical rows)
  j_bad <- function(x) {
    A <- diag(3)
    A[3, ] <- A[1, ]
    A
  }
  at <- c(1, 2, 0)
  sel <- NMAR:::el_select_variance_jacobian(
    equation_system_func = f,
    analytical_jac_func = j_bad,
    estimates = at,
    variance_jacobian = "auto"
  )
  expect_equal(sel$A_source, "numeric")
  expect_true(sel$jacobian_auto_rule %in% c("kappa_ratio_high", "rel_diff_high"))
})
