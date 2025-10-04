test_that("invert_jacobian: plain vs ridge vs pinv toggles work", {
# Well-conditioned case: plain inverse
  A_plain <- diag(3)
  res_plain <- NMAR:::invert_jacobian(A_plain)
  expect_equal(res_plain$invert_rule, "plain")
  expect_false(isTRUE(res_plain$used_pinv))
  expect_false(isTRUE(res_plain$used_ridge))

# Singular case: ridge fallback
  A_sing <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2) # rank-1
  res_ridge <- NMAR:::invert_jacobian(A_sing, variance_ridge = TRUE)
  expect_equal(res_ridge$invert_rule, "ridge")
  expect_false(isTRUE(res_ridge$used_pinv))
  expect_true(isTRUE(res_ridge$used_ridge))
  expect_true(is.finite(res_ridge$kappa) || is.infinite(res_ridge$kappa))

# Singular case: pseudoinverse
  res_pinv <- NMAR:::invert_jacobian(A_sing, variance_pseudoinverse = TRUE)
  expect_equal(res_pinv$invert_rule, "pinv")
  expect_true(isTRUE(res_pinv$used_pinv))
  expect_false(isTRUE(res_pinv$used_ridge))
})

test_that("invert_jacobian supports svd_tol pseudo-inverse", {
  skip_if_not_installed("MASS")
  A_sing <- matrix(c(1, 1, 1, 1), 2, 2)
  res_tol <- NMAR:::invert_jacobian(A_sing, variance_pseudoinverse = TRUE, svd_tol = 1e-6)
  expect_equal(res_tol$invert_rule, "pinv")
  expect_true(is.matrix(res_tol$inv))
  expect_equal(dim(res_tol$inv), dim(A_sing))
  expect_true(isTRUE(res_tol$used_pinv))
})
