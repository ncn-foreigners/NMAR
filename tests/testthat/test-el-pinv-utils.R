test_that("invert_jacobian uses pseudoinverse when requested", {
  A <- matrix(c(1, 1, 1, 1), 2, 2) # rank-1 singular
  expect_error(nmar:::invert_jacobian(A, variance_ridge = FALSE, variance_pseudoinverse = FALSE))
  res <- nmar:::invert_jacobian(A, variance_pseudoinverse = TRUE)
  expect_true(is.list(res) && is.matrix(res$inv))
  expect_true(isTRUE(res$used_pinv))
})
