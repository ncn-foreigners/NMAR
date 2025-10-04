test_that("logit family functions behave numerically", {
  fam <- NMAR:::logit_family()
  eta <- c(-10, -2, 0, 2, 10)
  p <- fam$linkinv(eta)
  expect_true(all(p > 0 & p < 1))
  expect_true(all(is.finite(fam$mu.eta(eta))))
  expect_true(all(is.finite(fam$d2mu.deta2(eta))))
  sc <- fam$score_eta(eta, 1)
  expect_equal(length(sc), length(eta))
})

test_that("choose_jacobian returns reasonable structure", {
# Simple quadratic function f(x) = Ax; Jacobian is A.
  A <- matrix(c(2, 0, 0, 3), 2, 2)
  f <- function(x) as.numeric(A %*% x)
  res <- NMAR:::choose_jacobian(analytic_fun = function(x) A, numeric_fun = function(x) numDeriv::jacobian(f, x), at = c(1, 1))
  expect_true(is.matrix(res$A))
  expect_true(res$source %in% c("analytic", "numeric"))
})

test_that("trim_weights caps and redistributes mass", {
  w <- c(5, 1, 1, 1)
  cap <- 2
  res <- NMAR:::trim_weights(w, cap)
  expect_equal(max(res$weights), cap)
  expect_equal(sum(res$weights), sum(w))
})

test_that("trim_weights warns when all mass is capped", {
  w <- c(5, 5)
  cap <- 1
  expect_warning(res <- NMAR:::trim_weights(w, cap), "reduced total weight")
  expect_equal(res$weights, rep(cap, length(w)))
})

test_that("enforce_nonneg_weights flags large negative weights", {
  w <- c(0.2, -1e-6, -0.3)
  res <- NMAR:::enforce_nonneg_weights(w, tol = 1e-3)
  expect_false(res$ok)
  expect_equal(res$weights, w)
})

test_that("enforce_nonneg_weights clips small negative weights", {
  w <- c(0.2, -1e-10)
  res <- NMAR:::enforce_nonneg_weights(w, tol = 1e-3)
  expect_true(res$ok)
  expect_equal(res$weights, c(0.2, 0))
})
