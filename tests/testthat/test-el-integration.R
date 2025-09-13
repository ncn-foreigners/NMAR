set.seed(2025)

test_that("EL engine runs and returns expected structure (data.frame)", {
  N <- 400
  X <- rnorm(N)
  Z <- rnorm(N)
  Y <- 2 + 0.5 * X + Z
  p <- plogis(-1.0 + 0.4 * scale(Y)[, 1])
  R <- runif(N) < p
  dat <- data.frame(Y_miss = Y, X = X)
  dat$Y_miss[!R] <- NA

  eng <- el_engine(variance_method = "delta", auxiliary_means = c(X = 0))
  fml <- list(outcome = ~Y_miss, covariates_outcome = ~X, covariates_missingness = ~NULL)
  res <- nmar(formula = fml, data = dat, engine = eng)

  expect_s3_class(res, "nmar_result_el")
  expect_true(isTRUE(res$converged))
  expect_true(is.numeric(res$y_hat))
  expect_true(is.numeric(res$se))

  est <- res$y_hat
  expect_equal(as.numeric(est), res$y_hat)

  ci <- nmar:::confint.nmar_result_el(res)
  expect_true(is.matrix(ci) && nrow(ci) == 1)
  expect_equal(rownames(ci), "Y_miss")

  w <- weights(res)
  expect_true(is.numeric(w) || length(w) == 0)
})
