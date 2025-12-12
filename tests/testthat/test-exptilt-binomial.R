test_that("exptilt works with binomial outcome", {
  set.seed(42)
  n <- 120
  x1 <- rnorm(n)
# generate latent outcome and observed outcome with missingness depending on y
  lin <- -0.3 + 0.7 * x1
  y_true <- rbinom(n, 1, plogis(lin))
  respond <- rbinom(n, 1, plogis(0.5 + 0.4 * y_true))
  if (all(respond == 1)) respond[sample.int(n, 1)] <- 0
  y_obs <- ifelse(respond == 1, y_true, NA)
  dat <- data.frame(Y = y_obs, x1 = x1)

  eng <- exptilt_engine(
    y_dens = 'binomial',
    family = 'probit',
    control = list(maxit = 6),
    stopping_threshold = 0.01
  )

  res <- nmar(formula = Y ~ x1, data = dat, engine = eng)

  expect_s3_class(res, "nmar_result")
  expect_type(res, "list")
# primary estimand stored in y_hat
  expect_true(is.finite(res$y_hat))
# standard error should be numeric (may be NA if variance_method='none')
  expect_true(is.numeric(res$se) || is.na(res$se))
})
