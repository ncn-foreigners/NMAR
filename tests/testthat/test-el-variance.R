test_that("bootstrap and delta SE differ on typical df case", {
  set.seed(3201)
  N <- 500
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-1 + 0.4 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_
  fml <-  Y_miss ~ X

  fit_d <- nmar(formula = fml, data = df, engine = el_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = FALSE))
  fit_b <- nmar(formula = fml, data = df, engine = el_engine(auxiliary_means = c(X = 0), variance_method = "bootstrap", bootstrap_reps = 50, standardize = FALSE, suppress_warnings = TRUE))

  expect_true(fit_d$converged && fit_b$converged)
  expect_true(is.finite(fit_d$se) && is.finite(fit_b$se))
  expect_false(isTRUE(all.equal(fit_d$se, fit_b$se)))
})

test_that("bootstrap path runs under finite trim without delta warning assertion", {
  set.seed(3202)
  N <- 400
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-1 + 0.6 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_
  fml <- Y_miss ~ X

  prom2 <- testthat::evaluate_promise(
    nmar(formula = fml, data = df, engine = el_engine(auxiliary_means = c(X = 0), variance_method = "bootstrap", trim_cap = 5, bootstrap_reps = 20, standardize = FALSE, suppress_warnings = TRUE))
  )
  expect_true(length(prom2$warnings) == 0 || !any(grepl("Delta method variance is not recommended", prom2$warnings)))
})
