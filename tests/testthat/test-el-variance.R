test_that("bootstrap and delta SE differ on typical df case", {
  set.seed(3201)
  df <- make_iid_nmar(n = 500, alpha = 0.4, seed = 3201)
  fml <- Y_miss ~ X

  fit_d <- suppressWarnings(nmar(formula = fml, data = df, engine = make_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = FALSE)))
  fit_b <- nmar(formula = fml, data = df, engine = make_engine(auxiliary_means = c(X = 0), variance_method = "bootstrap", bootstrap_reps = 50, standardize = FALSE))

  expect_true(fit_d$converged && fit_b$converged)
  expect_true(is.na(fit_d[['se']]) || is.finite(fit_d[['se']]))
  expect_true(is.na(fit_b[['se']]) || is.finite(fit_b[['se']]))
  if (isTRUE(is.finite(fit_b[['se']]))) {
    expect_false(isTRUE(all.equal(fit_d[['se']], fit_b[['se']])))
  }
})

test_that("bootstrap path runs under finite trim without delta warning assertion", {
  set.seed(3202)
  df <- make_iid_nmar(n = 400, alpha = 0.6, seed = 3202)
  fml <- Y_miss ~ X

  prom2 <- testthat::evaluate_promise(
    nmar(formula = fml, data = df, engine = make_engine(auxiliary_means = c(X = 0), variance_method = "bootstrap", trim_cap = 5, bootstrap_reps = 20, standardize = FALSE))
  )
  expect_true(length(prom2$warnings) == 0 || !any(grepl("Delta method variance is not recommended", prom2$warnings)))
})
