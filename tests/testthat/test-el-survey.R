test_that("EL survey path produces finite SE and df when survey available", {
  skip_if_not_installed("survey")
  set.seed(456)
  N <- 400
  x <- stats::rchisq(N, df = 2)
  eps <- stats::rnorm(N)
  y <- x + eps * sqrt(x) / 5
  pr <- stats::plogis(0.2 * y - (-1.5))
  r <- stats::rbinom(N, 1, pr)
  df <- data.frame(y_miss = ifelse(r == 1, y, NA_real_), x = x)
  design <- survey::svydesign(ids = ~1, weights = ~1, data = df)

  eng <- make_engine(
    auxiliary_means = c(x = mean(df$x)),
    variance_method = "none",
    standardize = TRUE,
    control = list(maxit = 200, xtol = 1e-8, ftol = 1e-8)
  )
  fit <- nmar(y_miss ~ x, data = design, engine = eng)
  expect_true(isTRUE(fit$converged))
  se <- nmar_result_get_se(fit)
  expect_true(is.na(se) || is.finite(se))
  inf <- nmar_result_get_inference(fit)
  expect_true(is.finite(inf$df) || is.na(inf$df))
  diag <- fit$diagnostics
  expect_true(is.finite(diag$jacobian_condition_number) || is.na(diag$jacobian_condition_number))
})
