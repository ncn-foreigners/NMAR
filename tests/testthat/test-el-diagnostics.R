test_that("EL diagnostics: constraints and Jacobian quality", {
  set.seed(123)
  N <- 600
  x <- stats::rchisq(N, df = 2)
  eps <- stats::rnorm(N)
  y <- x + eps * sqrt(x) / 5
  pr <- stats::plogis(0.2 * y - (-2))
  r <- stats::rbinom(N, 1, pr)
  df <- data.frame(y_miss = ifelse(r == 1, y, NA_real_), x = x)

  eng <- make_engine(
    auxiliary_means = c(x = mean(df$x)),
    variance_method = "delta",
    standardize = TRUE,
    solver_args = list(global = "dbldog"),
    control = list(maxit = 200, xtol = 1e-8, ftol = 1e-8)
  )
  fit <- nmar(y_miss ~ x, data = df, engine = eng)
  expect_true(isTRUE(fit$converged))
  diag <- fit$diagnostics
# Constraint sums (equations) should be near zero
  expect_true(is.finite(diag$constraint_sum_W))
  expect_lt(abs(diag$constraint_sum_W), 1e-5)
  if (!is.null(diag$constraint_sum_aux) && length(diag$constraint_sum_aux)) {
    expect_true(all(is.finite(diag$constraint_sum_aux)))
    expect_lt(max(abs(diag$constraint_sum_aux)), 1e-5)
  }
# Analytic vs numeric Jacobian should agree at the solution
  if (is.finite(diag$jacobian_rel_diff)) {
    expect_lt(diag$jacobian_rel_diff, 1e-3)
  }
})
