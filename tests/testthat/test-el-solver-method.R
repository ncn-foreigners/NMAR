test_that("solver_method routing is honored in diagnostics", {
  set.seed(11)
  n <- 40
  Y_full <- rbinom(n, 1, 0.5)
  R <- rbinom(n, 1, 0.7)
  Y_miss <- Y_full
  Y_miss[R == 0] <- NA
  df <- data.frame(Y = Y_miss, X = rnorm(n))
# Force Newton only
  fit_newton <- NMAR:::el.data.frame(df, Y ~ X,
    standardize = TRUE,
    variance_method = "none",
    solver_method = "newton",
    solver_jacobian = "analytic"
  )
  expect_true(fit_newton$converged)
  expect_equal(fit_newton$diagnostics$solver_method, "Newton")

# Force Broyden only
  fit_broyden <- NMAR:::el.data.frame(df, Y ~ X,
    standardize = TRUE,
    variance_method = "none",
    solver_method = "broyden",
    solver_jacobian = "none",
    solver_args = list(global = "dbldog", xscalm = "fixed"),
    control = list(global = "qline")
  )
  expect_true(fit_broyden$converged)
  expect_equal(fit_broyden$diagnostics$solver_method, "Broyden")
# solver_args should take precedence over control-provided top-level args
  expect_equal(fit_broyden$diagnostics$nleqslv_global, "dbldog")
  expect_equal(fit_broyden$diagnostics$nleqslv_xscalm, "fixed")
})

test_that("solver_args take precedence for survey designs as well", {
  skip_if_not_installed("survey")
  set.seed(7)
  df <- make_iid_nmar(n = 80, alpha = 0.4, seed = 7)
  w <- runif(nrow(df), 0.5, 2)
  des <- survey::svydesign(ids = ~1, weights = ~w, data = transform(df, w = w))
  fit <- NMAR:::el.survey.design(des, Y_miss ~ X,
    standardize = TRUE,
    variance_method = "none",
    solver_method = "broyden",
    solver_jacobian = "none",
    solver_args = list(global = "dbldog", xscalm = "auto"),
    control = list(global = "qline")
  )
  expect_true(fit$converged)
  expect_equal(fit$diagnostics$solver_method, "Broyden")
  expect_equal(fit$diagnostics$nleqslv_global, "dbldog")
  expect_equal(fit$diagnostics$nleqslv_xscalm, "auto")
})
