test_that("solver presets converge and record effective top-level args", {
  df <- make_iid_nmar(n = 120, alpha = 0.4, seed = 8201)
  fml <- Y_miss ~ X
  table <- list(
    list(name = "newton-dbldog", solver_method = "newton",
         solver_jacobian = "analytic",
         solver_args = list(global = "dbldog"),
         control = list(xtol = 1e-8, ftol = 1e-8, maxit = 100),
         expect_method = "Newton", expect_global = "dbldog", expect_xscalm = NA_character_),
    list(name = "broyden-fixed", solver_method = "broyden",
         solver_jacobian = "none",
         solver_args = list(xscalm = "fixed"),
         control = list(maxit = 200),
         expect_method = "Broyden", expect_global = NA_character_, expect_xscalm = "fixed"),
    list(name = "auto-qline", solver_method = "auto",
         solver_jacobian = "analytic",
         solver_args = list(global = "qline"),
         control = list(maxit = 100),
         expect_method = NULL, expect_global = "qline", expect_xscalm = NA_character_)
  )

  for (row in table) {
    eng <- make_engine(
      variance_method = "none",
      auxiliary_means = c(X = 0),
      solver_method = row$solver_method,
      solver_args = row$solver_args,
      control = row$control,
      standardize = TRUE
    )
    fit <- nmar(fml, data = df, engine = eng)
    expect_true(fit$converged, info = row$name)
# Check recorded method (for auto, allow either Newton or Broyden)
    if (!is.null(row$expect_method)) {
      expect_equal(fit$diagnostics$solver_method, row$expect_method, info = row$name)
    } else {
      expect_true(fit$diagnostics$solver_method %in% c("Newton", "Broyden"), info = row$name)
    }
# Check top-level nleqslv args recorded
    if (is.na(row$expect_global)) {
      expect_true(is.na(fit$diagnostics$nleqslv_global) || is.null(fit$diagnostics$nleqslv_global))
    } else {
      expect_equal(fit$diagnostics$nleqslv_global, row$expect_global, info = row$name)
    }
    if (is.na(row$expect_xscalm)) {
      expect_true(is.na(fit$diagnostics$nleqslv_xscalm) || is.null(fit$diagnostics$nleqslv_xscalm))
    } else {
      expect_equal(fit$diagnostics$nleqslv_xscalm, row$expect_xscalm, info = row$name)
    }
  }
})
