test_that("Coverage harness re-draws missingness across replicates", {
  skip_on_cran()
  if (Sys.getenv("NMAR_RUN_INTEGRATION", unset = "0") != "1") skip("integration-only")
  if (!file.exists("inst/simulations/el_solver_compare.R")) skip("harness missing")
  source("inst/simulations/el_solver_compare.R")
  out <- run_coverage_targeted_multi(B = 12L, boot_reps = 60L)
  expect_true(all(out$var_emp > 0, na.rm = TRUE))
})
