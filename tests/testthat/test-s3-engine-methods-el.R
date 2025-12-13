test_that("engine_name/engine_config/format work for EL engines", {
  eng <- el_engine(
    variance_method = "none",
    standardize = TRUE,
    trim_cap = Inf,
    family = "logit"
  )

  expect_equal(engine_name(eng), "empirical_likelihood")

  cfg <- engine_config(eng)
  expect_true(is.list(cfg))
  expect_false(inherits(cfg, "nmar_engine"))

  s <- format(eng)
  expect_true(is.character(s) && length(s) == 1L)
  expect_match(s, "Empirical Likelihood \\(EL\\)")
  expect_match(s, "standardize=TRUE")
  expect_match(s, "variance_method=none")
  expect_match(s, "family=logit")
  expect_match(s, "trim_cap=Inf")
})

test_that("engine print formatting covers bootstrap, trim_cap and start fields", {
  eng_boot <- el_engine(
    variance_method = "bootstrap",
    bootstrap_reps = 12,
    trim_cap = Inf,
    start = list(beta = c("(Intercept)" = 0), W = 0.5, lambda = c(X = 0))
  )

  out <- capture.output(print(eng_boot))
  expect_true(any(grepl("^NMAR engine: Empirical Likelihood \\(EL\\)$", out)))
  expect_true(any(grepl("variance_method:", out)))
  expect_true(any(grepl("bootstrap.*reps", out)))
  expect_true(any(grepl("trim_cap:", out)))
  expect_true(any(grepl("^start:", out)))

  s <- format(eng_boot)
  expect_match(s, "variance_method=bootstrap\\(12\\)")
})
