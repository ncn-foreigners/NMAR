test_that("el_resolve_auxiliaries works for data.frame with level drops", {
  set.seed(1)
  n <- 20
  full <- data.frame(
    Y = rnorm(n),
    f = factor(sample(c("A", "B"), n, replace = TRUE))
  )
  aux_formula <- ~ f - 1
  auxiliary_design_full <- model.matrix(aux_formula, data = full)
  respondent_mask <- full$f == "A"
  out <- NMAR:::el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = NULL)
  expect_true(is.matrix(out$auxiliary_design))
  expect_true(length(out$means) >= 1)
  expect_setequal(colnames(out$auxiliary_design), names(out$means))
  expect_equal(nrow(out$auxiliary_design), sum(respondent_mask))
})

test_that("el_resolve_auxiliaries computes design-weighted means for survey.design", {
  skip_if_not_installed("survey")
  set.seed(2)
  n <- 30
  full <- data.frame(
    Y = rnorm(n),
    f = factor(sample(c("A", "B"), n, replace = TRUE)),
    w = runif(n, 0.5, 2)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = full)
  aux_formula <- ~ f - 1
  auxiliary_design_full <- model.matrix(aux_formula, data = des$variables)
  respondent_mask <- des$variables$f == "A"
  out <- NMAR:::el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = NULL, weights_full = weights(des))
  mm_full <- model.matrix(aux_formula, data = full)
  mu_expected <- as.numeric(colSums(mm_full * full$w) / sum(full$w))
  names(mu_expected) <- colnames(mm_full)
# Restrict to returned columns
  mu_expected <- mu_expected[colnames(out$auxiliary_design)]
  expect_equal(unname(out$means), unname(mu_expected), tolerance = 1e-12)
})

test_that("el_resolve_auxiliaries warns on extra names in auxiliary_means and ignores them", {
  set.seed(3)
  n <- 25
  full <- data.frame(
    Y = rnorm(n),
    X = rnorm(n)
  )
  auxiliary_design_full <- model.matrix(~ X - 1, data = full)
  respondent_mask <- rep(TRUE, n)
  aux_means_supplied <- c(X = 0, EXTRA = 123)
  out <- expect_warning(
    NMAR:::el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = aux_means_supplied),
    regexp = "Ignoring unused names in 'auxiliary_means'",
    fixed = FALSE
  )
  expect_equal(names(out$means), colnames(out$auxiliary_design))
  expect_false("EXTRA" %in% names(out$means))
})
