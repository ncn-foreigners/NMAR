test_that("el_resolve_auxiliaries works for data.frame with level drops", {
  set.seed(1)
  n <- 20
  full <- data.frame(
    Y = rnorm(n),
    f = factor(sample(c("A", "B"), n, replace = TRUE))
  )
  resp <- subset(full, f == "A")
  aux_formula <- ~ f - 1
  aux_resp <- model.matrix(aux_formula, data = resp)
  aux_full <- model.matrix(aux_formula, data = full)
  out <- NMAR:::el_resolve_auxiliaries(aux_resp, aux_full, auxiliary_means = NULL)
  expect_true(is.matrix(out$matrix))
  expect_true(length(out$means) >= 1)
  expect_setequal(colnames(out$matrix), names(out$means))
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
  resp <- subset(full, f == "A")
  aux_formula <- ~ f - 1
  aux_resp <- model.matrix(aux_formula, data = resp)
  aux_full <- model.matrix(aux_formula, data = des$variables)
  out <- NMAR:::el_resolve_auxiliaries(aux_resp, aux_full, auxiliary_means = NULL, weights_full = weights(des))
  mm_full <- model.matrix(aux_formula, data = full)
  mu_expected <- as.numeric(colSums(mm_full * full$w) / sum(full$w))
  names(mu_expected) <- colnames(mm_full)
# Restrict to returned columns
  mu_expected <- mu_expected[colnames(out$matrix)]
  expect_equal(unname(out$means), unname(mu_expected), tolerance = 1e-12)
})
