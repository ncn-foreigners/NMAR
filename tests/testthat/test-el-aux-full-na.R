test_that("el_resolve_auxiliaries rejects NA in full auxiliary data when means missing", {
  auxiliary_design_full <- rbind(
    c(X = 1),
    c(X = NA_real_),
    c(X = 0)
  )
  respondent_mask <- c(TRUE, FALSE, TRUE)
  expect_error(
    NMAR:::el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = NULL),
    regexp = "Auxiliary variables contain NA values",
    fixed = TRUE
  )
})

test_that("el_resolve_auxiliaries allows NA among nonrespondents when auxiliary_means supplied", {
  auxiliary_design_full <- rbind(
    c(X = 1),
    c(X = NA_real_), # nonrespondent with missing X
    c(X = 0)
  )
  respondent_mask <- c(TRUE, FALSE, TRUE)
  out <- NMAR:::el_resolve_auxiliaries(
    auxiliary_design_full,
    respondent_mask = respondent_mask,
    auxiliary_means = c(X = 0.5)
  )
  expect_equal(nrow(out$auxiliary_design), sum(respondent_mask))
  expect_equal(colnames(out$auxiliary_design), "X")
  expect_equal(out$means[["X"]], 0.5)
})
