test_that("el_resolve_auxiliaries rejects NA in full auxiliary data", {
  auxiliary_design_full <- rbind(
    c(X = 1),
    c(X = NA_real_),
    c(X = 0)
  )
  respondent_mask <- c(TRUE, FALSE, TRUE)
  expect_error(
    NMAR:::el_resolve_auxiliaries(auxiliary_design_full, respondent_mask, auxiliary_means = NULL),
    "Auxiliary variables contain NA values",
    fixed = FALSE
  )
})
