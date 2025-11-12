test_that("el_resolve_auxiliaries rejects NA in full auxiliary data", {
  aux_resp <- cbind(X = rnorm(3))
  aux_full <- rbind(
    c(X = 1),
    c(X = NA_real_),
    c(X = 0)
  )
  expect_error(
    NMAR:::el_resolve_auxiliaries(aux_resp, aux_full, auxiliary_means = NULL),
    "Auxiliary variables contain NA values",
    fixed = FALSE
  )
})
