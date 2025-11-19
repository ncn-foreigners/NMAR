test_that("validate_and_apply_nmar_scaling handles weighted standardization", {
  Z <- cbind(`(Intercept)` = 1, x1 = c(1, 2, 3), x2 = c(2, 5, 6))
  weights <- c(1, 2, 3)
  scaling <- validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = FALSE,
    response_model_matrix_unscaled = Z,
    auxiliary_matrix_unscaled = matrix(nrow = nrow(Z), ncol = 0),
    mu_x_unscaled = NULL,
    weights = weights
  )
  Z_scaled <- scaling$response_model_matrix_scaled
  expect_equal(Z_scaled[, "(Intercept)"], rep(1, 3))
  weighted_mean_x1 <- sum(weights * Z_scaled[, "x1"]) / sum(weights)
  weighted_mean_x2 <- sum(weights * Z_scaled[, "x2"]) / sum(weights)
  expect_equal(weighted_mean_x1, 0, tolerance = 1e-12)
  expect_equal(weighted_mean_x2, 0, tolerance = 1e-12)
  weighted_sd <- function(col) sqrt(sum(weights * col^2) / sum(weights))
  expect_equal(weighted_sd(Z_scaled[, "x1"]), 1, tolerance = 1e-12)
  expect_equal(weighted_sd(Z_scaled[, "x2"]), 1, tolerance = 1e-12)
})

test_that("validate_and_apply_nmar_scaling warns on constant predictors", {
  Z <- cbind(`(Intercept)` = 1, x1 = c(1, 1, 1), x2 = c(0, 1, 2))
  expect_warning(
    scaling <- validate_and_apply_nmar_scaling(
      standardize = TRUE,
      has_aux = FALSE,
      response_model_matrix_unscaled = Z,
      auxiliary_matrix_unscaled = matrix(nrow = nrow(Z), ncol = 0),
      mu_x_unscaled = NULL
    ),
    "nearly constant",
    fixed = FALSE
  )
  Z_scaled <- scaling$response_model_matrix_scaled
  expect_equal(Z_scaled[, "x1"], rep(0, 3))
})

test_that("unscale_coefficients warns when recipe incomplete (identity fallback)", {
  coeffs <- c(`(Intercept)` = 1, x1 = 0.5, x2 = -0.2)
  vcov <- diag(c(0.04, 0.01, 0.09))
  recipe <- structure(list(x1 = list(mean = 0, sd = 1)), class = "nmar_scaling_recipe")
  expect_warning(
    res <- unscale_coefficients(coeffs, vcov, recipe),
    "missing entries",
    fixed = FALSE
  )
# Identity fallback should preserve inputs for missing item(s)
  expect_equal(res$coefficients, coeffs)
  expect_equal(res$vcov, vcov)
})
