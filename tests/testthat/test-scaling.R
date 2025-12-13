test_that("validate_and_apply_nmar_scaling handles weighted standardization", {
  Z <- cbind(`(Intercept)` = 1, x1 = c(1, 2, 3), x2 = c(2, 5, 6))
  weights <- c(1, 2, 3)
  scaling <- validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = FALSE,
    response_model_matrix_unscaled = Z,
    aux_matrix_unscaled = matrix(nrow = nrow(Z), ncol = 0),
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
      aux_matrix_unscaled = matrix(nrow = nrow(Z), ncol = 0),
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

test_that("scale_coefficients and unscale_coefficients are inverse transforms", {
  beta_unscaled <- c(`(Intercept)` = 2, x1 = 0.5, x2 = -0.3)
  recipe <- structure(list(
    x1 = list(mean = 1, sd = 2),
    x2 = list(mean = -2, sd = 4)
  ), class = "nmar_scaling_recipe")
  cols <- c("(Intercept)", "x1", "x2")

  beta_scaled <- scale_coefficients(beta_unscaled, recipe, cols)
  vcov_scaled <- diag(rep(1, length(cols)))

  roundtrip <- unscale_coefficients(beta_scaled, vcov_scaled, recipe)

  expect_equal(roundtrip$coefficients, beta_unscaled, tolerance = 1e-10)
})

test_that("scale_aux_multipliers preserves X' lambda under scaling", {
  set.seed(123)
  X_un <- matrix(rnorm(20), nrow = 5)
  colnames(X_un) <- c("x1", "x2", "x3", "x4")
  lambda_un <- c(x1 = 0.5, x2 = -0.2, x3 = 0.7, x4 = -0.1)
  mu_x_un <- colMeans(X_un)

  recipe <- create_nmar_scaling_recipe(X_un)
  X_scaled <- apply_nmar_scaling(X_un, recipe)
  lambda_scaled <- scale_aux_multipliers(lambda_un, recipe, colnames(X_un))

  mu_x_scaled <- vapply(names(mu_x_un), function(n) {
    (mu_x_un[[n]] - recipe[[n]]$mean) / recipe[[n]]$sd
  }, numeric(1))

  Xc_un <- sweep(X_un, 2, mu_x_un, "-")
  Xc_scaled <- sweep(X_scaled, 2, mu_x_scaled, "-")

  t_un <- as.numeric(Xc_un %*% lambda_un)
  t_scaled <- as.numeric(Xc_scaled %*% lambda_scaled)

  expect_equal(t_un, t_scaled, tolerance = 1e-10)
})

test_that("prepare_nmar_scaling matches validate_and_apply_nmar_scaling behavior", {
  Z_un <- cbind(`(Intercept)` = 1, z1 = c(1, 2, 3), z2 = c(4, 5, 6))
  X_un <- cbind(x1 = c(7, 8, 9), x2 = c(10, 11, 12))
  mu_x <- c(x1 = mean(X_un[, "x1"]), x2 = mean(X_un[, "x2"]))
  w <- c(1, 2, 3)

  sc_low <- prepare_nmar_scaling(Z_un, X_un, mu_x, standardize = TRUE, weights = w)
  sc_high <- validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = TRUE,
    response_model_matrix_unscaled = Z_un,
    aux_matrix_unscaled = X_un,
    mu_x_unscaled = mu_x,
    weights = w
  )

  expect_equal(sc_low$Z, sc_high$response_model_matrix_scaled, tolerance = 1e-10)
  expect_equal(sc_low$X, sc_high$auxiliary_matrix_scaled, tolerance = 1e-10)
  expect_equal(sc_low$mu_x, sc_high$mu_x_scaled, tolerance = 1e-10)
})

test_that("weight_mask produces same scaling as pre-masked weights", {
  Z_un <- cbind(`(Intercept)` = 1, x1 = c(1, 2, 3, 4), x2 = c(5, 6, 7, 8))
  w <- c(1, 2, 3, 4)
  mask <- c(TRUE, FALSE, TRUE, FALSE)

  sc_masked_weights <- validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = FALSE,
    response_model_matrix_unscaled = Z_un,
    aux_matrix_unscaled = matrix(nrow = nrow(Z_un), ncol = 0),
    mu_x_unscaled = NULL,
    weights = w * mask
  )

  sc_weight_mask <- validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = FALSE,
    response_model_matrix_unscaled = Z_un,
    aux_matrix_unscaled = matrix(nrow = nrow(Z_un), ncol = 0),
    mu_x_unscaled = NULL,
    weights = w,
    weight_mask = mask
  )

  expect_equal(
    sc_masked_weights$response_model_matrix_scaled,
    sc_weight_mask$response_model_matrix_scaled,
    tolerance = 1e-10
  )
})

test_that("numeric weight_mask acts as multipliers (not boolean)", {
  Z_un <- cbind(`(Intercept)` = 1, x1 = c(1, 2, 3, 4), x2 = c(5, 6, 7, 8))
  w <- c(1, 2, 3, 4)
  mask <- c(1, 10, 1, 10)

  sc_mask <- validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = FALSE,
    response_model_matrix_unscaled = Z_un,
    aux_matrix_unscaled = matrix(nrow = nrow(Z_un), ncol = 0),
    mu_x_unscaled = NULL,
    weights = w,
    weight_mask = mask
  )

  sc_mult <- validate_and_apply_nmar_scaling(
    standardize = TRUE,
    has_aux = FALSE,
    response_model_matrix_unscaled = Z_un,
    aux_matrix_unscaled = matrix(nrow = nrow(Z_un), ncol = 0),
    mu_x_unscaled = NULL,
    weights = w * mask
  )

  expect_equal(
    sc_mask$response_model_matrix_scaled,
    sc_mult$response_model_matrix_scaled,
    tolerance = 1e-10
  )
})

test_that("all-zero effective weights after weight_mask errors", {
  Z_un <- cbind(`(Intercept)` = 1, x1 = c(1, 2, 3))
  w <- c(1, 1, 1)
  mask <- c(FALSE, FALSE, FALSE)

  expect_error(
    validate_and_apply_nmar_scaling(
      standardize = TRUE,
      has_aux = FALSE,
      response_model_matrix_unscaled = Z_un,
      aux_matrix_unscaled = matrix(nrow = nrow(Z_un), ncol = 0),
      mu_x_unscaled = NULL,
      weights = w,
      weight_mask = mask
    ),
    regexp = "no positive entries",
    fixed = FALSE
  )
})
