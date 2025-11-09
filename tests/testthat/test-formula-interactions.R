test_that("formula interactions work end-to-end", {
  set.seed(123)
  n <- 50
  df <- data.frame(
    Y_miss = c(rnorm(30), rep(NA, 20)),
    X = rnorm(n),
    Z = rnorm(n)
  )

# Interaction with * (auxiliary side: X, Z, X:Z)
  eng <- el_engine(variance_method = "none")
  fit_interaction <- nmar(Y_miss ~ X * Z, data = df, engine = eng)

  expect_s3_class(fit_interaction, "nmar_result")
  expect_s3_class(fit_interaction, "nmar_result_el")
  expect_true(fit_interaction$converged)

# coef() returns response model params (delta ~ Y_miss only)
  beta <- coef(fit_interaction)
  expect_equal(length(beta), 2) # Intercept + Y_miss

# Check that interaction was used in auxiliaries
  aux_constraints <- fit_interaction$diagnostics$constraint_sum_aux
  expect_true("X:Z" %in% names(aux_constraints))
  expect_true("X" %in% names(aux_constraints))
  expect_true("Z" %in% names(aux_constraints))
})

test_that("polynomial formulas work via poly()", {
  set.seed(124)
  n <- 50
  df <- data.frame(
    Y_miss = c(rnorm(30), rep(NA, 20)),
    X = rnorm(n)
  )

  eng <- el_engine(variance_method = "none")
  fit_poly <- nmar(Y_miss ~ poly(X, 2), data = df, engine = eng)

  expect_s3_class(fit_poly, "nmar_result")
  expect_s3_class(fit_poly, "nmar_result_el")
  expect_true(fit_poly$converged)

# coef() returns response model params (delta ~ Y_miss only)
  beta <- coef(fit_poly)
  expect_equal(length(beta), 2) # Intercept + Y_miss

# Check that poly() was used in auxiliaries
  aux_constraints <- fit_poly$diagnostics$constraint_sum_aux
  expect_equal(length(aux_constraints), 2) # poly(X, 2) creates 2 columns
  expect_true(all(grepl("poly", names(aux_constraints))))
})

test_that("explicit interaction via : works", {
  set.seed(125)
  n <- 50
  df <- data.frame(
    Y_miss = c(rnorm(30), rep(NA, 20)),
    X = rnorm(n),
    Z = rnorm(n)
  )

  eng <- el_engine(variance_method = "none")
  fit_explicit <- nmar(Y_miss ~ X + Z + X:Z, data = df, engine = eng)

  expect_s3_class(fit_explicit, "nmar_result")
  expect_true(fit_explicit$converged)

# coef() returns response model params (delta ~ Y_miss only)
  beta <- coef(fit_explicit)
  expect_equal(length(beta), 2) # Intercept + Y_miss

# Check that explicit interaction was used in auxiliaries
  aux_constraints <- fit_explicit$diagnostics$constraint_sum_aux
  expect_true("X:Z" %in% names(aux_constraints))
  expect_true("X" %in% names(aux_constraints))
  expect_true("Z" %in% names(aux_constraints))
  expect_equal(length(aux_constraints), 3) # X, Z, X:Z
})

test_that("interactions work with partitioned formula", {
  set.seed(126)
  n <- 50
  df <- data.frame(
    Y_miss = c(rnorm(30), rep(NA, 20)),
    X1 = rnorm(n),
    X2 = rnorm(n),
    Z = rnorm(n)
  )

  eng <- el_engine(variance_method = "none")
  fit_partition <- nmar(Y_miss ~ X1 * X2 | Z, data = df, engine = eng)

  expect_s3_class(fit_partition, "nmar_result")
  expect_true(fit_partition$converged)

# Response model: delta ~ Y_miss + Z (with intercept)
  beta <- coef(fit_partition)
  expect_true("Z" %in% names(beta))
  expect_true("Y_miss" %in% names(beta))
  expect_equal(length(beta), 3) # Intercept + Y_miss + Z

# Check that interaction was used in auxiliaries
  aux_constraints <- fit_partition$diagnostics$constraint_sum_aux
  expect_true("X1:X2" %in% names(aux_constraints))
  expect_true("X1" %in% names(aux_constraints))
  expect_true("X2" %in% names(aux_constraints))
})

test_that("I() identity wrapper works in formulas", {
  set.seed(127)
  n <- 50
  df <- data.frame(
    Y_miss = c(rnorm(30), rep(NA, 20)),
    X = rnorm(n)
  )

  eng <- el_engine(variance_method = "none")
  fit_identity <- nmar(Y_miss ~ X + I(X^2), data = df, engine = eng)

  expect_s3_class(fit_identity, "nmar_result")
  expect_true(fit_identity$converged)

# coef() returns response model params (delta ~ Y_miss only)
  beta <- coef(fit_identity)
  expect_equal(length(beta), 2) # Intercept + Y_miss

# Check that I(X^2) was used in auxiliaries
  aux_constraints <- fit_identity$diagnostics$constraint_sum_aux
  expect_true("I(X^2)" %in% names(aux_constraints))
  expect_true("X" %in% names(aux_constraints))
  expect_equal(length(aux_constraints), 2) # X and I(X^2)
})
