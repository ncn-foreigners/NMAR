interaction_fixture <- function(seed = 123, n = 50) {
  set.seed(seed)
  data.frame(
    Y_miss = c(rnorm(ceiling(n * 0.6)), rep(NA_real_, floor(n * 0.4))),
    X = rnorm(n),
    Z = rnorm(n)
  )
}

test_that("auxiliary interactions and transforms propagate through nmar()", {
  df <- interaction_fixture()
  eng <- el_engine(variance_method = "none")

  fit_interaction <- nmar(Y_miss ~ X * Z, data = df, engine = eng)
  expect_s3_class(fit_interaction, "nmar_result_el")
  expect_true(fit_interaction$converged)
  expect_equal(length(coef(fit_interaction)), 2)
  expect_true(all(c("X", "Z", "X:Z") %in% names(fit_interaction$diagnostics$constraint_sum_aux)))

  fit_explicit <- nmar(Y_miss ~ X + Z + X:Z, data = df, engine = eng)
  expect_true(fit_explicit$converged)
  expect_equal(length(coef(fit_explicit)), 2)
  expect_setequal(names(fit_explicit$diagnostics$constraint_sum_aux), c("X", "Z", "X:Z"))

  fit_identity <- nmar(Y_miss ~ X + I(X^2), data = df, engine = eng)
  expect_true(fit_identity$converged)
  expect_equal(length(coef(fit_identity)), 2)
  expect_setequal(names(fit_identity$diagnostics$constraint_sum_aux), c("X", "I(X^2)"))

  fit_poly <- nmar(Y_miss ~ poly(X, 2), data = df, engine = eng)
  expect_true(fit_poly$converged)
  expect_equal(length(coef(fit_poly)), 2)
  poly_names <- names(fit_poly$diagnostics$constraint_sum_aux)
  expect_equal(length(poly_names), 2)
  expect_true(all(grepl("poly", poly_names)))
})

test_that("partitioned formulas propagate interactions to auxiliaries and response predictors", {
  df <- interaction_fixture(seed = 126)
  df$X1 <- df$X
  df$X2 <- rnorm(nrow(df))
  df$Z2 <- rnorm(nrow(df))
  eng <- el_engine(variance_method = "none")

  fit_partition <- nmar(Y_miss ~ X1 * X2 | Z2, data = df, engine = eng)
  expect_true(fit_partition$converged)
  beta <- coef(fit_partition)
  expect_true(all(c("Y_miss", "Z2") %in% names(beta)))
  expect_equal(length(beta), 3)
  expect_true(all(c("X1", "X2", "X1:X2") %in% names(fit_partition$diagnostics$constraint_sum_aux)))
})
