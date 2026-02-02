test_that("EL handles factor auxiliaries and response predictors; intercept policy enforced", {
  set.seed(101)
  n <- 80
  G <- factor(sample(letters[1:3], n, replace = TRUE))
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 2 + 0.5 * X + 0.2 * (as.numeric(G) - 1) + rnorm(n)
  p <- plogis(-0.3 + 0.4 * scale(Y)[, 1] + 0.2 * X)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z, G = G)
  df$Y_miss[!R] <- NA_real_

# Case 1: auxiliary_means supplied -> allow NA in factor auxiliary among nonrespondents
  df2 <- df
  nonresp <- which(is.na(df2$Y_miss))[1]
  df2$G[nonresp] <- NA
  eng <- el_engine(auxiliary_means = c(X = 0, Ga = 0, Gb = 0, Gc = 0), variance_method = "none")
# Use a model where RHS expands factor(G) in auxiliaries; extra auxiliary_means
# names (e.g., Ga) are ignored with a warning.
  expect_warning(
    fit <- nmar(Y_miss ~ X + G | Z, data = df2, engine = eng),
    regexp = "Ignoring unused names in 'auxiliary_means'",
    fixed = FALSE
  )
  expect_s3_class(fit, "nmar_result_el")

# Case 2: no auxiliary_means -> NA columns trigger error
  eng2 <- el_engine(auxiliary_means = NULL, variance_method = "none")
  expect_error(
    nmar(Y_miss ~ X + G | Z, data = df2, engine = eng2),
    regexp = "Auxiliary variables contain NA values",
    fixed = TRUE
  )
})
