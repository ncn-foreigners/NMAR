test_that("EL survey vs iid equivalence for iid-like design", {
  skip_if_not_installed("survey")

  set.seed(123)
  df <- make_iid_nmar(n = 300, alpha = 0.4, include_z = TRUE, seed = 1L)

  eng <- make_engine(
    variance_method = "none",
    auxiliary_means = c(X = mean(df$X), Z = mean(df$Z)),
    standardize = TRUE,
    trim_cap = Inf
  )

  fit_df <- nmar(Y_miss ~ X + Z | X, data = df, engine = eng, trace_level = 0)

  design <- survey::svydesign(ids = ~1, weights = ~1, data = df)
  fit_svy <- nmar(Y_miss ~ X + Z | X, data = design, engine = eng, trace_level = 0)

  expect_true(isTRUE(fit_df$converged))
  expect_true(isTRUE(fit_svy$converged))

  expect_equal(fit_df$y_hat, fit_svy$y_hat, tolerance = 1e-6)

  w_df <- as.numeric(weights(fit_df, scale = "probability"))
  w_svy <- as.numeric(weights(fit_svy, scale = "probability"))
  expect_equal(w_df, w_svy, tolerance = 1e-6)
})


test_that("survey EL constraints (including Wu strata auxiliaries) are near zero", {
  skip_if_not_installed("survey")

  set.seed(321)
  N <- 400
  strata <- factor(rep(letters[1:4], length.out = N))
  X <- stats::rnorm(N)
  Y <- stats::rnorm(N)
  R <- stats::rbinom(N, size = 1, prob = 0.8)

  df <- data.frame(
    Y_miss = ifelse(R == 1, Y, NA_real_),
    X = X,
    strata = strata
  )

  design <- survey::svydesign(ids = ~1, weights = ~1, strata = ~strata, data = df)

  aux_means <- c(X = mean(df$X))
  eng <- make_engine(
    variance_method = "none",
    auxiliary_means = aux_means,
    standardize = TRUE,
    trim_cap = Inf
  )

  fit <- nmar(Y_miss ~ X | X, data = design, engine = eng, trace_level = 0)
  expect_true(isTRUE(fit$converged))

  diag <- fit$diagnostics
  sumw <- sum(stats::weights(design))

  expect_lt(abs(diag$constraint_sum_W) / sumw, 1e-6)

  ca <- diag$constraint_sum_aux
  expect_true(is.numeric(ca))
  if (length(ca) > 0) {
    expect_lt(max(abs(ca) / sumw), 1e-6)
    strata_idx <- grepl("^strata_", names(ca))
    if (any(strata_idx)) {
      expect_lt(max(abs(ca[strata_idx]) / sumw), 1e-6)
    }
  }
})
