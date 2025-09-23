test_that("standardize vs manual scaling yields identical y_hat (df)", {
  set.seed(3403)
  N <- 200
  X <- rnorm(N)
  Y <- 2 + 0.5 * X + rnorm(N)
  p <- plogis(-1 + 0.4 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_

  # Auto scaling
  fit_auto <- nmar:::el.data.frame(df, Y_miss ~ X,
    response_predictors = NULL,
    auxiliary_means = c(X = 0), standardize = TRUE,
    variance_method = "delta"
  )
  expect_true(fit_auto$converged)

  # Manual scaling of response predictor (Y_miss) and auxiliary X, then unscale y_hat
  resp_df <- df[!is.na(df$Y_miss), ]
  Z_un <- model.matrix(~Y_miss, data = resp_df)
  X_un <- model.matrix(~ X - 1, data = resp_df)
  recipe <- nmar:::create_nmar_scaling_recipe(Z_un, X_un)
  df_scaled <- df
  # Scale both Y_miss and X in the data (affects response model and outcome)
  if ("Y_miss" %in% names(recipe)) {
    df_scaled$Y_miss <- (df_scaled$Y_miss - recipe$Y_miss$mean) / recipe$Y_miss$sd
  }
  if ("X" %in% names(recipe) && is.numeric(df_scaled$X)) {
    df_scaled$X <- (df_scaled$X - recipe$X$mean) / recipe$X$sd
  }
  aux_mean_scaled <- c(X = 0)
  if ("X" %in% names(recipe)) aux_mean_scaled["X"] <- (0 - recipe$X$mean) / recipe$X$sd

  fit_manual <- nmar:::el.data.frame(df_scaled, Y_miss ~ X,
    response_predictors = NULL,
    auxiliary_means = aux_mean_scaled, standardize = FALSE,
    variance_method = "delta"
  )
  expect_true(fit_manual$converged)
  # Unscale y_hat back to original Y_miss scale
  y_hat_unscaled <- fit_manual[['estimate']] * recipe$Y_miss$sd + recipe$Y_miss$mean
  expect_equal(fit_auto[['estimate']], y_hat_unscaled, tolerance = 1e-8)
})

test_that("multi-predictor scaling matches manual rescaling", {
  set.seed(1234)
  N <- 150
  X1 <- rnorm(N)
  X2 <- runif(N, -1, 1)
  Y <- 1.5 + 0.7 * X1 - 0.3 * X2 + rnorm(N)
  p <- plogis(-0.2 + 0.5 * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X1 = X1, X2 = X2)
  df[!R, "Y_miss"] <- NA_real_

  fit_auto <- nmar:::el.data.frame(df, Y_miss ~ X1 + X2,
    response_predictors = NULL,
    auxiliary_means = c(X1 = 0, X2 = 0), standardize = TRUE,
    variance_method = "delta"
  )
  expect_true(fit_auto$converged)

  resp_df <- df[!is.na(df$Y_miss), ]
  Z_un <- model.matrix(~ Y_miss + X1 + X2, data = resp_df)
  X_un <- model.matrix(~ X1 + X2 - 1, data = resp_df)
  recipe <- nmar:::create_nmar_scaling_recipe(Z_un, X_un)
  df_scaled <- df
  for (nm in names(recipe)) {
    if (nm %in% names(df_scaled)) {
      df_scaled[[nm]] <- (df_scaled[[nm]] - recipe[[nm]]$mean) / recipe[[nm]]$sd
    }
  }
  aux_means_scaled <- c(
    X1 = (0 - recipe$X1$mean) / recipe$X1$sd,
    X2 = (0 - recipe$X2$mean) / recipe$X2$sd
  )

  fit_manual <- nmar:::el.data.frame(df_scaled, Y_miss ~ X1 + X2,
    response_predictors = NULL,
    auxiliary_means = aux_means_scaled, standardize = FALSE,
    variance_method = "delta"
  )
  expect_true(fit_manual$converged)
  y_hat_manual <- fit_manual[["estimate"]]
  y_hat_unscaled <- y_hat_manual * recipe$Y_miss$sd + recipe$Y_miss$mean
  expect_equal(fit_auto[["estimate"]], y_hat_unscaled, tolerance = 1e-8)
})
