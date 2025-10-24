test_that("exptilt scaling yields near-invariant estimates (IID)", {
  skip_if(Sys.getenv("NMAR_RUN_INTEGRATION") != "1", "Set NMAR_RUN_INTEGRATION=1 to run integration tests.")

  set.seed(7317)
  n <- 800
  x1 <- rnorm(n)
  x2 <- runif(n, -1, 1)
  y_true <- 0.8 + 0.5 * x1 - 0.25 * x2 + rnorm(n, sd = 0.6)
  p_obs <- plogis(-0.3 + 0.5 * scale(y_true)[, 1])
  delta <- rbinom(n, 1, p_obs)
  y_obs <- y_true
  y_obs[delta == 0] <- NA_real_
  dat <- data.frame(y = y_obs, x1 = x1, x2 = x2)

# Fit with automatic standardization
  fit_std <- nmar(
    y ~ x1 + x2,
    data = dat,
    engine = exptilt_engine(
      standardize = TRUE,
      y_dens = "normal",
      variance_method = "delta",
      control = list(maxit = 30),
      stopping_threshold = 1e-6
    ),
    response_predictors = NULL
  )
  expect_true(fit_std$converged)

# Manual scaling recipe from respondents, then fit without standardize
  obs_mask <- !is.na(dat$y)
  Z_un <- model.matrix(~ y + x1 + x2, data = dat[obs_mask, ])
  X_un <- model.matrix(~ x1 + x2 - 1, data = dat[obs_mask, ])
  recipe <- NMAR:::create_nmar_scaling_recipe(Z_un, X_un)
  dat_scaled <- dat
  for (nm in names(recipe)) {
    dat_scaled[[nm]] <- (dat_scaled[[nm]] - recipe[[nm]]$mean) / recipe[[nm]]$sd
  }

  fit_raw <- nmar(
    y ~ x1 + x2,
    data = dat_scaled,
    engine = exptilt_engine(
      standardize = FALSE,
      y_dens = "normal",
      variance_method = "delta",
      control = list(maxit = 30),
      stopping_threshold = 1e-6
    ),
    response_predictors = NULL
  )
  expect_true(fit_raw$converged)

# Unscale the raw fit back to the original Y scale
  est_raw_unscaled <- as.numeric(fit_raw$y_hat * recipe$y$sd + recipe$y$mean)
  se_raw_unscaled <- as.numeric(fit_raw$se * recipe$y$sd)

# Compare with reasonable tolerances (iterative solver may land within small numerical deltas)
  est_rel_diff <- abs(fit_std$y_hat - est_raw_unscaled) / max(1, abs(fit_std$y_hat))
  se_rel_diff <- abs(fit_std$se - se_raw_unscaled) / fit_std$se

  expect_lt(est_rel_diff, 0.01) # < 1% relative difference in point estimate
  expect_lt(se_rel_diff, 0.20) # < 20% relative difference in SE
})
