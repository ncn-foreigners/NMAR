# Monte Carlo validation for EL estimator

simulate_df <- function(N, alpha, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- rnorm(N)
  Z <- rnorm(N)
  Y <- 2 + 0.5 * X + Z
  p <- plogis(-1.0 + alpha * scale(Y)[, 1])
  R <- runif(N) < p
  df <- data.frame(Y_miss = Y, X = X)
  df[!R, "Y_miss"] <- NA_real_
  df
}

fit_el_df <- function(df, standardize = FALSE) {
  NMAR::nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = NMAR::el_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = standardize)
  )
}

fit_el_df_bootstrap <- function(df, standardize = FALSE, bootstrap_reps = 200L) {
  NMAR::nmar(
    formula = Y_miss ~ X,
    data = df,
    engine = NMAR::el_engine(auxiliary_means = c(X = 0), variance_method = "bootstrap", bootstrap_reps = bootstrap_reps, standardize = standardize)
  )
}

run_mc_once <- function(N, alpha, i) {
  df <- simulate_df(N, alpha)
  fit <- try(fit_el_df(df), silent = TRUE)
  if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
    return(list(ok = FALSE))
  }
  list(ok = TRUE, y_hat = fit$estimate, se = fit$se)
}

run_mc <- function(N, alpha, reps = 20L, true_mean = 2.0, seed = 1L) {
  set.seed(seed)
  vals <- replicate(reps, run_mc_once(N, alpha, 0L), simplify = FALSE)
  ok <- vapply(vals, function(v) isTRUE(v$ok), logical(1))
  y <- vapply(vals[ok], function(v) v$y_hat, numeric(1))
  se <- vapply(vals[ok], function(v) v$se, numeric(1))
  conv <- sum(ok)
  out <- list(
    N = N, alpha = alpha, reps = reps, converged = conv,
    mean_y = if (conv) mean(y) else NA_real_,
    bias = if (conv) mean(y) - true_mean else NA_real_,
    emp_sd = if (conv) stats::sd(y) else NA_real_,
    mean_se = if (conv) mean(se) else NA_real_,
    sd_se_ratio = if (conv) (stats::sd(y) / max(mean(se), 1e-8)) else NA_real_
  )
  class(out) <- "el_mc_summary"
  out
}

run_mc_bootstrap <- function(N, alpha, reps = 20L, true_mean = 2.0, seed = 1L, bootstrap_reps = 200L) {
  set.seed(seed)
  vals <- replicate(reps, {
    df <- simulate_df(N, alpha)
    fit <- try(fit_el_df_bootstrap(df, bootstrap_reps = bootstrap_reps), silent = TRUE)
    if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
      list(ok = FALSE)
    } else {
      list(ok = TRUE, y_hat = fit$estimate, se = fit$se)
    }
  }, simplify = FALSE)
  ok <- vapply(vals, function(v) isTRUE(v$ok), logical(1))
  y <- vapply(vals[ok], function(v) v$y_hat, numeric(1))
  se <- vapply(vals[ok], function(v) v$se, numeric(1))
  conv <- sum(ok)
  out <- list(
    N = N, alpha = alpha, reps = reps, converged = conv,
    mean_y = if (conv) mean(y) else NA_real_,
    bias = if (conv) mean(y) - true_mean else NA_real_,
    emp_sd = if (conv) stats::sd(y) else NA_real_,
    mean_se = if (conv) mean(se) else NA_real_,
    sd_se_ratio = if (conv) (stats::sd(y) / max(mean(se), 1e-8)) else NA_real_
  )
  class(out) <- "el_mc_summary"
  out
}

run_mc_grid <- function(Ns = c(400L, 800L, 2000L), alphas = c(0.3, 0.5), reps_small = 30L, reps_large = 20L, n_large_threshold = 1500L, seed = 1L) {
  res <- list()
  set.seed(seed)
  k <- 1L
  for (N in Ns) {
    for (a in alphas) {
      reps <- if (N <= n_large_threshold) reps_small else reps_large
      cat(sprintf("Running MC: N=%d, alpha=%.2f, reps=%d...\n", N, a, reps))
      res[[k]] <- run_mc(N, a, reps = reps, seed = sample.int(1e6, 1))
      k <- k + 1L
    }
  }
  res
}

as.data.frame.el_mc_summary <- function(x, ...) {
  data.frame(
    N = x$N, alpha = x$alpha, reps = x$reps, converged = x$converged,
    mean_y = x$mean_y, bias = x$bias, emp_sd = x$emp_sd, mean_se = x$mean_se, sd_se_ratio = x$sd_se_ratio
  )
}

bind_mc <- function(lst) do.call(rbind, lapply(lst, as.data.frame))

## Survey variant
simulate_design <- function(N, alpha, gamma = 0.7, seed = NULL) {
  if (!requireNamespace("survey", quietly = TRUE)) stop("Install 'survey' to run survey MC.")
  if (!is.null(seed)) set.seed(seed)
  X <- rnorm(N)
  Z <- rnorm(N)
  Y <- 2 + 0.5 * X + Z
  p <- plogis(-1.0 + alpha * scale(Y)[, 1])
  R <- runif(N) < p
  w <- exp(gamma * scale(X)[, 1]); w <- as.numeric(w / mean(w))
  df <- data.frame(Y_miss = Y, X = X, w = w)
  df[!R, "Y_miss"] <- NA_real_
  survey::svydesign(ids = ~1, weights = ~w, data = df)
}

fit_el_svy <- function(design, standardize = FALSE) {
  NMAR::nmar(
    formula = Y_miss ~ X,
    data = design,
    engine = NMAR::el_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = standardize)
  )
}

fit_el_svy_bootstrap <- function(design, standardize = FALSE, bootstrap_reps = 200L) {
  if (!requireNamespace("survey", quietly = TRUE)) stop("Install 'survey' to run survey MC.")
  NMAR::nmar(
    formula = Y_miss ~ X,
    data = design,
    engine = NMAR::el_engine(auxiliary_means = c(X = 0), variance_method = "bootstrap", bootstrap_reps = bootstrap_reps, standardize = standardize)
  )
}

run_mc_survey_once <- function(N, alpha, gamma, i) {
  des <- simulate_design(N, alpha, gamma)
  fit <- try(fit_el_svy(des), silent = TRUE)
  if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
    return(list(ok = FALSE))
  }
  list(ok = TRUE, y_hat = fit$estimate, se = fit$se)
}

run_mc_survey <- function(N, alpha, gamma = 0.7, reps = 20L, true_mean = 2.0, seed = 1L) {
  if (!requireNamespace("survey", quietly = TRUE)) stop("Install 'survey' to run survey MC.")
  set.seed(seed)
  vals <- replicate(reps, run_mc_survey_once(N, alpha, gamma, 0L), simplify = FALSE)
  ok <- vapply(vals, function(v) isTRUE(v$ok), logical(1))
  y <- vapply(vals[ok], function(v) v$y_hat, numeric(1))
  se <- vapply(vals[ok], function(v) v$se, numeric(1))
  conv <- sum(ok)
  out <- list(
    N = N, alpha = alpha, gamma = gamma, reps = reps, converged = conv,
    mean_y = if (conv) mean(y) else NA_real_,
    bias = if (conv) mean(y) - true_mean else NA_real_,
    emp_sd = if (conv) stats::sd(y) else NA_real_,
    mean_se = if (conv) mean(se) else NA_real_,
    sd_se_ratio = if (conv) (stats::sd(y) / max(mean(se), 1e-8)) else NA_real_
  )
  class(out) <- "el_mc_summary_svy"
  out
}

run_mc_survey_bootstrap <- function(N, alpha, gamma = 0.7, reps = 20L, true_mean = 2.0, seed = 1L, bootstrap_reps = 200L) {
  if (!requireNamespace("survey", quietly = TRUE)) stop("Install 'survey' to run survey MC.")
  set.seed(seed)
  vals <- replicate(reps, {
    des <- simulate_design(N, alpha, gamma)
    fit <- try(fit_el_svy_bootstrap(design = des, bootstrap_reps = bootstrap_reps), silent = TRUE)
    if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
      list(ok = FALSE)
    } else {
      list(ok = TRUE, y_hat = fit$estimate, se = fit$se)
    }
  }, simplify = FALSE)
  ok <- vapply(vals, function(v) isTRUE(v$ok), logical(1))
  y <- vapply(vals[ok], function(v) v$y_hat, numeric(1))
  se <- vapply(vals[ok], function(v) v$se, numeric(1))
  conv <- sum(ok)
  out <- list(
    N = N, alpha = alpha, gamma = gamma, reps = reps, converged = conv,
    mean_y = if (conv) mean(y) else NA_real_,
    bias = if (conv) mean(y) - true_mean else NA_real_,
    emp_sd = if (conv) stats::sd(y) else NA_real_,
    mean_se = if (conv) mean(se) else NA_real_,
    sd_se_ratio = if (conv) (stats::sd(y) / max(mean(se), 1e-8)) else NA_real_
  )
  class(out) <- "el_mc_summary_svy"
  out
}

as.data.frame.el_mc_summary_svy <- function(x, ...) {
  data.frame(
    N = x$N, alpha = x$alpha, gamma = x$gamma, reps = x$reps, converged = x$converged,
    mean_y = x$mean_y, bias = x$bias, emp_sd = x$emp_sd, mean_se = x$mean_se, sd_se_ratio = x$sd_se_ratio
  )
}

bind_mc_svy <- function(lst) do.call(rbind, lapply(lst, as.data.frame))

if (identical(environment(), globalenv())) {
# Expanded grid and simple checks (mirroring earlier unit-test heuristics)
  grid <- run_mc_grid(Ns = c(400L, 800L, 2000L), alphas = c(0.3, 0.4, 0.5), reps_small = 30L, reps_large = 20L)
  tab <- bind_mc(grid)

# Heuristic checks
  bias_tol <- 0.4
  ratio_tol <- 0.5
  conv_min_frac <- 0.6
  tab$conv_frac <- tab$converged / tab$reps
  tab$ratio_diff <- abs(tab$emp_sd - tab$mean_se) / pmax(tab$mean_se, 1e-8)
  tab$check_bias <- abs(tab$bias) < bias_tol
  tab$check_ratio <- tab$ratio_diff < ratio_tol
  tab$check_conv <- tab$conv_frac >= conv_min_frac

  print(tab)

  cat("\nAnalysis:\n")
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    cat(sprintf(
      "N=%d, alpha=%.2f: conv=%d/%d (%.0f%%), bias=%.3f [%s], emp_sd=%.4f, mean_se=%.4f, |sd-se|/se=%.3f [%s]\n",
      r$N, r$alpha, r$converged, r$reps, 100 * r$conv_frac,
      r$bias, if (r$check_bias) "OK" else "WARN",
      r$emp_sd, r$mean_se, r$ratio_diff, if (r$check_ratio) "OK" else "WARN"
    ))
  }
  cat("\nGuidelines: bias should shrink with N; sd≈se under delta method indicates calibration. For strong nonresponse (alpha high) or trimming, prefer bootstrap variance.\n")

# Survey MC
  if (requireNamespace("survey", quietly = TRUE)) {
    cat("\nRunning survey MC (delta variance, design-based B)...\n")
    grid_svy <- list(
      run_mc_survey(400, 0.4, gamma = 0.7, reps = 20L, seed = 2L),
      run_mc_survey(800, 0.4, gamma = 0.7, reps = 20L, seed = 3L),
      run_mc_survey(2000, 0.4, gamma = 0.7, reps = 20L, seed = 4L)
    )
    tab_svy <- bind_mc_svy(grid_svy)
    print(tab_svy)
    cat("\nSurvey analysis: compare emp_sd vs mean_se; closer ratios indicate better delta calibration under design-based covariance.\n")

# IID bootstrap MC
    cat("\nRunning IID MC (bootstrap variance) ...\n")
    grid_boot <- list(
      run_mc_bootstrap(400, 0.4, reps = 10L, seed = 21L, bootstrap_reps = 100L),
      run_mc_bootstrap(800, 0.4, reps = 10L, seed = 22L, bootstrap_reps = 100L)
    )
    tab_boot <- bind_mc(grid_boot)
    print(tab_boot)
    cat("\nIID bootstrap analysis: ratios closer to 1 than delta indicate improved calibration.\n")

# Survey bootstrap MC via svrep
    if (requireNamespace("svrep", quietly = TRUE)) {
      cat("\nRunning survey MC (bootstrap variance via svrep) ...\n")
      grid_svy_boot <- list(
        run_mc_survey_bootstrap(400, 0.4, gamma = 0.7, reps = 5L, seed = 31L, bootstrap_reps = 40L)
      )
      tab_svy_boot <- bind_mc_svy(grid_svy_boot)
      print(tab_svy_boot)
      cat("\nSurvey bootstrap analysis: expect sd≈se substantially closer to 1 than with delta.\n")
    } else {
      cat("\nSkipping survey bootstrap MC: 'svrep' package not installed.\n")
    }
  } else {
    cat("\nSkipping survey MC: 'survey' package not installed.\n")
  }
}
