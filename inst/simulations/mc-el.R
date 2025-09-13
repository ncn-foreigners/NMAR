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
  nmar::nmar(
    formula = list(outcome = ~Y_miss, covariates_outcome = ~X, covariates_missingness = ~NULL),
    data = df,
    engine = el_engine(auxiliary_means = c(X = 0), variance_method = "delta", standardize = standardize)
  )
}

run_mc_once <- function(N, alpha, i) {
  df <- simulate_df(N, alpha)
  fit <- try(fit_el_df(df), silent = TRUE)
  if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
    return(list(ok = FALSE))
  }
  list(ok = TRUE, y_hat = fit$y_hat, se = fit$se)
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
}
