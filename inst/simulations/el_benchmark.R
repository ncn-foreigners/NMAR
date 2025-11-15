## Empirical Likelihood benchmarking scaffold
##
## This script exercises the EL engine across a range of data-generating
## settings and nleqslv controls, measuring convergence, accuracy, and time.
##
## Usage:
##   Sys.setenv(RUN_EL_BENCH = "1")
##   Rscript inst/simulations/el_benchmark.R
##
## Outputs:
##   - Prints a compact list of scenario/control summaries
##   - Optionally, you can adapt `write.csv` calls to persist detailed rows

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    try(devtools::load_all("."), silent = TRUE)
  }
# Attach package if available (prefer the local, loaded version)
  if ("NMAR" %in% .packages(all.available = TRUE)) {
    suppressWarnings(try(library(NMAR), silent = TRUE))
  }
  suppressWarnings(requireNamespace("survey", quietly = TRUE))
  suppressWarnings(requireNamespace("MASS", quietly = TRUE))
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# DGPs
r_dgp_paper <- function(N, model = 1L, a = -2, resp_family = c("logit", "probit")) {
  resp_family <- match.arg(resp_family)
  X <- stats::rchisq(N, df = 2)
  eps <- stats::rnorm(N)
  Y <- switch(
    as.integer(model),
    X + eps * sqrt(X) / 5,
    X + 0.05 * X^2 + eps * sqrt(X) / 5,
    1.05 + X + eps * sqrt(X) / 5,
    3 + X - 0.05 * X^2 + eps * sqrt(X) / 5,
    5 + eps * sqrt(X) / 5
  )
  p <- if (resp_family == "logit") stats::plogis(0.2 * Y - a) else stats::pnorm(0.2 * Y - a)
  R <- stats::rbinom(N, 1, p)
  data.frame(y_miss = ifelse(R == 1, Y, NA_real_), x = X, R = R, y_true = Y)
}

r_dgp_linear_hd <- function(N, p = 30, rho = 0.7, sparse = TRUE, resp_family = c("logit", "probit")) {
  resp_family <- match.arg(resp_family)
  idx <- 0:(p - 1)
  Sigma <- outer(idx, idx, function(i, j) rho^abs(i - j))
  X <- MASS::mvrnorm(N, rep(0, p), Sigma)
  if (sparse) {
    beta <- numeric(p); beta[1:min(5, p)] <- c(0.8, 0.6, 0.4, 0.2, -0.3)[seq_len(min(5, p))]
  } else {
    beta <- seq(0.5, 0.5 / max(1, p), length.out = p)
  }
  Y <- drop(X %*% beta) + stats::rnorm(N, 0, 1)
  eta <- -0.2 + 0.3 * Y + 0.2 * X[, 1]
  pR <- if (resp_family == "logit") stats::plogis(eta) else stats::pnorm(eta)
  R <- stats::rbinom(N, 1, pR)
  df <- as.data.frame(X); names(df) <- paste0("x", seq_len(p))
  df$y_miss <- ifelse(R == 1, Y, NA_real_)
  df$R <- R; df$y_true <- Y
  df
}

make_design <- function(df, gamma = 0, strata = FALSE) {
  stopifnot(requireNamespace("survey", quietly = TRUE))
  x1 <- if ("x" %in% names(df)) df$x else df$x1
  w <- exp(gamma * scale(x1)[, 1])
  w <- as.numeric(w / mean(w))
  if (isTRUE(strata)) {
    st <- cut(x1, breaks = 4)
    survey::svydesign(ids = ~1, weights = ~w, strata = ~st, data = df)
  } else {
    survey::svydesign(ids = ~1, weights = ~w, data = df)
  }
}

# Build auxiliary means possibly perturbed
aux_means_with_perturb <- function(df, aux_names, perturb = 0) {
  if (is.null(aux_names) || length(aux_names) == 0) return(NULL)
  mu <- vapply(aux_names, function(nm) mean(df[[nm]]), numeric(1))
  if (isTRUE(perturb > 0)) {
    set.seed(1L)
    sgn <- sample(c(-1, 1), length(mu), replace = TRUE)
    mu <- mu * (1 + perturb * sgn)
  }
  mu
}

# Controls grid
control_grid <- function() {
  expand.grid(
    global = I(list(NULL, "dbldog", "qline")),
    xscalm = I(list(NULL, "auto", "fixed")),
    xtol = c(1e-8, 1e-10),
    ftol = c(1e-8, 1e-10),
    maxit = c(100L, 200L),
    btol = I(list(NULL, 1e-6)),
    stepmax = I(list(NULL, 50L)),
    delta = I(list(NULL, 1e-3)),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
}

# One fit runner
run_fit <- function(obj, frm, aux_means = NULL, control, family = c("logit", "probit"),
                    variance_method = c("delta", "bootstrap"), variance_jacobian = c("auto", "analytic", "numeric"),
                    trim_cap = Inf, standardize = TRUE, bootstrap_reps = 200,
                    variance_ridge = FALSE, variance_pseudoinverse = FALSE) {
  family <- match.arg(family)
  variance_method <- match.arg(variance_method)
  variance_jacobian <- match.arg(variance_jacobian)
  eng <- NMAR::el_engine(
    control = Filter(Negate(is.null), control),
    variance_method = variance_method,
    variance_jacobian = variance_jacobian,
    solver_jacobian = "analytic",
    auxiliary_means = aux_means,
    standardize = standardize,
    trim_cap = trim_cap,
    bootstrap_reps = bootstrap_reps,
    variance_ridge = variance_ridge,
    variance_pseudoinverse = variance_pseudoinverse,
    suppress_warnings = TRUE,
    family = family
  )
  t0 <- proc.time()
  fit <- try(suppressMessages(NMAR::nmar(frm, data = obj, engine = eng)), silent = TRUE)
  elapsed <- as.numeric((proc.time() - t0)[3])
  if (inherits(fit, "try-error")) {
    return(list(ok = FALSE, time = elapsed))
  }
  di <- fit$diagnostics %||% list()
  df_val <- tryCatch({
    inf <- fit$inference
    if (is.list(inf)) as.numeric(inf$df) else NA_real_
  }, error = function(e) NA_real_)
  se_val <- tryCatch({
    v <- stats::vcov(fit)
    sqrt(as.numeric(v[1, 1]))
  }, error = function(e) NA_real_)
  list(
    ok = isTRUE(fit$converged),
    yhat = as.numeric(fit$estimate),
    se = se_val,
    df = df_val,
    time = elapsed,
    iter = di$solver_iterations %||% NA_real_,
    method = di$solver_method %||% NA_character_,
    max_eq = di$max_equation_residual %||% NA_real_,
    kappa = di$jacobian_condition_number %||% NA_real_,
    rel_diff = di$jacobian_rel_diff %||% NA_real_,
    jac_src = di$jacobian_source %||% NA_character_,
    min_den = di$min_denominator %||% NA_real_,
    frac_small = di$fraction_small_denominators %||% NA_real_,
    trimmed = di$trimmed_fraction %||% NA_real_,
    invert_rule = di$invert_rule %||% NA_character_,
    used_pinv = isTRUE(di$used_pseudoinverse),
    used_ridge = isTRUE(di$used_ridge),
    solve_time = di$solver_time %||% NA_real_,
    var_time = di$variance_time %||% NA_real_
  )
}

# Scenario runner
run_scenario <- function(N = 1000, trials = 200, scenario, ctrl_row, survey = FALSE,
                         variance = list(method = "delta", jac = "auto"),
                         trim_caps = list(Inf), standardize = TRUE,
                         family_fit = "logit", family_dgp = "logit",
                         bootstrap_reps = 200, aux_perturb = 0,
                         var_ridge = FALSE, var_pinv = FALSE) {
  set.seed(1L)
  ctl <- list(
    global  = ctrl_row$global[[1]],
    xscalm  = ctrl_row$xscalm[[1]],
    xtol    = ctrl_row$xtol,
    ftol    = ctrl_row$ftol,
    maxit   = ctrl_row$maxit,
    btol    = ctrl_row$btol[[1]],
    stepmax = ctrl_row$stepmax[[1]],
    delta   = ctrl_row$delta[[1]]
  )
  res <- vector("list", trials)
  for (t in seq_len(trials)) {
    if (scenario$type == "paper") {
      df <- r_dgp_paper(N, model = scenario$model, a = scenario$a, resp_family = family_dgp)
      if (isTRUE(scenario$aux)) {
        frm <- y_miss ~ x
        aux_means <- aux_means_with_perturb(df, c("x"), perturb = aux_perturb)
      } else {
        frm <- y_miss ~ 1
        aux_means <- NULL
      }
    } else {
      df <- r_dgp_linear_hd(N, p = scenario$p, rho = scenario$rho, sparse = scenario$sparse, resp_family = family_dgp)
      if (!is.null(scenario$auxk) && scenario$auxk > 0) {
        aux_names <- paste0("x", seq_len(min(scenario$auxk, scenario$p)))
        frm <- as.formula(paste("y_miss ~", paste(aux_names, collapse = " + ")))
        aux_means <- aux_means_with_perturb(df, aux_names, perturb = aux_perturb)
      } else {
        frm <- as.formula(paste("y_miss ~", paste0("x", seq_len(scenario$p), collapse = " + ")))
        aux_means <- NULL
      }
    }
    truth <- mean(df$y_true)
    obj <- if (survey) make_design(df, gamma = scenario$gamma %||% 0, strata = scenario$strata %||% FALSE) else df
    for (cap in trim_caps) {
      out <- run_fit(
        obj, frm, aux_means, control = ctl, family = family_fit,
        variance_method = variance$method, variance_jacobian = variance$jac,
        trim_cap = if (is.character(cap) && cap == "Inf") Inf else as.numeric(cap),
        standardize = standardize, bootstrap_reps = bootstrap_reps,
        variance_ridge = var_ridge, variance_pseudoinverse = var_pinv
      )
      out$trim_cap <- cap
      out$truth <- truth
      res[[length(res) + 1L]] <- out
    }
  }
  res
}

# Aggregation
summarize_results <- function(res_list) {
  get_num <- function(x, nm) {
    val <- x[[nm]]
    if (is.null(val)) return(NA_real_)
    as.numeric(val)[1]
  }
  get_lgl <- function(x, nm) {
    val <- x[[nm]]
    if (is.null(val)) return(NA)
    as.logical(val)[1]
  }
  ok <- vapply(res_list, function(x) get_lgl(x, "ok"), logical(1))
  time <- vapply(res_list, function(x) get_num(x, "time"), numeric(1))
  iter <- vapply(res_list, function(x) get_num(x, "iter"), numeric(1))
  max_eq <- vapply(res_list, function(x) get_num(x, "max_eq"), numeric(1))
  kappa <- vapply(res_list, function(x) get_num(x, "kappa"), numeric(1))
  trimmed <- vapply(res_list, function(x) get_num(x, "trimmed"), numeric(1))
  yhat <- vapply(res_list, function(x) get_num(x, "yhat"), numeric(1))
  truth <- vapply(res_list, function(x) get_num(x, "truth"), numeric(1))
  se <- vapply(res_list, function(x) get_num(x, "se"), numeric(1))
  dfv <- vapply(res_list, function(x) get_num(x, "df"), numeric(1))
  solve_time <- vapply(res_list, function(x) get_num(x, "solve_time"), numeric(1))
  var_time <- vapply(res_list, function(x) get_num(x, "var_time"), numeric(1))
  used_ridge <- vapply(res_list, function(x) get_num(x, "used_ridge"), numeric(1))
  used_pinv <- vapply(res_list, function(x) get_num(x, "used_pinv"), numeric(1))
  p90 <- function(x) stats::quantile(x, 0.9, na.rm = TRUE, names = FALSE)
# Coverage: t-based if df available, else normal
  coverage <- rep(NA_real_, length(yhat))
  for (i in seq_along(yhat)) {
    if (!is.finite(se[i])) next
    crit <- if (is.finite(dfv[i])) stats::qt(0.975, df = dfv[i]) else stats::qnorm(0.975)
    lo <- yhat[i] - crit * se[i]
    hi <- yhat[i] + crit * se[i]
    coverage[i] <- as.numeric(truth[i] >= lo & truth[i] <= hi)
  }
  list(
    conv_rate = mean(ok, na.rm = TRUE),
    med_time = stats::median(time, na.rm = TRUE),
    p90_time = p90(time),
    med_solve_time = stats::median(solve_time, na.rm = TRUE),
    med_var_time = stats::median(var_time, na.rm = TRUE),
    med_iter = stats::median(iter, na.rm = TRUE),
    med_max_eq = stats::median(max_eq, na.rm = TRUE),
    med_kappa = stats::median(kappa, na.rm = TRUE),
    med_trimmed = stats::median(trimmed, na.rm = TRUE),
    se_ok_rate = mean(is.finite(se), na.rm = TRUE),
    ridge_rate = mean(used_ridge == 1, na.rm = TRUE),
    pinv_rate = mean(used_pinv == 1, na.rm = TRUE),
    coverage_rate = mean(coverage, na.rm = TRUE),
    bias = mean(yhat - truth, na.rm = TRUE),
    rmse = sqrt(mean((yhat - truth)^2, na.rm = TRUE))
  )
}

# Example quick sweep
if (identical(Sys.getenv("RUN_EL_BENCH"), "1")) {
  grid <- control_grid()
  n_pick <- as.integer(Sys.getenv("BENCH_CONTROLS", unset = "4"))
  pick <- grid[seq_len(max(1L, min(n_pick, nrow(grid)))), , drop = FALSE]

# Scenario set (can be expanded)
  scenarios <- list(
    list(type = "paper", model = 1L, a = -2, aux = TRUE),
    list(type = "paper", model = 5L, a = -1, aux = FALSE),
    list(type = "hd", p = 30, rho = 0.7, sparse = TRUE),
    list(type = "hd", p = 60, rho = 0.7, sparse = FALSE)
  )
  trials <- as.integer(Sys.getenv("BENCH_TRIALS", unset = "100"))
  fam <- Sys.getenv("BENCH_FAMILY", unset = "logit")
  var_jac <- Sys.getenv("BENCH_VARJAC", unset = "auto")
  var_modes <- strsplit(Sys.getenv("BENCH_VAR_MODES", unset = "delta"), ",")[[1]]
  survey_flag <- as.logical(as.integer(Sys.getenv("BENCH_SURVEY", unset = "0")))
  trim_caps_env <- Sys.getenv("BENCH_TRIM_CAPS", unset = "Inf")
  trim_caps <- strsplit(trim_caps_env, ",")[[1]]
  aux_perturb <- as.numeric(Sys.getenv("BENCH_AUX_PERTURB", unset = "0"))
  family_dgp <- Sys.getenv("BENCH_RESP_DGP", unset = fam)
  family_fit <- Sys.getenv("BENCH_RESP_FIT", unset = fam)
  boot_reps <- as.integer(Sys.getenv("BENCH_BOOT_REPS", unset = "200"))
  var_ridge_vals <- strsplit(Sys.getenv("BENCH_VAR_RIDGE", unset = "FALSE"), ",")[[1]]
  var_pinv_vals <- strsplit(Sys.getenv("BENCH_VAR_PINV", unset = "FALSE"), ",")[[1]]
  out_path <- Sys.getenv("BENCH_OUTPUT", unset = "")

  results <- list(); idx <- 1L
  rows <- list(); rix <- 1L
  for (i in seq_len(nrow(pick))) {
    for (sc in scenarios) {
      for (vm in var_modes) for (vr in var_ridge_vals) for (vp in var_pinv_vals) {
        cat(sprintf("Running ctrl %d scenario %s (trials=%d, var=%s/%s, ridge=%s, pinv=%s, survey=%s, trims=%s) ...\n",
            i, sc$type, trials, vm, var_jac, vr, vp, survey_flag, paste(trim_caps, collapse = "/")))
        res <- run_scenario(N = 1000, trials = trials, scenario = sc, ctrl_row = pick[i, ], survey = survey_flag,
                            variance = list(method = vm, jac = var_jac),
                            trim_caps = as.list(trim_caps), standardize = TRUE,
                            family_fit = family_fit, family_dgp = family_dgp,
                            bootstrap_reps = boot_reps, aux_perturb = aux_perturb,
                            var_ridge = as.logical(as.integer(vr)), var_pinv = as.logical(as.integer(vp)))
        summ <- summarize_results(res)
# Flatten to a data.frame row
        to_str <- function(x) if (is.null(x)) "NULL" else as.character(x)
        row <- data.frame(
          global = to_str(pick$global[[i]]), xscalm = to_str(pick$xscalm[[i]]),
          xtol = pick$xtol[i], ftol = pick$ftol[i], maxit = pick$maxit[i],
          btol = to_str(pick$btol[[i]]), stepmax = to_str(pick$stepmax[[i]]), delta = to_str(pick$delta[[i]]),
          scenario_type = sc$type, model = sc$model %||% NA_integer_, p = sc$p %||% NA_integer_,
          rho = sc$rho %||% NA_real_, sparse = sc$sparse %||% NA, aux = sc$aux %||% NA, auxk = sc$auxk %||% NA,
          family_fit = family_fit, family_dgp = family_dgp, variance_method = vm, variance_jacobian = var_jac,
          survey = survey_flag, trials = trials, trim_caps = paste(trim_caps, collapse = ","),
          variance_ridge = vr, variance_pinv = vp,
          conv_rate = summ$conv_rate, med_time = summ$med_time, p90_time = summ$p90_time,
          med_solve_time = summ$med_solve_time, med_var_time = summ$med_var_time,
          med_iter = summ$med_iter, med_max_eq = summ$med_max_eq, med_kappa = summ$med_kappa,
          med_trimmed = summ$med_trimmed, se_ok_rate = summ$se_ok_rate,
          ridge_rate = summ$ridge_rate, pinv_rate = summ$pinv_rate,
          coverage_rate = summ$coverage_rate,
          bias = summ$bias, rmse = summ$rmse,
          stringsAsFactors = FALSE, check.names = FALSE
        )
        rows[[rix]] <- row; rix <- rix + 1L
        idx <- idx + 1L
      }
    }
  }
  df <- do.call(rbind, rows)
  print(df)
  if (nzchar(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(df, out_path, row.names = FALSE)
  }
}
