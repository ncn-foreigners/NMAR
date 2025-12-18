#' EL core helpers
#' @description Internal helpers for solving and post-processing the EL system.
#'   \code{el_run_solver()} orchestrates \code{nleqslv::nleqslv()} with a small,
#'   deterministic fallback ladder; \code{el_post_solution()} computes masses and
#'   the point estimate with
#'   denominator guards and optional trimming.
#' @name el_core_helpers
#' @keywords internal
NULL



#' Solver orchestration with staged policy
#'
#' @param equation_system_func Function mapping parameter vector to equation
#'   residuals.
#' @param analytical_jac_func Analytic Jacobian function; may be NULL if
#'   unavailable or when forcing Broyden.
#' @param init Numeric vector of initial parameter values.
#' @param final_control List passed to \code{nleqslv::nleqslv(control = ...)}.
#' @param top_args List of top-level \code{nleqslv::nleqslv} args (e.g., \code{global}, \code{xscalm}).
#' @param solver_method Character; one of "auto", "newton", or "broyden".
#' @param use_solver_jac Logical; whether to pass analytic Jacobian to Newton.
#' @param K_beta Integer; number of response model parameters.
#' @param K_aux Integer; number of auxiliary constraints.
#' @param respondent_weights Numeric vector of base sampling weights.
#' @param N_pop Numeric; population total (weighted when survey design).
#' @param trace_level Integer; verbosity level (0 silent, 1-3 increasingly verbose).
#'
#' @keywords internal
el_run_solver <- function(equation_system_func,
                          analytical_jac_func,
                          init,
                          final_control,
                          top_args,
                          solver_method,
                          use_solver_jac,
                          K_beta,
                          K_aux,
                          respondent_weights,
                          N_pop,
                          trace_level = 0) {
  verboser <- create_verboser(trace_level)
  solver_method_used <- "Newton"
# Curated defaults: Newton + analytic Jacobian, quadratic line search, auto scaling
  user_specified_global <- !is.null(top_args$global)
  user_specified_xscalm <- !is.null(top_args$xscalm)

  nl_args <- list(
    x = init,
    fn = equation_system_func,
    jac = if (use_solver_jac) analytical_jac_func else NULL,
    method = "Newton",
    control = final_control
  )
# Defaults if user did not specify: prefer robust line search and auto scaling
  if (!is.null(top_args$global)) nl_args$global <- top_args$global else nl_args$global <- "qline"
  if (!is.null(top_args$xscalm)) nl_args$xscalm <- top_args$xscalm else nl_args$xscalm <- "auto"

  nl_args_b <- nl_args
  nl_args_b$x <- init
  nl_args_b$jac <- NULL
  nl_args_b$method <- "Broyden"

  global_used <- nl_args$global %||% NULL
  xscalm_used <- nl_args$xscalm %||% NULL

  if (solver_method == "broyden") {
    broyden_control <- final_control
    if (!is.null(broyden_control$maxit) && is.finite(broyden_control$maxit) && broyden_control$maxit < 5) {
      broyden_control$maxit <- 50
    }
    nl_args_b$control <- broyden_control
    solution <- do.call(nleqslv::nleqslv, nl_args_b)
    solver_method_used <- "Broyden"
    global_used <- nl_args_b$global %||% NULL
    xscalm_used <- nl_args_b$xscalm %||% NULL
    return(list(solution = solution, method = solver_method_used, used_top = list(global = global_used, xscalm = xscalm_used)))
  }

  solution <- do.call(nleqslv::nleqslv, nl_args)

# Small pool of stochastic perturbed restarts, used only on solver failure.
# For reproducibility across runs, set a seed before calling nmar().
  if (any(is.na(solution$x)) || solution$termcd > 2) {
    verboser("  Initial attempt failed, trying perturbed starts...", level = 3)
    for (i in seq_len(3)) {
      verboser(sprintf("    Restart attempt %d/3", i), level = 3)
      init_beta_perturbed <- init[seq_len(K_beta)] + stats::rnorm(K_beta, mean = 0, sd = 0.5)
      W_seed <- sum(respondent_weights) / N_pop
      W_seed <- min(max(W_seed, 1e-12), 1 - 1e-12)
      z_seed <- stats::qlogis(W_seed)
# Preserve generic tail structure: parameters after (beta, z) are reset to 0
      tail_len <- max(length(init) - K_beta - 1L, 0L)
      tail_zeros <- if (tail_len > 0L) rep(0, tail_len) else numeric(0)
      init_perturbed <- c(init_beta_perturbed, z_seed, tail_zeros)
      nl_args$x <- init_perturbed
      solution <- do.call(nleqslv::nleqslv, nl_args)
      if (!any(is.na(solution$x)) && solution$termcd <= 2) {
        verboser(sprintf("    [OK] Restart %d converged successfully", i), level = 3, type = "result")
        break
      }
    }
  }

# Minimal, deterministic fallback ladder
  if (solver_method == "auto" && (any(is.na(solution$x)) || solution$termcd > 2)) {
# If the user did not choose a global strategy, try an alternative globalization before switching method
    if (!user_specified_global) {
      alt_args <- nl_args
      alt_args$global <- if (identical(nl_args$global, "qline")) "dbldog" else "qline"
      verboser(sprintf("  Newton failed; retry with global='%s'...", alt_args$global), level = 2)
      solution2 <- do.call(nleqslv::nleqslv, alt_args)
      if (!any(is.na(solution2$x)) && solution2$termcd <= 2) {
        solution <- solution2
        global_used <- alt_args$global %||% NULL
        xscalm_used <- alt_args$xscalm %||% NULL
      }
    }
  }

  if (solver_method == "auto" && (any(is.na(solution$x)) || solution$termcd > 2)) {
    verboser("  Newton method failed, falling back to Broyden...", level = 2)
    broyden_control <- final_control
    if (!is.null(broyden_control$maxit) && is.finite(broyden_control$maxit) && broyden_control$maxit < 5) {
      broyden_control$maxit <- 50
    }
    nl_args_b$control <- broyden_control
    solution <- do.call(nleqslv::nleqslv, nl_args_b)
    solver_method_used <- "Broyden"
    global_used <- nl_args_b$global %||% NULL
    xscalm_used <- nl_args_b$xscalm %||% NULL
  }

  list(solution = solution, method = solver_method_used, used_top = list(global = global_used, xscalm = xscalm_used))
}



el_post_solution <- function(estimates,
                             missingness_model_matrix_scaled,
                             missingness_model_matrix_unscaled,
                             auxiliary_matrix_scaled,
                             mu_x_scaled,
                             response_outcome,
                             family,
                             N_pop,
                             respondent_weights,
                             K_beta,
                             K_aux,
                             nmar_scaling_recipe,
                             standardize,
                             trim_cap,
                             X_centered = NULL,
                             lambda_W = NULL) {
  beta_hat_scaled <- estimates[1:K_beta]
  names(beta_hat_scaled) <- colnames(missingness_model_matrix_scaled)
  W_hat <- stats::plogis(estimates[K_beta + 1])
  W_hat <- min(max(W_hat, 1e-12), 1 - 1e-12)
# For survey designs the parameter vector includes lambda_W explicitly after z;
# for IID designs lambda_W is derived from (N_pop, sum(weights)).
  if (!is.null(lambda_W)) {
    lambda_hat <- if (K_aux > 0) estimates[(K_beta + 3):length(estimates)] else numeric(0)
  } else {
    lambda_hat <- if (K_aux > 0) estimates[(K_beta + 2):length(estimates)] else numeric(0)
  }
  eta_i_hat <- as.vector(missingness_model_matrix_scaled %*% beta_hat_scaled)
# Clip eta consistently with equations/Jacobian for diagnostic coherence
  ETA_CAP <- get_eta_cap()
  eta_i_hat_capped <- pmax(pmin(eta_i_hat, ETA_CAP), -ETA_CAP)
  w_i_hat <- family$linkinv(eta_i_hat_capped)
  if (is.null(lambda_W)) {
    C_const <- (N_pop / sum(respondent_weights)) - 1
    lambda_W_hat <- el_lambda_W(C_const, W_hat)
  } else {
    lambda_W_hat <- lambda_W
  }
  if (K_aux > 0) {
    if (is.null(X_centered)) {
# Fallback centering (should be provided by caller for efficiency)
      X_centered <- sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-")
    }
    Xc_lambda <- as.vector(X_centered %*% lambda_hat)
  } else {
    Xc_lambda <- 0
  }
  denom_floor <- nmar_get_el_denom_floor()
  dpack <- el_denominator(lambda_W_hat, W_hat, Xc_lambda, w_i_hat, denom_floor)
  masses <- el_masses(respondent_weights, dpack$denom, denom_floor, trim_cap)
  prob_mass <- masses$prob_mass
  y_hat <- el_mean(prob_mass, response_outcome)
  beta_hat_unscaled <- if (standardize) {
    unscale_coefficients(beta_hat_scaled, matrix(0, K_beta, K_beta), nmar_scaling_recipe)$coefficients
  } else {
    beta_hat_scaled
  }
  names(beta_hat_unscaled) <- colnames(missingness_model_matrix_unscaled)
  list(
    error = FALSE,
    y_hat = y_hat,
    weights = masses$mass_trimmed, # Store unnormalized (trimmed) EL masses as single source of truth
    mass_untrim = masses$mass_untrim,
    trimmed_fraction = masses$trimmed_fraction,
    beta_hat_scaled = beta_hat_scaled,
    beta_hat_unscaled = beta_hat_unscaled,
    W_hat = W_hat,
    lambda_hat = lambda_hat,
    eta_i_hat = eta_i_hat,
    w_i_hat = w_i_hat,
    denominator_hat = dpack$denom,
    lambda_W_hat = lambda_W_hat
  )
}

#' Build starting values for the EL solver (beta, z, lambda)
#' @keywords internal
el_build_start <- function(missingness_model_matrix_scaled,
                           auxiliary_matrix_scaled,
                           nmar_scaling_recipe,
                           start,
                           N_pop,
                           respondent_weights) {
  K_beta <- ncol(missingness_model_matrix_scaled)
  K_aux <- if (is.null(auxiliary_matrix_scaled)) 0L else ncol(auxiliary_matrix_scaled)

  init_beta <- rep(0, K_beta)
  names(init_beta) <- colnames(missingness_model_matrix_scaled)
  init_lambda <- if (K_aux > 0) setNames(rep(0, K_aux), colnames(auxiliary_matrix_scaled)) else numeric(0)
  W0 <- sum(respondent_weights) / N_pop
  W0 <- min(max(W0, 1e-12), 1 - 1e-12)
  init_z <- stats::qlogis(W0)

  if (!is.null(start) && is.list(start)) {
    if (!is.null(start$beta)) {
      init_beta <- scale_coefficients(start$beta, nmar_scaling_recipe, colnames(missingness_model_matrix_scaled))
    }
    if (!is.null(start$z)) {
      z_try <- as.numeric(start$z)[1]
      if (is.finite(z_try)) init_z <- z_try
    } else if (!is.null(start$W)) {
      W_try <- min(max(as.numeric(start$W)[1], 1e-12), 1 - 1e-12)
      init_z <- stats::qlogis(W_try)
    }
    if (K_aux > 0 && !is.null(start$lambda)) {
      init_lambda <- scale_aux_multipliers(start$lambda, nmar_scaling_recipe, colnames(auxiliary_matrix_scaled))
    }
  }

  init <- c(unname(init_beta), init_z, unname(init_lambda))
  list(init = init, init_beta = init_beta, init_z = init_z, init_lambda = init_lambda)
}

#' Prepare nleqslv top-level args and control
#' @keywords internal
el_prepare_nleqslv <- function(control) {
  control_top <- validate_nleqslv_top(extract_nleqslv_top(control))
  final_control <- modifyList(list(ftol = 1e-8, xtol = 1e-8, maxit = 200, trace = FALSE), control)
  final_control <- sanitize_nleqslv_control(final_control)
  list(top = control_top, control = final_control)
}

#' Compute diagnostics at the EL solution
#' @keywords internal
el_compute_diagnostics <- function(estimates,
                                   equation_system_func,
                                   analytical_jac_func,
                                   post,
                                   respondent_weights,
                                   auxiliary_matrix_scaled,
                                   K_beta,
                                   K_aux,
                                   X_centered) {
# Equation residuals and Jacobian conditioning
  eq_residuals <- tryCatch(equation_system_func(estimates), error = function(e) rep(NA_real_, length(estimates)))
  max_eq_resid <- suppressWarnings(max(abs(eq_residuals), na.rm = TRUE))
  A_condition <- tryCatch({ kappa(analytical_jac_func(estimates)) }, error = function(e) NA_real_)
  if (!is.finite(A_condition)) {
    A_condition <- NA_real_
  }

  is_survey_system <- length(estimates) == (K_beta + 2L + K_aux)
  constraint_sum_link <- if (is_survey_system && length(eq_residuals) == length(estimates)) {
    as.numeric(eq_residuals[length(eq_residuals)])
  } else {
    NA_real_
  }

  denom_floor <- nmar_get_el_denom_floor()
  denom_hat <- post$denominator_hat
  denom_stats <- list(
    min = suppressWarnings(min(denom_hat, na.rm = TRUE)),
    p_small = mean(denom_hat < 1e-6),
    p_floor = mean(denom_hat <= denom_floor)
  )

# Constraint summaries using untrimmed masses (m_i = d_i / D_i)
  cons <- tryCatch(
    constraint_summaries(post$w_i_hat, post$W_hat, post$mass_untrim, X_centered),
    error = function(e) {
      out <- list(constraint_sum_W = NA_real_, constraint_sum_aux = numeric(0))
      idx_W <- K_beta + 1L
      if (!is.null(eq_residuals) && length(eq_residuals) >= idx_W) {
        out$constraint_sum_W <- as.numeric(eq_residuals[idx_W])
        if (K_aux > 0 && length(eq_residuals) >= (idx_W + K_aux)) {
          aux_vals <- as.numeric(eq_residuals[(idx_W + 1):(idx_W + K_aux)])
          names(aux_vals) <- colnames(auxiliary_matrix_scaled)
          out$constraint_sum_aux <- aux_vals
        }
      }
      out
    }
  )

  sum_respondent_weights <- sum(respondent_weights)
  sum_unnormalized_weights_untrimmed <- sum(post$mass_untrim)
  normalization_ratio <- sum_unnormalized_weights_untrimmed / sum_respondent_weights

  denom_q <- tryCatch(stats::quantile(denom_hat, probs = c(0.01, 0.05, 0.5), na.rm = TRUE),
                      error = function(e) c(`1%` = NA_real_, `5%` = NA_real_, `50%` = NA_real_))
  denom_cnt_1e4 <- sum(denom_hat < 1e-4)
  weight_sum <- sum(post$weights)
  weight_share <- if (weight_sum > 0) sort(post$weights / weight_sum, decreasing = TRUE) else rep(NA_real_, length(post$weights))
  weight_max_share <- if (length(weight_share)) weight_share[1] else NA_real_
  weight_top5_share <- if (length(weight_share) >= 5) sum(weight_share[1:5]) else sum(weight_share, na.rm = TRUE)
  weight_ess <- if (weight_sum > 0) (weight_sum^2) / sum(post$weights^2) else NA_real_

  list(
    max_equation_residual = max_eq_resid,
    jacobian_condition_number = A_condition,
    denom_stats = denom_stats,
    constraint_sum_W = cons$constraint_sum_W,
    constraint_sum_aux = cons$constraint_sum_aux,
    constraint_sum_link = constraint_sum_link,
    normalization_ratio = normalization_ratio,
    sum_respondent_weights = sum_respondent_weights,
    sum_unnormalized_weights_untrimmed = sum_unnormalized_weights_untrimmed,
    denom_q = denom_q,
    denom_cnt_1e4 = denom_cnt_1e4,
    weight_max_share = weight_max_share,
    weight_top5_share = weight_top5_share,
    weight_ess = weight_ess,
    denom_floor = denom_floor
  )
}

#' Variance driver for EL (bootstrap or none)
#' @keywords internal
el_compute_variance <- function(y_hat,
                                full_data,
                                formula,
                                N_pop,
                                variance_method,
                                bootstrap_reps,
                                standardize,
                                trim_cap,
                                on_failure,
                                auxiliary_means,
                                control,
                                start,
                                family) {
  if (identical(variance_method, "none")) {
    return(list(se = NA_real_, message = "Variance skipped (variance_method='none')"))
  }
  if (identical(variance_method, "bootstrap")) {
    engine_args <- list(
      standardize = standardize,
      trim_cap = trim_cap,
      on_failure = on_failure,
      auxiliary_means = auxiliary_means,
      control = control,
      n_total = N_pop %||% NULL,
      start = start,
      family = family
    )
# Construct a local estimator closure to avoid cross-namespace calls
    est_closure <- function(data, formula, engine_args, ...) {
      engine_args$variance_method <- "none"
      eng <- do.call(el_engine, engine_args)
      nmar(formula = formula, data = data, engine = eng)
    }
    boot_args <- list(
      data = full_data,
      estimator_func = est_closure,
      point_estimate = y_hat,
      bootstrap_reps = bootstrap_reps,
      formula = formula,
      engine_args = engine_args
    )
    boot_try <- tryCatch(
      list(result = do.call(bootstrap_variance, boot_args), message = "Calculation successful"),
      error = function(e) list(result = NULL, message = paste("Bootstrap failed:", e$message))
    )
    if (!is.null(boot_try$result)) {
      return(list(se = boot_try$result$se, message = boot_try$message))
    } else {
      return(list(se = NA_real_, message = boot_try$message))
    }
  }
# Unknown method -> treat as none (defensive)
  list(se = NA_real_, message = sprintf("Unknown variance method '%s'", as.character(variance_method)))
}
