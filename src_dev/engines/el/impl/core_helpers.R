#' EL core helpers
#' @description Internal helpers for solving and post-processing the EL system.
#'   `el_run_solver()` orchestrates `nleqslv` with a small, deterministic fallback
#'   ladder; `el_post_solution()` computes masses and the point estimate with
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
#' @param final_control List passed to `nleqslv(control = ...)`.
#' @param top_args List of top-level `nleqslv` args (e.g., `global`, `xscalm`).
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
# Create verboser for this function
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

  if (solver_method == "broyden") {
    broyden_control <- final_control
    if (!is.null(broyden_control$maxit) && is.finite(broyden_control$maxit) && broyden_control$maxit < 5) {
      broyden_control$maxit <- 50
    }
    nl_args_b$control <- broyden_control
    solution <- do.call(nleqslv::nleqslv, nl_args_b)
    solver_method_used <- "Broyden"
    return(list(solution = solution, method = solver_method_used, used_top = list(global = top_args$global %||% NULL, xscalm = top_args$xscalm %||% NULL)))
  }

  solution <- do.call(nleqslv::nleqslv, nl_args)

# Small pool of perturbed restarts with the same configuration
  if (any(is.na(solution$x)) || solution$termcd > 2) {
    verboser("  Initial attempt failed, trying perturbed starts...", level = 3)
    for (i in seq_len(3)) {
      verboser(sprintf("    Restart attempt %d/3", i), level = 3)
      init_beta_perturbed <- init[seq_len(K_beta)] + stats::rnorm(K_beta, mean = 0, sd = 0.5)
      W_seed <- sum(respondent_weights) / N_pop
      W_seed <- min(max(W_seed, 1e-12), 1 - 1e-12)
      z_seed <- stats::qlogis(W_seed)
      init_perturbed <- c(init_beta_perturbed, z_seed, rep(0, K_aux))
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
      }
    }
  }

  if (solver_method == "auto" && (any(is.na(solution$x)) || solution$termcd > 2)) {
    verboser("  Newton method failed, falling back to Broyden...", level = 2)
    broyden_control <- final_control
    if (!is.null(broyden_control$maxit) && is.finite(broyden_control$maxit) && broyden_control$maxit < 5) {
      broyden_control$maxit <- 50
    }
# Prefer the same or more conservative global for Broyden
    nl_args_b$global <- nl_args$global %||% "qline"
    nl_args_b$xscalm <- nl_args$xscalm %||% "auto"
    nl_args_b$control <- broyden_control
    solution <- do.call(nleqslv::nleqslv, nl_args_b)
    solver_method_used <- "Broyden"
    if (!any(is.na(solution$x)) && solution$termcd <= 2) {
      verboser("  [OK] Broyden method converged successfully", level = 2, type = "result")
    }
  }
  list(solution = solution, method = solver_method_used, used_top = list(global = nl_args$global %||% NULL, xscalm = nl_args$xscalm %||% NULL))
}

#' Post-solution: denominators, masses, and mean (EL)
#'
#' Computes EL denominators and masses and the mean estimate after the solver
#' converges. The mean is calculated from normalized probability masses per
#' Qin, Leung, and Shao (2002, Eq. 11). When trimming is active, the masses are
#' capped and renormalized before computing the mean; diagnostics still report
#' constraint sums based on untrimmed masses. Denominator floors and the
#' guarding policy used in equations/Jacobian are applied consistently here for
#' diagnostic coherence.
#'
#' @param estimates Numeric vector (beta, z, lambda) at the solution.
#' @param response_model_matrix_scaled Scaled design matrix for the response model.
#' @param response_model_matrix_unscaled Unscaled design matrix for the response model.
#' @param auxiliary_matrix_scaled Scaled auxiliary matrix (or empty matrix).
#' @param mu_x_scaled Vector of population means for scaled auxiliaries (or NULL).
#' @param respondent_data Data frame of respondents.
#' @param outcome_var Character; outcome column name in respondent_data.
#' @param family Family object with linkinv and mu.eta.
#' @param N_pop Numeric; population total on the analysis scale.
#' @param respondent_weights Base weights for respondents.
#' @param K_beta,K_aux Integers; sizes of beta and lambda.
#' @param nmar_scaling_recipe Scaling recipe object for unscaling.
#' @param standardize Logical; whether coefficients need unscaling.
#' @param trim_cap Numeric; weight trimming cap (Inf = no trimming).
#'
#' @keywords internal
el_post_solution <- function(estimates,
                             response_model_matrix_scaled,
                             response_model_matrix_unscaled,
                             auxiliary_matrix_scaled,
                             mu_x_scaled,
                             respondent_data,
                             outcome_var,
                             family,
                             N_pop,
                             respondent_weights,
                             K_beta,
                             K_aux,
                             nmar_scaling_recipe,
                             standardize,
                             trim_cap,
                             X_centered = NULL) {
  beta_hat_scaled <- estimates[1:K_beta]
  names(beta_hat_scaled) <- colnames(response_model_matrix_scaled)
  W_hat <- stats::plogis(estimates[K_beta + 1])
  lambda_hat <- if (K_aux > 0) estimates[(K_beta + 2):length(estimates)] else numeric(0)
  eta_i_hat <- as.vector(response_model_matrix_scaled %*% beta_hat_scaled)
# Clip eta consistently with equations/Jacobian for diagnostic coherence
  ETA_CAP <- get_eta_cap()
  eta_i_hat_capped <- pmax(pmin(eta_i_hat, ETA_CAP), -ETA_CAP)
  w_i_hat <- family$linkinv(eta_i_hat_capped)
  C_const <- (N_pop / sum(respondent_weights)) - 1
  lambda_W_hat <- el_lambda_W(C_const, W_hat)
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
  y_hat <- el_mean(prob_mass, respondent_data[[outcome_var]])
  beta_hat_unscaled <- if (standardize) unscale_coefficients(beta_hat_scaled, matrix(0, K_beta, K_beta), nmar_scaling_recipe)$coefficients else beta_hat_scaled
  names(beta_hat_unscaled) <- colnames(response_model_matrix_unscaled)
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
