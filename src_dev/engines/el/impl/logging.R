#' Log a step banner line
#' @keywords internal
el_log_banner <- function(verboser, title) {
  verboser("============================================================", level = 1, type = "step")
  verboser(paste0("  ", title), level = 1, type = "step")
  verboser("============================================================", level = 1, type = "step")
}

#' Log a short trace message with the chosen level
#' @keywords internal
el_log_trace <- function(verboser, trace_level) {
  msg <- sprintf("Running with trace_level = %d", trace_level)
  if (trace_level < 3) {
    msg <- paste0(msg, sprintf(" | For more detail, use trace_level = %d", trace_level + 1))
  }
  verboser(msg, level = 1)
}

#' Log data prep summary
#' @keywords internal
el_log_data_prep <- function(verboser, outcome_var, family_name,
                             K_beta, K_aux, aux_names,
                             standardize, is_survey,
                             N_pop, n_resp_weighted) {
  verboser("", level = 1)
  verboser("-- DATA PREPARATION --", level = 1)
  response_rate <- (n_resp_weighted / N_pop) * 100
  verboser(sprintf("  Total weighted size:      %.1f", N_pop), level = 1)
  verboser(sprintf("  Respondents (weighted):   %.1f (%.1f%%)", n_resp_weighted, response_rate), level = 1)

  verboser("", level = 2)
  verboser("-- MODEL SPECIFICATION --", level = 2)
  verboser(sprintf("  Outcome variable:         %s", outcome_var), level = 2)
  verboser(sprintf("  Response model family:    %s", family_name %||% "<unknown>"), level = 2)
  verboser(sprintf("  Response predictors:      %d", K_beta), level = 2)
  if (K_aux > 0) {
    verboser(sprintf("  Auxiliary constraints:    %d", K_aux), level = 2)
    if (length(aux_names) > 0) verboser(sprintf("  Auxiliary variables:      %s", paste(aux_names, collapse = ", ")), level = 2)
  } else {
    verboser("  Auxiliary constraints:    (none)", level = 2)
  }
  verboser(sprintf("  Standardization:          %s", if (standardize) "enabled" else "disabled"), level = 2)
  verboser(sprintf("  Data type:                %s", if (is_survey) "survey.design" else "data.frame"), level = 2)
}

#' Log solver configuration
#' @keywords internal
el_log_solver_config <- function(verboser, control_top, final_control) {
  verboser("", level = 1)
  verboser("-- NONLINEAR SOLVER --", level = 1)
  verboser("  Method:                   Newton with analytic Jacobian", level = 2)
  verboser(sprintf("  Global strategy:          %s", control_top$global %||% "qline"), level = 2)
  verboser(sprintf("  Max iterations:           %d", final_control$maxit), level = 2)
  verboser(sprintf("  Function tolerance:       %.2e", final_control$ftol), level = 2)
  verboser(sprintf("  Parameter tolerance:      %.2e", final_control$xtol), level = 2)
}

#' Log starting values
#' @keywords internal
el_log_start_values <- function(verboser, init_beta, init_z, init_lambda) {
  verboser("", level = 3)
  verboser("  Starting values:", level = 3)
  if (length(init_beta)) verboser(sprintf("    beta (response model):  %s", paste(sprintf("%.4f", init_beta), collapse = ", ")), level = 3)
  verboser(sprintf("    z (logit response rate): %.4f", init_z), level = 3)
  if (length(init_lambda)) verboser(sprintf("    lambda_x (auxiliary):   %s", paste(sprintf("%.4f", init_lambda), collapse = ", ")), level = 3)
}

#' Log a short solver progress note
#' @keywords internal
el_log_solving <- function(verboser) {
  verboser("", level = 1)
  verboser("Solving stacked system...", level = 1)
}

#' Log solver termination status
#' @keywords internal
el_log_solver_result <- function(verboser, converged_success, solution, elapsed) {
  verboser("", level = 1)
  verboser(if (converged_success) "[OK] Solver converged successfully" else "[FAILED] Solver failed to converge", level = 1, type = "result")
  verboser("", level = 2)
  verboser(sprintf("  Termination code:         %d (%s)", solution$termcd, solution$message), level = 2)
  verboser(sprintf("  Iterations:               %d", if (!is.null(solution$iter)) solution$iter else NA), level = 2)
  verboser(sprintf("  Solver time:              %.3f seconds", elapsed), level = 2)
}

#' Log weight diagnostics
#' @keywords internal
el_log_weight_diagnostics <- function(verboser, W_hat, weights, trimmed_fraction) {
  verboser("", level = 2)
  verboser("-- EL MASS DIAGNOSTICS --", level = 2)
  weight_sum <- sum(weights)
  verboser(sprintf("  Estimated response rate:  %.4f", W_hat), level = 2)
  verboser(sprintf("  Weight sum (trimmed):     %.1f", weight_sum), level = 2)
  verboser(sprintf("  Trimmed fraction:         %.2f%%", trimmed_fraction * 100), level = 2)
  wr <- range(weights)
  verboser(sprintf("  Weight range:             [%.4f, %.4f]", wr[1], wr[2]), level = 2)
}

#' Log detailed diagnostics
#' @keywords internal
el_log_detailed_diagnostics <- function(verboser, beta_hat_unscaled, W_hat, lambda_W_hat, lambda_hat, denominator_hat) {
  verboser("", level = 3)
  verboser("-- DETAILED DIAGNOSTICS --", level = 3)
  if (length(beta_hat_unscaled)) {
    verboser(sprintf("  beta (response model, unscaled):"), level = 3)
    for (i in seq_along(beta_hat_unscaled)) {
      nm <- if (!is.null(names(beta_hat_unscaled))) names(beta_hat_unscaled)[i] else paste0("beta", i)
      verboser(sprintf("    %-25s %.6f", nm, beta_hat_unscaled[i]), level = 3)
    }
  }
  verboser(sprintf("  W (response rate):        %.6f", W_hat), level = 3)
  verboser(sprintf("  lambda_W (response multiplier): %.6f", lambda_W_hat), level = 3)
  if (length(lambda_hat)) {
    verboser(sprintf("  lambda_x (auxiliary multipliers):"), level = 3)
    for (i in seq_along(lambda_hat)) {
      nm <- if (!is.null(names(lambda_hat))) names(lambda_hat)[i] else paste0("lambda", i)
      verboser(sprintf("    %-25s %.6f", nm, lambda_hat[i]), level = 3)
    }
  }
  verboser(sprintf("  Denominator min:          %.6e", min(denominator_hat, na.rm = TRUE)), level = 3)
  verboser(sprintf("  Denominator median:       %.6f", stats::median(denominator_hat, na.rm = TRUE)), level = 3)
}

#' Log variance header and result
#' @keywords internal
el_log_variance_header <- function(verboser, variance_method, bootstrap_reps) {
  if (identical(variance_method, "none")) return(invisible(NULL))
  verboser("", level = 1)
  verboser("-- VARIANCE ESTIMATION --", level = 1)
  if (identical(variance_method, "bootstrap")) {
    verboser(sprintf("  Method:                   %s (reps = %d)", variance_method, bootstrap_reps), level = 2)
  } else {
    verboser(sprintf("  Method:                   %s", variance_method), level = 2)
  }
}

#' @keywords internal
el_log_variance_result <- function(verboser, se, elapsed) {
  verboser("", level = 2)
  if (is.finite(se)) {
    verboser(sprintf("  Standard error:           %.6f", se), level = 2, type = "result")
    verboser(sprintf("  Variance time:            %.3f seconds", elapsed), level = 2)
  } else {
    verboser("  Standard error:           NA (estimation failed)", level = 2, type = "result")
  }
}

#' Log final summary
#' @keywords internal
el_log_final <- function(verboser, y_hat, se) {
  verboser("", level = 1)
  el_log_banner(verboser, "EMPIRICAL LIKELIHOOD ESTIMATION COMPLETED")
  verboser("", level = 1)
  verboser(sprintf("  Estimate (y_hat):         %.6f", y_hat), level = 1, type = "result")
  if (is.finite(se)) {
    verboser(sprintf("  Standard error:           %.6f", se), level = 1, type = "result")
    verboser(sprintf("  95%% CI:                   [%.6f, %.6f]", y_hat - 1.96 * se, y_hat + 1.96 * se), level = 1, type = "result")
  } else {
    verboser("  Standard error:           NA", level = 1, type = "result")
  }
  verboser("", level = 1)
}
