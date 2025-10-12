#' Core Empirical Likelihood Estimator
#'
#' Implements the core computational engine for empirical likelihood estimation
#' under nonignorable nonresponse, including parameter solving, variance calculation,
#' and diagnostic computation.
#'
#' @param full_data Data frame or survey design object containing all units.
#' @param respondent_data Data frame containing only responding units.
#' @param respondent_weights Numeric vector of base sampling weights for respondents.
#' @param N_pop Numeric. Total population size (weighted if survey design).
#' @param internal_formula List of internal formulas for outcome, response, and auxiliary models.
#' @param auxiliary_means Named numeric vector of known population means.
#' @param standardize Logical. Whether to standardize predictors during estimation.
#' @param trim_cap Numeric. Upper bound for empirical likelihood weight trimming.
#' @param control List of control parameters for the nonlinear equation solver.
#' @param compute_score_variance_func Function to compute covariance of score totals.
#' @param on_failure Character. Action when solver fails: "return" or "error".
#' @param family List. Link function specification (typically logit).
#' @param variance_method Character. Variance estimation method.
#' @param bootstrap_reps Integer. Number of bootstrap replications.
#' @param variance_jacobian Character. Jacobian method for variance computation.
#' @param solver_jacobian Character. Jacobian method for equation solving.
#' @param variance_pseudoinverse Logical. Use pseudoinverse for singular matrices.
#' @param variance_ridge Logical. Use ridge regularization for variance.
#' @param user_args List. Original user arguments for bootstrap replication.
#' @param ... Additional arguments passed to the solver.
#'
#' @return List containing estimation results, diagnostics, and metadata.
#'
#' @details
#' This function implements the complete empirical likelihood estimation algorithm:
#' 1. Data preparation and scaling
#' 2. Equation system construction
#' 3. Multi-stage nonlinear equation solving with robust fallback strategies
#' 4. Weight computation and trimming
#' 5. Variance estimation (delta method or bootstrap)
#' 6. Comprehensive diagnostic computation
#'
#' @keywords internal
el_estimator_core <- function(full_data, respondent_data, respondent_weights, N_pop,
                              internal_formula, auxiliary_means, standardize,
                              trim_cap, control, compute_score_variance_func,
                              on_failure, family = logit_family(),
                              variance_method, bootstrap_reps,
                              variance_jacobian = c("auto", "analytic", "numeric"),
                              solver_jacobian = c("auto", "analytic", "none"),
                              solver_method = c("auto", "newton", "broyden"),
                              solver_args = list(),
                              variance_pseudoinverse = FALSE,
                              variance_ridge = FALSE,
                              user_args, ...) {
  variance_jacobian <- match.arg(variance_jacobian)
  solver_jacobian <- match.arg(solver_jacobian)
  solver_method <- match.arg(solver_method)

# 0. Setup
  force(family)
  outcome_var <- all.vars(internal_formula$outcome)[1]

# 1. Data Preparation
  has_aux <- !is.null(internal_formula$auxiliary)

  response_model_formula <- update(internal_formula$response, NULL ~ .)
  response_model_matrix_unscaled <- model.matrix(response_model_formula, data = respondent_data)

  if (has_aux) {
# Build auxiliary design on respondents
    auxiliary_matrix_unscaled <- model.matrix(internal_formula$auxiliary, data = respondent_data)

# If user did not supply population moments, infer from the full sample
    if (is.null(auxiliary_means)) {
      if (inherits(full_data, "survey.design")) {
        mm_full <- model.matrix(internal_formula$auxiliary, data = full_data$variables)
        w_full <- as.numeric(weights(full_data))
# Align columns with respondent design (drop absent levels in respondents)
        common_cols <- intersect(colnames(auxiliary_matrix_unscaled), colnames(mm_full))
        if (length(common_cols) == 0L) {
# No overlap; disable auxiliaries gracefully
          auxiliary_matrix_unscaled <- matrix(nrow = nrow(respondent_data), ncol = 0)
          mu_x_unscaled <- NULL
          has_aux <- FALSE
        } else {
          mm_full <- mm_full[, common_cols, drop = FALSE]
# Weighted means of model-matrix columns
          mu <- as.numeric(colSums(mm_full * w_full) / sum(w_full))
          names(mu) <- colnames(mm_full)
# Reorder to respondent design
          auxiliary_matrix_unscaled <- auxiliary_matrix_unscaled[, common_cols, drop = FALSE]
          mu_x_unscaled <- mu[colnames(auxiliary_matrix_unscaled)]
          message("No auxiliary_means supplied; inferred from design-weighted sample and treated as fixed.")
        }
      } else {
        mm_full <- model.matrix(internal_formula$auxiliary, data = full_data)
        common_cols <- intersect(colnames(auxiliary_matrix_unscaled), colnames(mm_full))
        if (length(common_cols) == 0L) {
          auxiliary_matrix_unscaled <- matrix(nrow = nrow(respondent_data), ncol = 0)
          mu_x_unscaled <- NULL
          has_aux <- FALSE
        } else {
          mm_full <- mm_full[, common_cols, drop = FALSE]
          mu <- colMeans(mm_full)
# Ensure naming and order
          auxiliary_matrix_unscaled <- auxiliary_matrix_unscaled[, common_cols, drop = FALSE]
          mu_x_unscaled <- as.numeric(mu[colnames(auxiliary_matrix_unscaled)])
          names(mu_x_unscaled) <- colnames(auxiliary_matrix_unscaled)
          message("No auxiliary_means supplied; inferred from sample and treated as fixed.")
        }
      }
    } else {
# Use provided means; restrict/reorder to match respondent design
      keep <- intersect(names(auxiliary_means), colnames(auxiliary_matrix_unscaled))
      auxiliary_matrix_unscaled <- auxiliary_matrix_unscaled[, keep, drop = FALSE]
      mu_x_unscaled <- as.numeric(auxiliary_means[keep])
      names(mu_x_unscaled) <- keep
      if (ncol(auxiliary_matrix_unscaled) == 0L) {
        has_aux <- FALSE
        mu_x_unscaled <- NULL
      }
    }
  } else {
    auxiliary_matrix_unscaled <- matrix(nrow = nrow(respondent_data), ncol = 0)
    mu_x_unscaled <- NULL
  }

# 2. Scaling
  scaling_result <- validate_and_apply_nmar_scaling(
    standardize = standardize,
    has_aux = has_aux,
    response_model_matrix_unscaled = response_model_matrix_unscaled,
    auxiliary_matrix_unscaled = auxiliary_matrix_unscaled,
    mu_x_unscaled = mu_x_unscaled,
    weights = respondent_weights
  )
  nmar_scaling_recipe <- scaling_result$nmar_scaling_recipe
  response_model_matrix_scaled <- scaling_result$response_model_matrix_scaled
  auxiliary_matrix_scaled <- scaling_result$auxiliary_matrix_scaled
  mu_x_scaled <- scaling_result$mu_x_scaled

# 3. Build Solver Components
  n_resp_weighted <- sum(respondent_weights)
  equation_system_func <- el_build_equation_system(
    family = family, response_model_matrix = response_model_matrix_scaled, auxiliary_matrix = auxiliary_matrix_scaled,
    respondent_weights = respondent_weights, N_pop = N_pop, n_resp_weighted = n_resp_weighted, mu_x_scaled = mu_x_scaled
  )
  analytical_jac_func <- build_el_jacobian(
    family = family, response_model_matrix = response_model_matrix_scaled, auxiliary_matrix = auxiliary_matrix_scaled,
    respondent_weights = respondent_weights, N_pop = N_pop, n_resp_weighted = n_resp_weighted, mu_x_scaled = mu_x_scaled
  )

# 4. Solve for Parameters with Automated Multi-Start Strategy
  K_beta <- ncol(response_model_matrix_scaled)
  K_aux <- ncol(auxiliary_matrix_scaled)
  init_beta <- rep(0, K_beta)
  init <- c(init_beta, 0, rep(0, K_aux))
  final_control <- modifyList(list(ftol = 1e-8, xtol = 1e-8, maxit = 100), control)
  top_args_ctrl <- extract_nleqslv_top(control)
  top_args <- merge_nleqslv_top(solver_args, top_args_ctrl)
  final_control <- sanitize_nleqslv_control(final_control)
  use_solver_jac <- solver_jacobian %in% c("auto", "analytic") && !is.null(analytical_jac_func)

# Instrumentation for timing and solver used
  t_solve_start <- proc.time()[[3]]
  solver_out <- el_run_solver(
    equation_system_func = equation_system_func,
    analytical_jac_func = analytical_jac_func,
    init = init,
    final_control = final_control,
    top_args = top_args,
    solver_method = solver_method,
    use_solver_jac = use_solver_jac,
    K_beta = K_beta,
    K_aux = K_aux,
    respondent_weights = respondent_weights,
    N_pop = N_pop
  )
  solution <- solver_out$solution
  solver_method_used <- solver_out$method
  nleqslv_global_used <- if (!is.null(solver_out$used_top$global)) solver_out$used_top$global else NA_character_
  nleqslv_xscalm_used <- if (!is.null(solver_out$used_top$xscalm)) solver_out$used_top$xscalm else NA_character_

  if (any(is.na(solution$x)) || solution$termcd > 2) {
    if (on_failure == "error") {
      stop(convergenceError(paste("Solver failed to converge after multiple attempts:", solution$message)))
    } else {
      return(list(
        converged = FALSE,
        diagnostics = list(convergence_code = solution$termcd, message = solution$message),
        nmar_scaling_recipe = nmar_scaling_recipe
      ))
    }
  }

# 5. Post-processing and Point Estimate Calculation
  estimates <- solution$x
  solver_time <- proc.time()[[3]] - t_solve_start
  post <- el_post_solution(
    estimates = estimates,
    response_model_matrix_scaled = response_model_matrix_scaled,
    response_model_matrix_unscaled = response_model_matrix_unscaled,
    auxiliary_matrix_scaled = auxiliary_matrix_scaled,
    mu_x_scaled = mu_x_scaled,
    respondent_data = respondent_data,
    outcome_var = outcome_var,
    family = family,
    N_pop = N_pop,
    respondent_weights = respondent_weights,
    K_beta = K_beta,
    K_aux = K_aux,
    nmar_scaling_recipe = nmar_scaling_recipe,
    standardize = standardize,
    trim_cap = trim_cap
  )
  if (post$error) {
    if (on_failure == "error") {
      stop(convergenceError(post$message))
    } else {
      return(list(
        converged = FALSE,
        diagnostics = list(convergence_code = -1, message = post$message),
        nmar_scaling_recipe = nmar_scaling_recipe
      ))
    }
  }
  y_hat <- post$y_hat
  p_i <- post$weights
  beta_hat_unscaled <- post$beta_hat_unscaled
  W_hat <- post$W_hat
  lambda_hat <- post$lambda_hat
  eta_i_hat <- post$eta_i_hat
  w_i_hat <- post$w_i_hat
  denominator_hat <- post$denominator_hat
  lambda_W_hat <- post$lambda_W_hat

# 6. Diagnostics at Solution
  eq_residuals <- tryCatch(equation_system_func(estimates), error = function(e) rep(NA_real_, length(estimates)))
  max_eq_resid <- suppressWarnings(max(abs(eq_residuals), na.rm = TRUE))
  jac_sel <- el_select_variance_jacobian(
    equation_system_func = equation_system_func,
    analytical_jac_func = analytical_jac_func,
    estimates = estimates,
    variance_jacobian = variance_jacobian
  )
  A_source <- jac_sel$A_source
  A_condition <- jac_sel$A_condition
  A_diff_norm <- jac_sel$A_diff_norm
  jacobian_auto_rule <- jac_sel$jacobian_auto_rule
  denom_stats <- list(min = suppressWarnings(min(denominator_hat, na.rm = TRUE)), p_small = mean(denominator_hat < 1e-6))
  p_untrim <- respondent_weights / denominator_hat
  cons <- constraint_summaries(w_i_hat, W_hat, p_untrim, if (K_aux > 0) sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-") else NULL)
  constraint_eqW_sum <- cons$constraint_sum_W
  constraint_aux_sum <- cons$constraint_sum_aux

# 7. Conditional Variance Calculation
  se_y_hat <- NA
  vcov_unscaled <- NA
  vcov_message <- "Calculation successful"
  suppress_var_warn <- isTRUE(user_args$suppress_warnings)
  if (variance_method == "delta" && is.finite(trim_cap) && !suppress_var_warn) {
    warning("Delta method variance is not recommended with weight trimming. Consider variance_method = 'bootstrap'.", call. = FALSE)
  }

  t_var_start <- proc.time()[[3]]
  if (variance_method == "bootstrap") {
    user_args_internal <- user_args
    user_args_internal$variance_method <- "none"
    user_args_internal$suppress_warnings <- TRUE
    boot_args <- c(list(data = full_data, estimator_func = el, point_estimate = y_hat, bootstrap_reps = bootstrap_reps), user_args_internal)
    boot_try <- tryCatch(list(result = do.call(bootstrap_variance, boot_args), message = "Calculation successful"), error = function(e) list(result = NULL, message = paste("Bootstrap failed:", e$message)))
    vcov_message <- boot_try$message
    if (!is.null(boot_try$result)) {
      se_y_hat <- boot_try$result$se
    } else {
      se_y_hat <- NA_real_
    }
  } else if (variance_method == "delta") {
    var_out <- el_variance_delta(
      equation_system_func = equation_system_func,
      analytical_jac_func = analytical_jac_func,
      estimates = estimates,
      variance_jacobian = variance_jacobian,
      family = family,
      response_model_matrix_scaled = response_model_matrix_scaled,
      response_model_matrix_unscaled = response_model_matrix_unscaled,
      auxiliary_matrix_scaled = auxiliary_matrix_scaled,
      mu_x_scaled = mu_x_scaled,
      eta_i_hat = eta_i_hat,
      w_i_hat = w_i_hat,
      W_hat = W_hat,
      denominator_hat = denominator_hat,
      lambda_W_hat = lambda_W_hat,
      full_data = full_data,
      compute_score_variance_func = compute_score_variance_func,
      respondent_weights = respondent_weights,
      N_pop = N_pop,
      n_resp_weighted = sum(respondent_weights),
      trim_cap = trim_cap,
      outcome_vec = respondent_data[[outcome_var]],
      K_beta = K_beta,
      standardize = standardize,
      nmar_scaling_recipe = nmar_scaling_recipe,
      variance_ridge = variance_ridge,
      variance_pseudoinverse = variance_pseudoinverse
    )
    se_y_hat <- var_out$se_y_hat
    vcov_unscaled <- var_out$vcov_unscaled
    vcov_message <- var_out$vcov_message
    used_pseudoinverse <- var_out$diag$used_pseudoinverse
    used_ridge <- var_out$diag$used_ridge
    invert_rule <- var_out$diag$invert_rule
# overwrite diagnostics from selection to ensure consistency (no change if not recalculated)
    A_condition <- var_out$diag$A_condition
    A_source <- var_out$diag$A_source
    A_diff_norm <- var_out$diag$A_diff_norm
    jacobian_auto_rule <- var_out$diag$jacobian_auto_rule
    auto_var_rule <- var_out$diag$variance_auto_rule
  } else if (variance_method == "none") {
# Skip variance calculation entirely
    se_y_hat <- NA_real_
    vcov_unscaled <- matrix(NA_real_, K_beta, K_beta, dimnames = list(colnames(response_model_matrix_unscaled), colnames(response_model_matrix_unscaled)))
    vcov_message <- "Variance skipped (variance_method='none')"
  }
  variance_time <- proc.time()[[3]] - t_var_start

# 8. Return Final Results List
  return(list(
    y_hat = y_hat, se = se_y_hat, weights = p_i,
    coefficients = list(response_model = beta_hat_unscaled, nuisance = list(W_hat = W_hat, lambda_x = lambda_hat)),
    vcov = vcov_unscaled, converged = TRUE,
    diagnostics = list(
      convergence_code = solution$termcd,
      message = solution$message,
      vcov_message = vcov_message,
      trimmed_fraction = post$trimmed_fraction,
      solver_jacobian = if (use_solver_jac) "analytic" else "none",
      solver_method = solver_method_used,
      nleqslv_global = nleqslv_global_used,
      nleqslv_xscalm = nleqslv_xscalm_used,
      solver_iterations = if (!is.null(solution$iter)) solution$iter else NA,
      solver_time = solver_time,
      variance_time = variance_time,
      reparam_W = "logit",
      max_equation_residual = max_eq_resid,
      jacobian_condition_number = A_condition,
      jacobian_source = A_source,
      jacobian_rel_diff = A_diff_norm,
      jacobian_auto_rule = jacobian_auto_rule,
      min_denominator = denom_stats$min,
      fraction_small_denominators = denom_stats$p_small,
      constraint_sum_W = constraint_eqW_sum,
      constraint_sum_aux = constraint_aux_sum,
      used_pseudoinverse = if (exists("used_pseudoinverse")) used_pseudoinverse else FALSE,
      used_ridge = if (exists("used_ridge")) used_ridge else FALSE,
      invert_rule = if (exists("invert_rule")) invert_rule else NA_character_,
      variance_auto_rule = if (exists("auto_var_rule")) auto_var_rule else NA_character_
    ),
    nmar_scaling_recipe = nmar_scaling_recipe, fitted_values = drop(w_i_hat)
  ))
}
