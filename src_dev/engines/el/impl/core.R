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
                              variance_pseudoinverse = FALSE,
                              variance_ridge = FALSE,
                              user_args, ...) {
  variance_jacobian <- match.arg(variance_jacobian)
  solver_jacobian <- match.arg(solver_jacobian)

  # 0. Setup
  force(family)
  outcome_var <- all.vars(internal_formula$outcome)[1]

  # 1. Data Preparation
  has_aux <- !is.null(internal_formula$auxiliary)
  if (has_aux && is.null(auxiliary_means)) {
    message("An `auxiliary` formula was created but no `auxiliary_means` were given. Ignoring auxiliary information.")
    has_aux <- FALSE
  }

  response_model_formula <- update(internal_formula$response, NULL ~ .)
  response_model_matrix_unscaled <- model.matrix(response_model_formula, data = respondent_data)
  auxiliary_matrix_unscaled <- if (has_aux) model.matrix(internal_formula$auxiliary, data = respondent_data) else matrix(nrow = nrow(respondent_data), ncol = 0)
  mu_x_unscaled <- if (has_aux) auxiliary_means else NULL

  # 2. Scaling
  scaling_result <- validate_and_apply_nmar_scaling(
    standardize = standardize,
    has_aux = has_aux,
    response_model_matrix_unscaled = response_model_matrix_unscaled,
    auxiliary_matrix_unscaled = auxiliary_matrix_unscaled,
    mu_x_unscaled = mu_x_unscaled
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
  final_control <- modifyList(list(ftol = 1e-10, xtol = 1e-10, maxit = 100), control)
  use_solver_jac <- solver_jacobian %in% c("auto", "analytic") && !is.null(analytical_jac_func)

  # STAGE 1: Newton method
  solution <- nleqslv::nleqslv(
    x = init,
    fn = equation_system_func,
    jac = if (use_solver_jac) analytical_jac_func else NULL,
    method = "Newton",
    control = final_control, ...
  )

  # STAGE 2: Perturbation-based recovery
  if (any(is.na(solution$x)) || solution$termcd > 2) {
    for (i in seq_len(3)) {
      init_beta_perturbed <- init_beta + rnorm(K_beta, mean = 0, sd = 0.5)
      # Parameterization uses z = logit(W). Seed z near observed response rate.
      W_seed <- sum(respondent_weights) / N_pop
      W_seed <- min(max(W_seed, 1e-12), 1 - 1e-12)
      z_seed <- stats::qlogis(W_seed)
      init_perturbed <- c(init_beta_perturbed, z_seed, rep(0, K_aux))
      solution <- nleqslv::nleqslv(
        x = init_perturbed,
        fn = equation_system_func,
        jac = if (use_solver_jac) analytical_jac_func else NULL,
        method = "Newton",
        control = final_control, ...
      )
      if (!any(is.na(solution$x)) && solution$termcd <= 2) break
    }
  }

  # STAGE 3: Broyden fallback
  if (any(is.na(solution$x)) || solution$termcd > 2) {
    broyden_control <- final_control
    if (!is.null(broyden_control$maxit) && is.finite(broyden_control$maxit) && broyden_control$maxit < 5) {
      broyden_control$maxit <- 50
    }
    solution <- nleqslv::nleqslv(
      x = init,
      fn = equation_system_func,
      jac = NULL,
      method = "Broyden",
      control = broyden_control, ...
    )
  }

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
  beta_hat_scaled <- estimates[1:K_beta]
  names(beta_hat_scaled) <- colnames(response_model_matrix_scaled)
  W_hat <- plogis(estimates[K_beta + 1])
  lambda_hat <- if (K_aux > 0) estimates[(K_beta + 2):length(estimates)] else numeric(0)
  eta_i_hat <- as.vector(response_model_matrix_scaled %*% beta_hat_scaled)
  w_i_hat <- family$linkinv(eta_i_hat)

  lambda_W_hat <- ((N_pop / sum(respondent_weights)) - 1) / (1 - W_hat)
  denominator_hat <- 1 + lambda_W_hat * (w_i_hat - W_hat)
  if (K_aux > 0) denominator_hat <- denominator_hat + as.vector(sweep(auxiliary_matrix_scaled, 2, mu_x_scaled, "-") %*% lambda_hat)
  p_i_untrimmed <- respondent_weights / denominator_hat

  TOL <- 1e-8
  min_w <- suppressWarnings(min(p_i_untrimmed, na.rm = TRUE))
  if (is.finite(min_w) && min_w < -TOL) {
    msg <- paste0(
      "Negative EL weights produced (min = ", round(min_w, 6), "). This often indicates that the auxiliary means are ",
      "inconsistent with the sample data, or the model has failed to converge to a valid solution."
    )
    if (on_failure == "error") {
      stop(convergenceError(msg))
    } else {
      return(list(
        converged = FALSE,
        diagnostics = list(convergence_code = -1, message = msg),
        nmar_scaling_recipe = nmar_scaling_recipe
      ))
    }
  }
  if (any(p_i_untrimmed < 0)) p_i_untrimmed[p_i_untrimmed < 0] <- 0

  trim_results <- trim_weights(p_i_untrimmed, cap = trim_cap)
  p_i <- trim_results$weights
  y_hat <- sum(p_i * respondent_data[[outcome_var]]) / sum(p_i)

  beta_hat_unscaled <- if (standardize) unscale_coefficients(beta_hat_scaled, matrix(0, K_beta, K_beta), nmar_scaling_recipe)$coefficients else beta_hat_scaled
  names(beta_hat_unscaled) <- colnames(response_model_matrix_unscaled)

  # 6. Diagnostics at Solution
  eq_residuals <- tryCatch(equation_system_func(estimates), error = function(e) rep(NA_real_, length(estimates)))
  max_eq_resid <- suppressWarnings(max(abs(eq_residuals), na.rm = TRUE))
  jac_choice <- choose_jacobian(
    analytic_fun = analytical_jac_func,
    numeric_fun = function(x) numDeriv::jacobian(func = equation_system_func, x = x),
    at = estimates
  )
  A_matrix_analytic <- jac_choice$A_analytic
  A_matrix_numeric <- jac_choice$A_numeric
  A_diff_norm <- jac_choice$rel_diff
  A_condition_numeric <- tryCatch(if (!is.null(A_matrix_numeric)) kappa(A_matrix_numeric) else NA_real_, error = function(e) NA_real_)
  A_condition_analytic <- tryCatch(if (!is.null(A_matrix_analytic)) kappa(A_matrix_analytic) else NA_real_, error = function(e) NA_real_)
  jacobian_auto_rule <- "default"
  if (variance_jacobian == "numeric") {
    A_matrix_var <- A_matrix_numeric
    A_condition <- A_condition_numeric
    A_source <- "numeric"
  } else if (variance_jacobian == "analytic") {
    A_matrix_var <- A_matrix_analytic
    A_condition <- A_condition_analytic
    A_source <- "analytic"
  } else {
    # auto: prefer analytic when available, unless quality gates suggest numeric
    REL_THR <- 1e-3
    KAPPA_RATIO_THR <- 10
    if (!is.null(A_matrix_analytic)) {
      # Start with analytic as default
      A_matrix_var <- A_matrix_analytic
      A_condition <- A_condition_analytic
      A_source <- "analytic"
      # If numeric also available, evaluate quality gates
      if (!is.null(A_matrix_numeric)) {
        # Gate 1: large relative difference between analytic and numeric
        if (is.finite(A_diff_norm) && A_diff_norm > REL_THR) {
          A_matrix_var <- A_matrix_numeric
          A_condition <- A_condition_numeric
          A_source <- "numeric"
          jacobian_auto_rule <- "rel_diff_high"
        } else {
          # Gate 2: analytic much worse conditioned than numeric
          if (is.finite(A_condition_analytic) && is.finite(A_condition_numeric) &&
            A_condition_analytic > KAPPA_RATIO_THR * A_condition_numeric) {
            A_matrix_var <- A_matrix_numeric
            A_condition <- A_condition_numeric
            A_source <- "numeric"
            jacobian_auto_rule <- "kappa_ratio_high"
          }
        }
      }
    } else if (!is.null(A_matrix_numeric)) {
      A_matrix_var <- A_matrix_numeric
      A_condition <- A_condition_numeric
      A_source <- "numeric"
    } else {
      # Fallback to whatever choose_jacobian provided
      A_matrix_var <- jac_choice$A
      A_source <- jac_choice$source
      A_condition <- jac_choice$kappa
    }
  }
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

 if (variance_method == "bootstrap") {
    user_args_internal <- user_args
    user_args_internal$variance_method <- "delta"
    user_args_internal$suppress_warnings <- TRUE
    boot_args <- c(list(data = full_data, estimator_func = el, point_estimate = y_hat, bootstrap_reps = bootstrap_reps), user_args_internal)
    boot_try <- tryCatch(list(result = do.call(bootstrap_variance, boot_args), message = "Calculation successful"), error = function(e) list(result = NULL, message = paste("Bootstrap failed:", e$message)))
    vcov_message <- boot_try$message
    if (!is.null(boot_try$result)) se_y_hat <- boot_try$result$se
    else se_y_hat <- NA_real_
  } else {
    vcov_unscaled <- matrix(NA, K_beta, K_beta, dimnames = list(colnames(response_model_matrix_unscaled), colnames(response_model_matrix_unscaled)))
    vcov_try <- tryCatch(
      {
        A_matrix <- if (!is.null(A_matrix_var)) A_matrix_var else if (!is.null(analytical_jac_func)) analytical_jac_func(estimates) else numDeriv::jacobian(func = equation_system_func, x = estimates)
        res <- el_compute_delta_variance(
          A_matrix = A_matrix,
          family = family,
          response_model_matrix_scaled = response_model_matrix_scaled,
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
          estimates = estimates,
          variance_ridge = variance_ridge,
          variance_pseudoinverse = variance_pseudoinverse
        )
        list(result = res, message = "Calculation successful")
      },
      error = function(e) list(result = NULL, message = e$message)
    )
    vcov_message <- vcov_try$message
    vcov_result <- vcov_try$result
    need_var_fb <- (is.null(vcov_result) || !is.finite(as.numeric(vcov_result$var_y_hat))) && (!is.null(analytical_jac_func)) && (variance_jacobian == "numeric" || (variance_jacobian == "auto" && A_source == "numeric"))
    if (need_var_fb) {
      vcov_try_fb <- tryCatch(
        {
          A_matrix_alt <- analytical_jac_func(estimates)
          res <- el_compute_delta_variance(
            A_matrix = A_matrix_alt,
            family = family,
            response_model_matrix_scaled = response_model_matrix_scaled,
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
            estimates = estimates,
            variance_ridge = variance_ridge,
            variance_pseudoinverse = variance_pseudoinverse
          )
          list(result = res, message = "Calculation successful")
        },
        error = function(e) list(result = vcov_result, message = vcov_message)
      )
      vcov_result <- vcov_try_fb$result
      vcov_message <- vcov_try_fb$message
    }
    if (!is.null(vcov_result)) {
      var_y_hat_val <- as.numeric(vcov_result$var_y_hat)
      if (is.finite(var_y_hat_val)) se_y_hat <- sqrt(pmax(var_y_hat_val, 0))
      vcov_beta_scaled <- vcov_result$vcov_matrix_sandwich_scaled[1:K_beta, 1:K_beta, drop = FALSE]
      vcov_unscaled <- if (standardize) unscale_coefficients(beta_hat_scaled, vcov_beta_scaled, nmar_scaling_recipe)$vcov else vcov_beta_scaled
      used_pseudoinverse <- isTRUE(vcov_result$used_pseudoinverse)
      used_ridge <- isTRUE(vcov_result$used_ridge)
      invert_rule <- if (!is.null(vcov_result$invert_rule)) vcov_result$invert_rule else NA_character_
   }
    if (!is.finite(se_y_hat)) se_y_hat <- NA_real_
  }

  # 8. Return Final Results List
  return(list(
    y_hat = y_hat, se = se_y_hat, weights = p_i,
    coefficients = list(response_model = beta_hat_unscaled, nuisance = list(W_hat = W_hat, lambda_x = lambda_hat)),
    vcov = vcov_unscaled, converged = TRUE,
    diagnostics = list(
      convergence_code = solution$termcd,
      message = solution$message,
      vcov_message = vcov_message,
      trimmed_fraction = trim_results$trimmed_fraction,
      solver_jacobian = if (use_solver_jac) "analytic" else "none",
      solver_method = if (!is.null(solution$method)) solution$method else NA,
      solver_iterations = if (!is.null(solution$iter)) solution$iter else NA,
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
      invert_rule = if (exists("invert_rule")) invert_rule else NA_character_
    ),
    nmar_scaling_recipe = nmar_scaling_recipe, fitted_values = drop(w_i_hat)
  ))
}
