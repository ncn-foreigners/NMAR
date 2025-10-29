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
#' @param user_args List. Original user arguments for bootstrap replication.
#' @param ... Additional arguments passed to the solver.
#'
#' @return List containing estimation results, diagnostics, and metadata.
#'
#' @details
#' Implements the complete empirical likelihood pipeline under nonignorable
#' nonresponse with optional auxiliary moment constraints. The solver uses the
#' full stacked system in \eqn{(\beta, z, \lambda_x)} with an analytic Jacobian
#' and globalization via `nleqslv`. Numerical safeguards (denominator
#' positivity guards, predictor standardization, and stable linear algebra)
#' are applied. After solving, weights are constructed with denominator guards
#' and optional trimming; delta variance is computed by a sandwich formula
#' evaluated at the solution or via bootstrap when requested.
#'
#' Steps
#' 1. Data preparation and scaling
#' 2. Equation system construction
#' 3. Nonlinear equation solving with safeguards
#' 4. Weight construction and optional trimming
#' 5. Variance estimation
#' 6. Diagnostic computation
#'
#' @keywords internal
el_estimator_core <- function(full_data, respondent_data, respondent_weights, N_pop,
                              internal_formula, auxiliary_means, standardize,
                              trim_cap, control, compute_score_variance_func,
                              on_failure, family = logit_family(),
                              variance_method, bootstrap_reps,
                              user_args, start = NULL, ...) {

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
# Auxiliary means inferred from design-weighted sample and treated as fixed
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
# Auxiliary means inferred from sample and treated as fixed
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

# Heuristic inconsistency check for user-supplied auxiliary means:
# Compare provided means to respondents' sample means; if far (|z|>threshold),
# warn and record diagnostics.
  aux_inconsistency_max_z <- NA_real_
  aux_inconsistency_cols <- character(0)
  if (has_aux && !is.null(auxiliary_means)) {
    aux_df <- model.matrix(internal_formula$auxiliary, data = respondent_data)
# Drop columns with near-zero variance to avoid division by zero
    bad_sd <- apply(aux_df, 2, function(col) stats::sd(col) < 1e-8)
    if (any(bad_sd)) aux_df <- aux_df[, !bad_sd, drop = FALSE]
    if (ncol(aux_df) > 0) {
      sample_means <- colMeans(aux_df)
      sample_sds <- apply(aux_df, 2, stats::sd)
      means_provided <- auxiliary_means[colnames(aux_df)]
      z_scores <- abs((means_provided - sample_means) / pmax(sample_sds, 1e-8))
      thr <- getOption("nmar.el_aux_z_threshold", 8)
      if (!is.numeric(thr) || length(thr) != 1L || !is.finite(thr) || thr <= 0) thr <- 8
      aux_inconsistency_max_z <- suppressWarnings(max(z_scores, na.rm = TRUE))
      flag <- is.finite(aux_inconsistency_max_z) && aux_inconsistency_max_z > thr
      if (flag) {
        aux_inconsistency_cols <- names(z_scores)[which(is.finite(z_scores) & z_scores > thr)]
        warning(sprintf(
          "Auxiliary means appear far from respondents' support (max |z| = %.2f, threshold = %.2f). Proceeding; see diagnostics.",
          aux_inconsistency_max_z, thr
        ), call. = FALSE)
      }
    }
  }
# Build full stacked builders
  equation_system_func <- el_build_equation_system(
    family = family, response_model_matrix = response_model_matrix_scaled, auxiliary_matrix = auxiliary_matrix_scaled,
    respondent_weights = respondent_weights, N_pop = N_pop, n_resp_weighted = n_resp_weighted, mu_x_scaled = mu_x_scaled
  )
  analytical_jac_func <- build_el_jacobian(
    family = family, response_model_matrix = response_model_matrix_scaled, auxiliary_matrix = auxiliary_matrix_scaled,
    respondent_weights = respondent_weights, N_pop = N_pop, n_resp_weighted = n_resp_weighted, mu_x_scaled = mu_x_scaled
  )
# 4. Solve for Parameters with a safeguarded Newton method
  K_beta <- ncol(response_model_matrix_scaled)
  K_aux <- ncol(auxiliary_matrix_scaled)
# Optional user-provided starting values on the original scale
  init_beta <- rep(0, K_beta)
  names(init_beta) <- colnames(response_model_matrix_scaled)
  init_lambda <- rep(0, K_aux)
  names(init_lambda) <- colnames(auxiliary_matrix_scaled)
# Initialize z at the empirical response rate
  init_z <- {
    W0 <- sum(respondent_weights) / N_pop
    W0 <- min(max(W0, 1e-12), 1 - 1e-12)
    stats::qlogis(W0)
  }
  if (!is.null(start) && is.list(start)) {
    if (!is.null(start$beta)) {
# Warn on unmatched coefficient names
      beta_names <- names(start$beta)
      if (!is.null(beta_names)) {
        unknown_beta <- setdiff(beta_names, colnames(response_model_matrix_scaled))
        if (length(unknown_beta) > 0) {
          warning(sprintf("Start 'beta' contains names not in the response model matrix and will be ignored: %s",
                          paste(unknown_beta, collapse = ", ")),
                  call. = FALSE)
        }
      }
      init_beta <- scale_coefficients(start$beta, nmar_scaling_recipe, colnames(response_model_matrix_scaled))
    }
    if (!is.null(start$z)) {
      init_z <- as.numeric(start$z)[1]
      if (!is.finite(init_z)) init_z <- 0
    } else if (!is.null(start$W)) {
      W0 <- min(max(as.numeric(start$W)[1], 1e-12), 1 - 1e-12)
      init_z <- stats::qlogis(W0)
    }
    if (K_aux > 0 && !is.null(start$lambda)) {
# Warn on unmatched lambda names
      lambda_names <- names(start$lambda)
      if (!is.null(lambda_names)) {
        unknown_lambda <- setdiff(lambda_names, colnames(auxiliary_matrix_scaled))
        if (length(unknown_lambda) > 0) {
          warning(sprintf("Start 'lambda' contains names not in the auxiliary design and will be ignored: %s",
                          paste(unknown_lambda, collapse = ", ")),
                  call. = FALSE)
        }
      }
      init_lambda <- scale_aux_multipliers(start$lambda, nmar_scaling_recipe, colnames(auxiliary_matrix_scaled))
    } else if (K_aux == 0 && !is.null(start$lambda)) {
      warning("Start 'lambda' provided, but the current model has no auxiliary constraints; ignoring.", call. = FALSE)
    }
  }
  init <- c(unname(init_beta), init_z, unname(init_lambda))
# Split user control into top-level nleqslv args (global/xscalm) and control list
  control_top <- validate_nleqslv_top(extract_nleqslv_top(control))
  final_control <- modifyList(list(ftol = 1e-8, xtol = 1e-8, maxit = 100, trace = FALSE), control)
  final_control <- sanitize_nleqslv_control(final_control)

# Instrumentation for timing and solver used
  t_solve_start <- proc.time()[[3]]
  {
# Solve full stacked equation system with analytic Jacobian
  solver_out <- el_run_solver(
    equation_system_func = equation_system_func,
    analytical_jac_func = analytical_jac_func,
    init = init,
    final_control = final_control,
    top_args = control_top,
    solver_method = "auto",
    use_solver_jac = TRUE,
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
      stop(convergenceError(paste("Solver failed to converge:", solution$message)))
    } else {
      return(list(
        converged = FALSE,
        diagnostics = list(
          convergence_code = solution$termcd,
          message = solution$message,
          aux_inconsistency_max_z = aux_inconsistency_max_z,
          aux_inconsistency_cols = aux_inconsistency_cols
        ),
        nmar_scaling_recipe = nmar_scaling_recipe
      ))
    }
  }
  }

# 5. Post-processing and Point Estimate Calculation
  estimates <- solution$x
  beta_hat_scaled <- estimates[1:K_beta]
  lambda_hat <- if (K_aux > 0) estimates[(K_beta + 2):(K_beta + 1 + K_aux)] else numeric(0)
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
  A_condition <- tryCatch({ kappa(analytical_jac_func(estimates)) }, error = function(e) NA_real_)
  denom_stats <- list(min = suppressWarnings(min(denominator_hat, na.rm = TRUE)), p_small = mean(denominator_hat < 1e-6))
  p_untrim <- respondent_weights / denominator_hat
  Xc_centered_diag <- NULL
  if (K_aux > 0) {
    mu_match <- as.numeric(mu_x_scaled[colnames(auxiliary_matrix_scaled)])
    Xc_centered_diag <- sweep(auxiliary_matrix_scaled, 2, mu_match, "-")
  }
  cons <- tryCatch(
    constraint_summaries(w_i_hat, W_hat, p_untrim, Xc_centered_diag),
    error = function(e) {
# Fallback: derive constraint sums directly from the stacked equations at the estimate
      eq_vals <- tryCatch(equation_system_func(estimates), error = function(e2) NULL)
      out <- list(constraint_sum_W = NA_real_, constraint_sum_aux = numeric(0))
      if (!is.null(eq_vals)) {
        idx_W <- K_beta + 1L
        out$constraint_sum_W <- as.numeric(eq_vals[idx_W])
        if (K_aux > 0 && length(eq_vals) >= (idx_W + K_aux)) {
          aux_vals <- as.numeric(eq_vals[(idx_W + 1):(idx_W + K_aux)])
          names(aux_vals) <- colnames(auxiliary_matrix_scaled)
          out$constraint_sum_aux <- aux_vals
        }
      }
      out
    }
  )
  constraint_eqW_sum <- cons$constraint_sum_W
  constraint_aux_sum <- cons$constraint_sum_aux

# Additional diagnostics: denominators and weight concentration
  denom_q <- tryCatch(stats::quantile(denominator_hat, probs = c(0.01, 0.05, 0.5), na.rm = TRUE), error = function(e) c(`1%` = NA_real_, `5%` = NA_real_, `50%` = NA_real_))
  denom_cnt_1e4 <- sum(denominator_hat < 1e-4)
  weight_sum <- sum(p_i)
  weight_share <- if (weight_sum > 0) sort(p_i / weight_sum, decreasing = TRUE) else rep(NA_real_, length(p_i))
  weight_max_share <- if (length(weight_share)) weight_share[1] else NA_real_
  weight_top5_share <- if (length(weight_share) >= 5) sum(weight_share[1:5]) else sum(weight_share, na.rm = TRUE)
  weight_ess <- if (weight_sum > 0) (weight_sum^2) / sum(p_i^2) else NA_real_

# 7. Conditional Variance Calculation
  se_y_hat <- NA
  vcov_unscaled <- NA
  vcov_message <- "Calculation successful"
  if (variance_method == "delta" && is.finite(trim_cap)) {
    warning("Delta method variance is not recommended with weight trimming. Consider variance_method = 'bootstrap'.", call. = FALSE)
  }

  t_var_start <- proc.time()[[3]]
# Defaults for variance diagnostics
  diag_grad_source <- NA_character_
  diag_var_yhat <- NA_real_
  diag_var_anal2 <- NA_real_
  diag_grad_l1 <- NA_real_
  diag_sigma_min_eig <- NA_real_
  diag_B_min_eig <- NA_real_

  if (variance_method == "bootstrap") {
    user_args_internal <- user_args
# Construct engine args for replicate fits using the same settings
    engine_args <- list(
      standardize = standardize,
      trim_cap = trim_cap,
      on_failure = on_failure,
      auxiliary_means = auxiliary_means,
      control = control,
      n_total = NULL,
      start = start,
      family = family
    )
# Create a method-local bootstrap estimator closure to avoid non-agnostic hooks
    est_closure <- function(data, formula, engine_args, ...) {
      engine_args$variance_method <- "none"
      eng <- do.call(NMAR::el_engine, engine_args)
      NMAR::nmar(formula = formula, data = data, engine = eng)
    }
    boot_args <- list(
      data = full_data,
      estimator_func = est_closure,
      point_estimate = y_hat,
      bootstrap_reps = bootstrap_reps,
      formula = user_args_internal$formula %||% internal_formula$outcome,
      engine_args = engine_args
    )
# Execute bootstrap
    boot_try <- tryCatch(
      list(result = do.call(bootstrap_variance, boot_args), message = "Calculation successful"),
      error = function(e) list(result = NULL, message = paste("Bootstrap failed:", e$message))
    )
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
      nmar_scaling_recipe = nmar_scaling_recipe
    )
    se_y_hat <- var_out$se_y_hat
    vcov_unscaled <- var_out$vcov_unscaled
    vcov_message <- var_out$vcov_message
    A_condition <- var_out$diag$A_condition
    diag_grad_source <- var_out$diag$grad_source %||% NA_character_
    diag_var_yhat <- var_out$diag$var_y_hat_val %||% NA_real_
    diag_var_anal2 <- var_out$diag$var_anal2 %||% NA_real_
    diag_grad_l1 <- var_out$diag$grad_l1 %||% NA_real_
    diag_sigma_min_eig <- var_out$diag$sigma_min_eig %||% NA_real_
    diag_B_min_eig <- var_out$diag$B_min_eig %||% NA_real_
# Guard against pathological near-zero or non-finite delta variance
    if (!is.finite(se_y_hat) || se_y_hat < 1e-8) {
      warning("Delta variance returned non-finite or near-zero standard error; returning NA. Consider variance_method='bootstrap'.", call. = FALSE)
      se_y_hat <- NA_real_
      vcov_unscaled <- matrix(NA_real_, K_beta, K_beta, dimnames = list(colnames(response_model_matrix_unscaled), colnames(response_model_matrix_unscaled)))
      vcov_message <- "Delta variance near-zero or non-finite; returned NA."
    }
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
      solver_method = solver_method_used,
      nleqslv_global = nleqslv_global_used,
      nleqslv_xscalm = nleqslv_xscalm_used,
      solver_iterations = if (!is.null(solution$iter)) solution$iter else NA,
      solver_time = solver_time,
      variance_time = variance_time,
      reparam_W = "logit",
      max_equation_residual = max_eq_resid,
      jacobian_condition_number = A_condition,
      aux_inconsistency_max_z = aux_inconsistency_max_z,
      aux_inconsistency_cols = aux_inconsistency_cols,
      grad_source = diag_grad_source,
      var_y_hat_val = diag_var_yhat,
      var_anal2 = diag_var_anal2,
      grad_l1 = diag_grad_l1,
      sigma_min_eig = diag_sigma_min_eig,
      B_min_eig = diag_B_min_eig,
      min_denominator = denom_stats$min,
      fraction_small_denominators = denom_stats$p_small,
      denom_q01 = denom_q[[1]],
      denom_q05 = denom_q[[2]],
      denom_median = denom_q[[3]],
      denom_count_lt_1e4 = denom_cnt_1e4,
      weight_max_share = weight_max_share,
      weight_top5_share = weight_top5_share,
      weight_ess = weight_ess,
      constraint_sum_W = constraint_eqW_sum,
      constraint_sum_aux = constraint_aux_sum
    ),
    nmar_scaling_recipe = nmar_scaling_recipe, fitted_values = drop(w_i_hat)
  ))
}
