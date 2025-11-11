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
#' Orchestrates EL estimation for NMAR following Qin, Leung, and Shao (2002).
#' The stacked system in \eqn{(\beta, z, \lambda_x)} with \eqn{z = \logit(W)} is
#' solved by \code{nleqslv} using an analytic Jacobian. Numerical safeguards are
#' applied consistently across equations, Jacobian, and post-solution weights:
#' bounded linear predictors, probability clipping in ratios, and a small floor
#' on denominators \code{D_i(\theta)} with an active-set mask in derivatives.
#' After solving, unnormalized masses \code{d_i/D_i(\theta)} are formed, optional
#' trimming may be applied (with normalization only for reporting), and optional
#' variance is computed via bootstrap. The analytical delta method for EL is not
#' implemented.
#'
#' @keywords internal
el_estimator_core <- function(full_data, respondent_data, respondent_weights, N_pop,
                              internal_formula, auxiliary_means, standardize,
                              trim_cap, control,
                              on_failure, family = logit_family(),
                              variance_method, bootstrap_reps,
                              user_args, start = NULL, trace_level = 0, ...) {

# 0. Setup
  force(family)
  outcome_var <- all.vars(internal_formula$outcome)[1]

# Create verboser for verbose output
  verboser <- create_verboser(trace_level)
# Start banner and trace message
  el_log_banner(verboser, "EMPIRICAL LIKELIHOOD ESTIMATION STARTED")
  el_log_trace(verboser, trace_level)

# 1. Data Preparation
  has_aux <- !is.null(internal_formula$auxiliary)

  response_model_formula <- update(internal_formula$response, NULL ~ .)
  response_model_matrix_unscaled <- model.matrix(response_model_formula, data = respondent_data, na.action = stats::na.pass)

# Resolve auxiliaries and their population means in a single place
  aux_res <- el_resolve_auxiliaries(
    full_data = full_data,
    respondent_data = respondent_data,
    aux_formula = internal_formula$auxiliary,
    auxiliary_means = auxiliary_means
  )
  auxiliary_matrix_unscaled <- aux_res$matrix
  mu_x_unscaled <- aux_res$means
  has_aux <- isTRUE(aux_res$has_aux)

# Data summary
  n_resp_weighted <- sum(respondent_weights)
  K_beta <- ncol(response_model_matrix_unscaled)
  K_aux <- if (has_aux) ncol(auxiliary_matrix_unscaled) else 0
  el_log_data_prep(
    verboser = verboser,
    outcome_var = outcome_var,
    family_name = family$name %||% "<unknown>",
    K_beta = K_beta,
    K_aux = K_aux,
    aux_names = if (K_aux > 0) colnames(auxiliary_matrix_unscaled) else character(0),
    standardize = standardize,
    is_survey = inherits(full_data, "survey.design"),
    N_pop = N_pop,
    n_resp_weighted = n_resp_weighted
  )

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

# Heuristic inconsistency check for user-supplied auxiliary means
  aux_inconsistency_max_z <- NA_real_
  aux_inconsistency_cols <- character(0)
  if (has_aux && !is.null(auxiliary_means)) {
    thr <- getOption("nmar.el_aux_z_threshold", 8)
    if (!is.numeric(thr) || length(thr) != 1L || !is.finite(thr) || thr <= 0) thr <- 8
    chk <- el_check_aux_inconsistency(respondent_data, internal_formula$auxiliary, provided_means = auxiliary_means, threshold = thr)
    aux_inconsistency_max_z <- chk$max_z
    aux_inconsistency_cols <- chk$cols
    if (is.finite(aux_inconsistency_max_z) && aux_inconsistency_max_z > thr) {
      warning(sprintf(
        "Auxiliary means appear far from respondents' support (max |z| = %.2f, threshold = %.2f). Proceeding; see diagnostics.",
        aux_inconsistency_max_z, thr
      ), call. = FALSE)
    }
  }
# Build full stacked builders
  equation_system_func <- el_build_equation_system(
    family = family, response_model_matrix = response_model_matrix_scaled, auxiliary_matrix = auxiliary_matrix_scaled,
    respondent_weights = respondent_weights, N_pop = N_pop, n_resp_weighted = n_resp_weighted, mu_x_scaled = mu_x_scaled
  )
  analytical_jac_func <- el_build_jacobian(
    family = family, response_model_matrix = response_model_matrix_scaled, auxiliary_matrix = auxiliary_matrix_scaled,
    respondent_weights = respondent_weights, N_pop = N_pop, n_resp_weighted = n_resp_weighted, mu_x_scaled = mu_x_scaled
  )
# 4. Solve for Parameters with a safeguarded Newton method
  K_beta <- ncol(response_model_matrix_scaled)
  K_aux <- ncol(auxiliary_matrix_scaled)

# Build starts using a single helper
  st <- el_build_start(
    response_model_matrix_scaled = response_model_matrix_scaled,
    auxiliary_matrix_scaled = auxiliary_matrix_scaled,
    nmar_scaling_recipe = nmar_scaling_recipe,
    start = start,
    N_pop = N_pop,
    respondent_weights = respondent_weights
  )
  init_beta <- st$init_beta
  init_z <- st$init_z
  init_lambda <- st$init_lambda
  init <- st$init

# Prepare nleqslv args and control
  prep <- el_prepare_nleqslv(control)
  control_top <- prep$top
  final_control <- prep$control

# Solver configuration and starts
  el_log_solver_config(verboser, control_top, final_control)
  el_log_start_values(verboser, init_beta, init_z, init_lambda)

# Instrumentation for timing and solver used
  t_solve_start <- proc.time()[[3]]

# Solving status message
  el_log_solving(verboser)
  {

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
    N_pop = N_pop,
    trace_level = trace_level
  )
  solution <- solver_out$solution
  solver_method_used <- solver_out$method
  nleqslv_global_used <- if (!is.null(solver_out$used_top$global)) solver_out$used_top$global else NA_character_
  nleqslv_xscalm_used <- if (!is.null(solver_out$used_top$xscalm)) solver_out$used_top$xscalm else NA_character_

# Convergence report
  converged_success <- !(any(is.na(solution$x)) || solution$termcd > 2)
  el_log_solver_result(verboser, converged_success, solution, proc.time()[[3]] - t_solve_start)

  if (any(is.na(solution$x)) || solution$termcd > 2) {
    if (on_failure == "error") {
      stop(
        "Empirical likelihood solver failed to converge: ", solution$message,
        "\n  Try increasing iterations with control = list(maxit = 500) or use on_failure = 'return' for diagnostics.",
        call. = FALSE
      )
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
# Precompute centered auxiliary design once for reuse
  Xc_centered_diag <- NULL
  if (K_aux > 0) {
    mu_match <- as.numeric(mu_x_scaled[colnames(auxiliary_matrix_scaled)])
    Xc_centered_diag <- sweep(auxiliary_matrix_scaled, 2, mu_match, "-")
  }

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
    trim_cap = trim_cap,
    X_centered = Xc_centered_diag
  )
  if (post$error) {
    if (on_failure == "error") {
      stop(post$message, call. = FALSE)
    } else {
      return(list(
        converged = FALSE,
        diagnostics = list(convergence_code = -1, message = post$message),
        nmar_scaling_recipe = nmar_scaling_recipe
      ))
    }
  }
  y_hat <- post$y_hat
  w_unnorm_trimmed <- post$weights # Unnormalized (trimmed) EL masses: d_i/D_i
  beta_hat_unscaled <- post$beta_hat_unscaled
  W_hat <- post$W_hat
  lambda_hat <- post$lambda_hat
  eta_i_hat <- post$eta_i_hat
  w_i_hat <- post$w_i_hat
  denominator_hat <- post$denominator_hat
  lambda_W_hat <- post$lambda_W_hat

# Weight and detail diagnostics
  el_log_weight_diagnostics(verboser, W_hat, w_unnorm_trimmed, post$trimmed_fraction)
  el_log_detailed_diagnostics(verboser, beta_hat_unscaled, W_hat, lambda_W_hat, lambda_hat, denominator_hat)

# 6. Diagnostics at Solution
  diag_pack <- el_compute_diagnostics(
    estimates = estimates,
    equation_system_func = equation_system_func,
    analytical_jac_func = analytical_jac_func,
    post = post,
    respondent_weights = respondent_weights,
    auxiliary_matrix_scaled = auxiliary_matrix_scaled,
    K_beta = K_beta,
    K_aux = K_aux,
    X_centered = Xc_centered_diag
  )
  constraint_eqW_sum <- diag_pack$constraint_sum_W
  constraint_aux_sum <- diag_pack$constraint_sum_aux
  sum_respondent_weights <- diag_pack$sum_respondent_weights
  sum_unnormalized_weights_untrimmed <- diag_pack$sum_unnormalized_weights_untrimmed
  normalization_ratio <- diag_pack$normalization_ratio

# Check if trimming has materially altered the mass distribution
  if (is.finite(post$trimmed_fraction) && post$trimmed_fraction > 0) {
    sum_w_trimmed <- sum(w_unnorm_trimmed)
    trimming_mass_shift <- abs(sum_w_trimmed / sum_respondent_weights - 1)

    if (trimming_mass_shift > 0.05) {
      warning(
        sprintf(
          "Trimming altered the EL mass by %.1f%% (>5%% threshold).\n",
          trimming_mass_shift * 100
        ),
        "Constraints remain solved but the identity sum(d_i/D_i)=sum(d_i) no longer holds exactly for trimmed weights.\n",
        "Bootstrap variance is recommended when trimming is active.",
        call. = FALSE
      )
    }
  }

# Additional diagnostics
  denom_q <- diag_pack$denom_q
  denom_cnt_1e4 <- diag_pack$denom_cnt_1e4
  weight_max_share <- diag_pack$weight_max_share
  weight_top5_share <- diag_pack$weight_top5_share
  weight_ess <- diag_pack$weight_ess

# 7. Conditional Variance Calculation
  se_y_hat <- NA
  vcov_unscaled <- NA
  vcov_message <- "Calculation successful"

# Variance estimation header
  el_log_variance_header(verboser, variance_method, bootstrap_reps)

  t_var_start <- proc.time()[[3]]
  diag_grad_source <- NA_character_
  diag_var_yhat <- NA_real_
  diag_var_anal2 <- NA_real_
  diag_grad_l1 <- NA_real_
  diag_sigma_min_eig <- NA_real_
  diag_B_min_eig <- NA_real_
  vcov_unscaled <- NA
  vcov_message <- "Calculation successful"

  var_out <- el_compute_variance(
    y_hat = y_hat,
    full_data = full_data,
    internal_formula = internal_formula,
    user_args = user_args,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    standardize = standardize,
    trim_cap = trim_cap,
    on_failure = on_failure,
    auxiliary_means = auxiliary_means,
    control = control,
    start = start,
    family = family
  )
  vcov_message <- var_out$message
  se_y_hat <- var_out$se
  if (identical(variance_method, "none")) {
    vcov_unscaled <- matrix(NA_real_, K_beta, K_beta, dimnames = list(colnames(response_model_matrix_unscaled), colnames(response_model_matrix_unscaled)))
  }
  variance_time <- proc.time()[[3]] - t_var_start
  el_log_variance_result(verboser, se_y_hat, variance_time)

# Final banner and summary
  el_log_final(verboser, y_hat, se_y_hat)

# 8. Return Final Results List
  return(list(
    y_hat = y_hat, se = se_y_hat, weights = w_unnorm_trimmed,
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
      max_equation_residual = diag_pack$max_equation_residual,
      jacobian_condition_number = diag_pack$jacobian_condition_number,
      aux_inconsistency_max_z = aux_inconsistency_max_z,
      aux_inconsistency_cols = aux_inconsistency_cols,
      grad_source = diag_grad_source,
      var_y_hat_val = diag_var_yhat,
      var_anal2 = diag_var_anal2,
      grad_l1 = diag_grad_l1,
      sigma_min_eig = diag_sigma_min_eig,
      B_min_eig = diag_B_min_eig,
      min_denominator = diag_pack$denom_stats$min,
      fraction_small_denominators = diag_pack$denom_stats$p_small,
      denom_q01 = denom_q[[1]],
      denom_q05 = denom_q[[2]],
      denom_median = denom_q[[3]],
      denom_count_lt_1e4 = denom_cnt_1e4,
      denom_floor = diag_pack$denom_floor,
      denom_floor_hits = diag_pack$denom_stats$p_floor,
      weight_max_share = weight_max_share,
      weight_top5_share = weight_top5_share,
      weight_ess = weight_ess,
      constraint_sum_W = constraint_eqW_sum,
      constraint_sum_aux = constraint_aux_sum,
      sum_respondent_weights = sum_respondent_weights,
      sum_unnormalized_weights_untrimmed = sum_unnormalized_weights_untrimmed,
      normalization_ratio = normalization_ratio,
      max_constraint_residual = max(abs(constraint_eqW_sum), if (length(constraint_aux_sum) > 0) max(abs(constraint_aux_sum)) else 0, na.rm = TRUE)
    ),
    nmar_scaling_recipe = nmar_scaling_recipe, fitted_values = drop(w_i_hat)
  ))
}
