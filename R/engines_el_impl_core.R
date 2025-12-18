#' Core Empirical Likelihood Estimator
#'
#' Implements the core computational engine for empirical likelihood estimation
#' under nonignorable nonresponse, including parameter solving, variance calculation,
#' and diagnostic computation.
#'
#' @param missingness_design Respondent-side missingness (response) model design matrix (intercept + predictors).
#' @param aux_matrix Auxiliary design matrix on respondents (may have zero columns).
#' @param aux_means Named numeric vector of auxiliary population means (aligned to columns of \code{aux_matrix}).
#' @param auxiliary_means Named numeric vector of known population means supplied by the user (optional; used for diagnostics).
#' @param respondent_weights Numeric vector of respondent weights aligned with \code{missingness_design} rows.
#' @param analysis_data Data object used for logging and variance (survey designs supply the design object).
#' @param outcome_expr Character string identifying the outcome expression displayed in outputs.
#' @param N_pop Population size on the analysis scale.
#' @param formula Original model formula used for estimation.
#' @param standardize Logical. Whether to standardize predictors during estimation.
#' @param trim_cap Numeric. Upper bound for empirical likelihood weight trimming.
#' @param control List of control parameters for the nonlinear equation solver.
#' @param on_failure Character. Action when solver fails: "return" or "error".
#' @param family List. Link function specification (typically logit).
#' @param variance_method Character. Variance estimation method.
#' @param bootstrap_reps Integer. Number of bootstrap replications.
#'
#' @return List containing estimation results, diagnostics, and metadata.
#'
#' @details
#' Orchestrates EL estimation for NMAR following Qin, Leung, and Shao (2002).
#' For \code{data.frame} inputs (IID setting) the stacked system in
#' \eqn{(\beta, z, \lambda_x)} with \eqn{z = \logit(W)} is solved by
#' \code{nleqslv::nleqslv()} using an analytic Jacobian. For \code{survey.design} inputs a
#' design-weighted analogue in \eqn{(\beta, z, \lambda_W, \lambda_x)} is solved
#' with an analytic Jacobian when the response family supplies second
#' derivatives, or with numeric/Broyden Jacobians otherwise. Numerical
#' safeguards are applied consistently across equations, Jacobian, and
#' post-solution weights: bounded linear predictors, probability clipping in
#' ratios, and a small floor on denominators \eqn{D_i(\theta)} with an
#' active-set mask in derivatives. After solving, unnormalized masses
#' \eqn{d_i/D_i(\theta)} are formed, optional trimming may be applied (with
#' normalization only for reporting), and optional variance is computed via
#' bootstrap when \code{variance_method = "bootstrap"}.
#'
#' @keywords internal
el_estimator_core <- function(missingness_design,
                              aux_matrix, aux_means,
                              respondent_weights,
                              analysis_data,
                              outcome_expr,
                              N_pop, formula,
                              standardize, trim_cap, control,
                              on_failure, family = logit_family(),
                              variance_method, bootstrap_reps,
                              start = NULL, trace_level = 0,
                              auxiliary_means = NULL) {

# Setup
  force(family)
  if (!is.matrix(missingness_design)) {
    stop("Internal error: missingness_design must be a matrix.", call. = FALSE)
  }
  if (is.null(aux_matrix)) {
    aux_matrix <- matrix(nrow = nrow(missingness_design), ncol = 0)
  }
  if (is.null(aux_means)) {
    aux_means <- numeric(0)
  }
  has_aux <- ncol(aux_matrix) > 0

  verboser <- create_verboser(trace_level)
  el_log_banner(verboser, "EMPIRICAL LIKELIHOOD ESTIMATION STARTED")
  el_log_trace(verboser, trace_level)

# Unscaled matrices and outcome extraction
  missingness_model_matrix_unscaled <- missingness_design
  aux_matrix_unscaled <- aux_matrix
  mu_x_unscaled <- aux_means
  K_aux <- if (has_aux) ncol(aux_matrix_unscaled) else 0
# Outcome is carried as a dedicated predictor column in the missingness design.
  if (!is.null(outcome_expr) && outcome_expr %in% colnames(missingness_model_matrix_unscaled)) {
    response_outcome <- missingness_model_matrix_unscaled[, outcome_expr]
  } else {
    stop("Internal error: outcome column not found in missingness design.", call. = FALSE)
  }

  n_resp_weighted <- sum(respondent_weights)
  K_beta <- ncol(missingness_model_matrix_unscaled)
  el_log_data_prep(
    verboser = verboser,
    outcome_var = outcome_expr,
    family_name = family$name %||% "<unknown>",
    K_beta = K_beta,
    K_aux = K_aux,
    auxiliary_names = if (K_aux > 0) colnames(aux_matrix_unscaled) else character(0),
    standardize = standardize,
    is_survey = inherits(analysis_data, "survey.design"),
    N_pop = N_pop,
    n_resp_weighted = n_resp_weighted
  )

# Scaling / standardization
  scaling_result <- validate_and_apply_nmar_scaling(
    standardize = standardize,
    has_aux = has_aux,
    response_model_matrix_unscaled = missingness_model_matrix_unscaled,
    aux_matrix_unscaled = aux_matrix_unscaled,
    mu_x_unscaled = mu_x_unscaled,
    weights = respondent_weights
  )
  nmar_scaling_recipe <- scaling_result$nmar_scaling_recipe
  missingness_model_matrix_scaled <- scaling_result$response_model_matrix_scaled
  auxiliary_matrix_scaled <- scaling_result$auxiliary_matrix_scaled
  mu_x_scaled <- scaling_result$mu_x_scaled

  n_resp_weighted <- sum(respondent_weights)

# Heuristic inconsistency check for user-supplied auxiliary means
  auxiliary_inconsistency_max_z <- NA_real_
  auxiliary_inconsistency_cols <- character(0)
  if (has_aux && !is.null(auxiliary_means)) {
    thr <- getOption("nmar.el_aux_z_threshold", 8)
    if (!is.numeric(thr) || length(thr) != 1L || !is.finite(thr) || thr <= 0) thr <- 8
    chk <- el_check_auxiliary_inconsistency_matrix(aux_matrix_unscaled, provided_means = auxiliary_means)
    auxiliary_inconsistency_max_z <- chk$max_z
    auxiliary_inconsistency_cols <- chk$cols
    if (is.finite(auxiliary_inconsistency_max_z) && auxiliary_inconsistency_max_z > thr) {
      warning(sprintf(
        "Auxiliary means appear far from respondents' support (max |z| = %.2f, threshold = %.2f). Proceeding; see diagnostics.",
        auxiliary_inconsistency_max_z, thr
      ), call. = FALSE)
    }
  }
# Build equation system and analytic Jacobian (if available)
  is_survey_design <- inherits(analysis_data, "survey.design")
  if (is_survey_design) {
    equation_system_func <- el_build_equation_system_survey(
      family = family,
      missingness_model_matrix = missingness_model_matrix_scaled,
      auxiliary_matrix = auxiliary_matrix_scaled,
      respondent_weights = respondent_weights,
      N_pop = N_pop,
      n_resp_weighted = n_resp_weighted,
      mu_x_scaled = mu_x_scaled
    )
    analytical_jac_func <- el_build_jacobian_survey(
      family = family,
      missingness_model_matrix = missingness_model_matrix_scaled,
      auxiliary_matrix = auxiliary_matrix_scaled,
      respondent_weights = respondent_weights,
      N_pop = N_pop,
      n_resp_weighted = n_resp_weighted,
      mu_x_scaled = mu_x_scaled
    )
  } else {
    equation_system_func <- el_build_equation_system(
      family = family,
      missingness_model_matrix = missingness_model_matrix_scaled,
      auxiliary_matrix = auxiliary_matrix_scaled,
      respondent_weights = respondent_weights,
      N_pop = N_pop,
      n_resp_weighted = n_resp_weighted,
      mu_x_scaled = mu_x_scaled
    )
    analytical_jac_func <- el_build_jacobian(
      family = family,
      missingness_model_matrix = missingness_model_matrix_scaled,
      auxiliary_matrix = auxiliary_matrix_scaled,
      respondent_weights = respondent_weights,
      N_pop = N_pop,
      n_resp_weighted = n_resp_weighted,
      mu_x_scaled = mu_x_scaled
    )
  }
# Solve the stacked estimating equations
  K_beta <- ncol(missingness_model_matrix_scaled)
  K_aux <- ncol(auxiliary_matrix_scaled)

# Starting values: (beta, z, lambda_x); lambda_W is explicit for survey designs.
  st <- el_build_start(
    missingness_model_matrix_scaled = missingness_model_matrix_scaled,
    auxiliary_matrix_scaled = auxiliary_matrix_scaled,
    nmar_scaling_recipe = nmar_scaling_recipe,
    start = start,
    N_pop = N_pop,
    respondent_weights = respondent_weights
  )
  init_beta <- st$init_beta
  init_z <- st$init_z
  init_lambda <- st$init_lambda
  if (is_survey_design) {
    init_lambda_W <- 0
    init <- c(unname(init_beta), init_z, init_lambda_W, unname(init_lambda))
  } else {
    init <- st$init
  }

  prep <- el_prepare_nleqslv(control)
  control_top <- prep$top
  final_control <- prep$control

  el_log_solver_config(verboser, control_top, final_control)
  el_log_start_values(verboser, init_beta, init_z, init_lambda)

  t_solve_start <- proc.time()[[3]]

  el_log_solving(verboser)
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
          auxiliary_inconsistency_max_z = auxiliary_inconsistency_max_z,
          auxiliary_inconsistency_cols = auxiliary_inconsistency_cols
        ),
        nmar_scaling_recipe = nmar_scaling_recipe
      ))
    }
  }

# Post-processing and point estimate
  estimates <- solution$x
  beta_hat_scaled <- estimates[1:K_beta]
  if (is_survey_design) {
    z_idx <- K_beta + 1L
    lambda_W_idx <- K_beta + 2L
    lambda_start_idx <- K_beta + 3L
    lambda_hat <- if (K_aux > 0) estimates[lambda_start_idx:(lambda_start_idx + K_aux - 1L)] else numeric(0)
    lambda_W_solved <- estimates[lambda_W_idx]
  } else {
    z_idx <- K_beta + 1L
    lambda_hat <- if (K_aux > 0) estimates[(K_beta + 2):(K_beta + 1 + K_aux)] else numeric(0)
    lambda_W_solved <- NULL
  }
  solver_time <- proc.time()[[3]] - t_solve_start
  Xc_centered_diag <- NULL
  if (K_aux > 0) {
    mu_match <- as.numeric(mu_x_scaled[colnames(auxiliary_matrix_scaled)])
    Xc_centered_diag <- sweep(auxiliary_matrix_scaled, 2, mu_match, "-")
  }

  post <- el_post_solution(
    estimates = estimates,
    missingness_model_matrix_scaled = missingness_model_matrix_scaled,
    missingness_model_matrix_unscaled = missingness_model_matrix_unscaled,
    auxiliary_matrix_scaled = auxiliary_matrix_scaled,
    mu_x_scaled = mu_x_scaled,
    response_outcome = response_outcome,
    family = family,
    N_pop = N_pop,
    respondent_weights = respondent_weights,
    K_beta = K_beta,
    K_aux = K_aux,
    nmar_scaling_recipe = nmar_scaling_recipe,
    standardize = standardize,
    trim_cap = trim_cap,
    X_centered = Xc_centered_diag,
    lambda_W = lambda_W_solved
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

  el_log_weight_diagnostics(verboser, W_hat, w_unnorm_trimmed, post$trimmed_fraction)
  el_log_detailed_diagnostics(verboser, beta_hat_unscaled, W_hat, lambda_W_hat, lambda_hat, denominator_hat)

# Diagnostics at the solution
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
  constraint_link_sum <- diag_pack$constraint_sum_link
  sum_respondent_weights <- diag_pack$sum_respondent_weights
  sum_unnormalized_weights_untrimmed <- diag_pack$sum_unnormalized_weights_untrimmed
  normalization_ratio <- diag_pack$normalization_ratio

# Check if trimming forced a loss of total mass (cap too tight)
  if (is.finite(post$trimmed_fraction) && post$trimmed_fraction > 0) {
    sum_w_trimmed <- sum(w_unnorm_trimmed)
    sum_w_untrimmed <- sum_unnormalized_weights_untrimmed
    trimming_mass_shift <- if (is.finite(sum_w_untrimmed) && sum_w_untrimmed > 0) {
      abs(sum_w_trimmed / sum_w_untrimmed - 1)
    } else {
      NA_real_
    }

    if (is.finite(trimming_mass_shift) && trimming_mass_shift > 0.05) {
      warning(
        sprintf(
          "Trimming changed the total unnormalized EL mass by %.1f%% (>5%% threshold).\n",
          trimming_mass_shift * 100
        ),
        "This typically indicates trim_cap is too tight to preserve total mass.\n",
        "Bootstrap variance is recommended when trimming is active.",
        call. = FALSE
      )
    }
  }

  denom_q <- diag_pack$denom_q
  denom_cnt_1e4 <- diag_pack$denom_cnt_1e4
  weight_max_share <- diag_pack$weight_max_share
  weight_top5_share <- diag_pack$weight_top5_share
  weight_ess <- diag_pack$weight_ess

# Variance calculation (bootstrap or none)
  se_y_hat <- NA
  vcov_unscaled <- NULL
  vcov_message <- "Calculation successful"

  el_log_variance_header(verboser, variance_method, bootstrap_reps)

  t_var_start <- proc.time()[[3]]

  var_out <- el_compute_variance(
    y_hat = y_hat,
    full_data = analysis_data,
    formula = formula,
    N_pop = N_pop,
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
  variance_time <- proc.time()[[3]] - t_var_start
  el_log_variance_result(verboser, se_y_hat, variance_time)

  el_log_final(verboser, y_hat, se_y_hat)

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
      auxiliary_inconsistency_max_z = auxiliary_inconsistency_max_z,
      auxiliary_inconsistency_cols = auxiliary_inconsistency_cols,
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
      constraint_sum_link = constraint_link_sum,
      sum_respondent_weights = sum_respondent_weights,
      sum_unnormalized_weights_untrimmed = sum_unnormalized_weights_untrimmed,
      normalization_ratio = normalization_ratio,
      max_constraint_residual = max(
        abs(constraint_eqW_sum),
        if (length(constraint_aux_sum) > 0) max(abs(constraint_aux_sum)) else 0,
        if (is.finite(constraint_link_sum)) abs(constraint_link_sum) else 0,
        na.rm = TRUE
      )
    ),
    nmar_scaling_recipe = nmar_scaling_recipe, fitted_values = drop(w_i_hat)
  ))
}
