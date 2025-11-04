#' @importFrom nleqslv nleqslv
#' @importFrom stats as.formula coef dnorm dgamma sd setNames
#' @exportS3Method exptilt data.frame
exptilt.data.frame <- function(data, formula,
                               auxiliary_means = NULL,
                               standardize = TRUE,
                               prob_model_type = c("logit", "probit"),
                               y_dens = c("auto", "normal", "lognormal", "exponential"),
                               variance_method = c("delta", "bootstrap"),
                               bootstrap_reps = 10,
                               control = list(),
                               stopping_threshold = 1,
                               on_failure = c("return", "error"),
                               supress_warnings = FALSE,
                               design_weights = NULL,
                               survey_design = NULL,
                               trace_level = 0,
                               sample_size = 2000,
                               ...) {
  prob_model_type <- match.arg(prob_model_type)
  y_dens <- match.arg(y_dens)
  variance_method <- match.arg(variance_method)
  on_failure <- match.arg(on_failure)

  outcome_var <- all.vars(formula[[2]])[1]
# Split RHS into auxiliaries and response predictors via `|` if present
  rhs <- formula[[3]]
  aux_expr <- rhs
  resp_expr <- NULL
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    aux_expr <- rhs[[2L]]
    resp_expr <- rhs[[3L]]
  }
  aux_vars <- unique(all.vars(aux_expr))
  if (length(aux_vars) == 0) aux_vars <- character()
  resp_pred <- if (is.null(resp_expr)) character() else unique(all.vars(resp_expr))
# Response-only predictors; the outcome enters the response model automatically elsewhere
  response_predictors <- unique(resp_pred)

  required_cols <- unique(c(outcome_var, aux_vars, response_predictors))
  data_subset <- data[, required_cols, drop = FALSE]

  if (is.null(survey_design)) {
    data <- as.data.frame(data)
    data_subset <- as.data.frame(data_subset)
  }

# Stratified sampling for memory optimization
# Preserve respondent/non-respondent ratio
  n_total <- nrow(data_subset)
  respondent_idx <- which(!is.na(data_subset[, outcome_var]))
  nonrespondent_idx <- which(is.na(data_subset[, outcome_var]))
  n_resp <- length(respondent_idx)
  n_nonresp <- length(nonrespondent_idx)

  sampling_performed <- FALSE
  original_n_total <- n_total
  original_n_resp <- n_resp
  original_n_nonresp <- n_nonresp

  if (n_total > sample_size) {
# Calculate stratified sample sizes
    resp_ratio <- n_resp / n_total
    nonresp_ratio <- n_nonresp / n_total

    n_resp_sample <- round(sample_size * resp_ratio)
    n_nonresp_sample <- sample_size - n_resp_sample

# Ensure we have at least some observations from each stratum
    n_resp_sample <- max(1, min(n_resp_sample, n_resp))
    n_nonresp_sample <- max(1, min(n_nonresp_sample, n_nonresp))

# Sample from each stratum
    sampled_resp_idx <- sample(respondent_idx, n_resp_sample, replace = FALSE)
    sampled_nonresp_idx <- sample(nonrespondent_idx, n_nonresp_sample, replace = FALSE)

# Combine sampled indices
    sampled_idx <- c(sampled_resp_idx, sampled_nonresp_idx)

# Update data_subset
    data_subset <- data_subset[sampled_idx, , drop = FALSE]

    sampling_performed <- TRUE
  }

  model <- list(
    data = data_subset,
    required_cols = required_cols,
    col_y = outcome_var,
    cols_y_observed = aux_vars,
    cols_delta = response_predictors,
    prob_model_type = prob_model_type,
    y_dens = y_dens,
    stopping_threshold = stopping_threshold,
    auxiliary_means = auxiliary_means,
    standardize = standardize,
    control = control,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,

    supress_warnings = supress_warnings,
    design_weights = design_weights,
    design = survey_design,
    formula = formula,
    call = match.call(),
    trace_level = trace_level,

# Sampling information
    sample_size = sample_size,
    sampling_performed = sampling_performed,
    original_n_total = original_n_total,
    original_n_resp = original_n_resp,
    original_n_nonresp = original_n_nonresp
  )

# Create verboser
  model$verboser <- create_verboser(trace_level)

  class(model) <- "nmar_exptilt"

  model$is_survey <- !is.null(survey_design)
  model$family <- if (prob_model_type == "logit") {
    logit_family()
  } else {
    probit_family()
  }


  model$original_params <- model # for bootstrap purposes, to re-run exptilt_fit_model

  exptilt_fit_model(data_subset, model, on_failure = on_failure, ...)
}

exptilt_fit_model <- function(data, model, on_failure = c("return", "error"), ...) {

  on_failure <- match.arg(on_failure)
  model$x <- data
  model$is_survey <- isTRUE(model$is_survey)
  if (is.null(model$standardize)) {
    model$standardize <- TRUE
  }

# Initialize verboser if not already present
  if (is.null(model$verboser)) {
    model$verboser <- create_verboser(trace_level = 0)
  }

  verboser <- model$verboser

  verboser("============================================================", level = 1, type = "step")
  verboser("  EXPTILT ESTIMATION STARTED", level = 1, type = "step")
  verboser("============================================================", level = 1, type = "step")

# Show trace level info
  trace_msg <- sprintf("Running with trace_level = %d", model$trace_level)
  if (model$trace_level < 3) {
    trace_msg <- paste0(trace_msg, sprintf(" | For more detailed output, use trace_level = %d", model$trace_level + 1), ". Available trace_level = c(1,2,3)")
  }
  verboser(trace_msg, level = 1)

# Show formula at level 1
  formula_str <- deparse(model$formula)
  verboser(sprintf("Formula: %s", formula_str), level = 1)

  model$cols_required <- colnames(model$x)

# model$x_1 <- model$x[!is.na(model$x[,model$col_y]),,drop=FALSE] #observed
# model$x_0 <- model$x[is.na(model$x[,model$col_y]),,drop=FALSE] #unobserved
# model$y_1 <- model$x_1[,model$col_y,drop=TRUE] #observed y
#
# model$x_for_y_obs <- model$x_1[,model$cols_y_observed,drop=FALSE]
# model$x_for_y_unobs <- model$x_0[,model$cols_y_observed,drop=FALSE]

  has_aux <- length(model$cols_y_observed) > 0 && !is.null(model$auxiliary_means)
  filtered_aux_means <- if (has_aux) {
    aux_names <- intersect(names(model$auxiliary_means), model$cols_y_observed)
    model$auxiliary_means[aux_names]
  } else {
    NULL
  }

  respondent_mask <- !is.na(model$x[, model$col_y])
  model$respondent_mask <- respondent_mask

# Level 1: Show basic data info
  n_total <- length(respondent_mask)
  n_resp <- sum(respondent_mask)
  n_nonresp <- n_total - n_resp
  pct_nonresp <- 100 * n_nonresp / n_total

  verboser("", level = 1)
  verboser("-- DATA SUMMARY --", level = 1)
  verboser(sprintf("  Total observations:      %d", n_total), level = 1)
  verboser(sprintf("  Respondents:             %d (%.1f%%)", n_resp, 100 - pct_nonresp), level = 1)
  verboser(sprintf("  Non-respondents:         %d (%.1f%%)", n_nonresp, pct_nonresp), level = 1)

# Display sampling information if sampling was performed
  if (isTRUE(model$sampling_performed)) {
    verboser("", level = 1)
    verboser("-- STRATIFIED SAMPLING --", level = 1)
    verboser(sprintf("  Original sample size:    %d", model$original_n_total), level = 1)
    verboser(sprintf("    Respondents:           %d (%.1f%%)",
                     model$original_n_resp,
                     100 * model$original_n_resp / model$original_n_total), level = 1)
    verboser(sprintf("    Non-respondents:       %d (%.1f%%)",
                     model$original_n_nonresp,
                     100 * model$original_n_nonresp / model$original_n_total), level = 1)
    verboser(sprintf("  Sampled size:            %d (sample_size=%d)", n_total, model$sample_size), level = 1)
    verboser(sprintf("    Respondents:           %d (%.1f%%)", n_resp, 100 - pct_nonresp), level = 1)
    verboser(sprintf("    Non-respondents:       %d (%.1f%%)", n_nonresp, pct_nonresp), level = 1)
    verboser("  Ratio preserved:         stratified sampling", level = 1)
  }

  verboser("", level = 2)
  verboser("-- MODEL SPECIFICATION --", level = 2)
  verboser(sprintf("  Outcome variable:        %s", model$col_y), level = 2)

  if (length(model$cols_y_observed) > 0) {
    aux_str <- paste(model$cols_y_observed, collapse = ", ")
    verboser(sprintf("  Auxiliary variables:     %s", aux_str), level = 2)
  } else {
    verboser("  Auxiliary variables:     (none)", level = 2)
  }

  if (length(model$cols_delta) > 0) {
    miss_str <- paste(model$cols_delta, collapse = ", ")
    verboser(sprintf("  Missingness predictors:  %s", miss_str), level = 2)
  } else {
    verboser("  Missingness predictors:  (intercept only)", level = 2)
  }

  verboser(sprintf("  Response model:          %s", model$prob_model_type), level = 2)
  verboser(sprintf("  Outcome density:         %s", model$y_dens), level = 2)
  verboser(sprintf("  Standardization:         %s", if (model$standardize) "enabled" else "disabled"), level = 2)

  scaling_weights <- model$design_weights
# When scaling we only want respondent rows to contribute; pass a mask so the
# helper can zero-out nonrespondents without reallocating weights
  weight_mask <- if (length(respondent_mask) == nrow(model$x)) respondent_mask else NULL

  scaling_result <- validate_and_apply_nmar_scaling(
    standardize = model$standardize,
    has_aux = has_aux,
    response_model_matrix_unscaled = model$x[, c(model$col_y, model$cols_delta), drop = FALSE],
    auxiliary_matrix_unscaled = model$x[, model$cols_y_observed, drop = FALSE],
    mu_x_unscaled = filtered_aux_means,
    weights = scaling_weights,
    weight_mask = weight_mask
  )

  if (model$standardize) {
    verboser("  OK Variables standardized", level = 2)
  }

  model$nmar_scaling_recipe <- scaling_result$nmar_scaling_recipe
  response_model_matrix_scaled <- scaling_result$response_model_matrix_scaled
  auxiliary_matrix_scaled <- scaling_result$auxiliary_matrix_scaled
  mu_x_scaled <- scaling_result$mu_x_scaled

  model$x_1 <- response_model_matrix_scaled[respondent_mask, , drop = FALSE] # observed
  model$x_0 <- response_model_matrix_scaled[!respondent_mask, , drop = FALSE] # unobserved
  model$y_1 <- if (nrow(model$x_1)) model$x_1[, model$col_y, drop = TRUE] else numeric(0) # observed y
  model$x_for_y_obs <- auxiliary_matrix_scaled[respondent_mask, , drop = FALSE]
  model$x_for_y_unobs <- auxiliary_matrix_scaled[!respondent_mask, , drop = FALSE]

# Normalize design weights to respondent length
# - For non-survey: unit weights for respondents
# - For survey: subset analysis weights to respondent rows to avoid length mismatches
  if (!isTRUE(model$is_survey)) {
    model$design_weights <- rep(1, nrow(model$x_1))
  } else {
    if (!is.null(model$design_weights)) {
      if (length(model$design_weights) == nrow(model$x)) {
        model$design_weights <- model$design_weights[respondent_mask]
      } else if (length(model$design_weights) != nrow(model$x_1)) {
        stop(sprintf(
          "design_weights length (%d) must match respondents (%d) or full data (%d).",
          length(model$design_weights), nrow(model$x_1), nrow(model$x)
        ), call. = FALSE)
      }
    } else {
      model$design_weights <- rep(1, nrow(model$x_1))
    }
  }


# Track the current scale of feature matrices. We fit f1(.) on the scaled
# space and compute EM steps there. After unscaling coefficients for
# presentation we flip this flag to FALSE so downstream density evaluations
# can re-apply the same scaling recipe to inputs as needed
  model$features_are_scaled <- TRUE


# Smarter initialization: Use data-driven heuristics
# (GLM doesn't work well because unobserved Y values are unknown)

model$theta <- stats::runif(length(model$cols_delta) + 2, -0.1, 0.1)
names(model$theta) <- c("(Intercept)", model$cols_delta, model$col_y)

  verboser("", level = 2)
  verboser("-- PARAMETER INITIALIZATION --", level = 2)
  verboser(sprintf("  Number of parameters:    %d", length(model$theta)), level = 2)

# Pretty print theta with meaningful names at level 3
  if (model$verboser("", level = 3, type = "detail") != invisible(NULL) || TRUE) {
    verboser("  Initial values:", level = 3)
    theta_names <- names(model$theta)
    theta_labels <- character(length(theta_names))
    for (i in seq_along(theta_names)) {
      if (theta_names[i] == "(Intercept)") {
        theta_labels[i] <- sprintf("    %-25s = %7.4f  [response model intercept]",
                                   theta_names[i], model$theta[i])
      } else if (theta_names[i] == model$col_y) {
        theta_labels[i] <- sprintf("    %-25s = %7.4f  [outcome effect on response]",
                                   theta_names[i], model$theta[i])
      } else {
        theta_labels[i] <- sprintf("    %-25s = %7.4f  [missingness predictor]",
                                   theta_names[i], model$theta[i])
      }
    }
    for (label in theta_labels) {
      verboser(label, level = 3)
    }
  }

  verboser("", level = 1)
  verboser("-- CONDITIONAL DENSITY ESTIMATION --", level = 1)
  dens_response <- generate_conditional_density(model)
  model$density_fun <- dens_response$density_function
  model$density_num_of_coefs <- dens_response$num_of_coefs
  model$chosen_y_dens <- dens_response$chosen_distribution

  verboser(sprintf("  Selected distribution:   %s", model$chosen_y_dens), level = 1)
  verboser(sprintf("  Density parameters:      %d", model$density_num_of_coefs), level = 2)

# Show density model details at level 3
  if (!is.null(dens_response$density_model)) {
    verboser("  Fitted density model summary:", level = 3)
    verboser("", level = 3, obj = summary(dens_response$density_model))
  }

  if (!is.null(dens_response$aic_comparison) && length(dens_response$aic_comparison) > 1) {
    verboser("  AIC comparison (lower is better):", level = 3)
    for (dist_name in names(dens_response$aic_comparison)) {
      aic_val <- dens_response$aic_comparison[dist_name]
      marker <- if (dist_name == model$chosen_y_dens) " [selected]" else ""
      verboser(sprintf("    %s: %.2f%s", dist_name, aic_val, marker), level = 3)
    }
  }

# model$O_matrix_nieobs <- generate_Odds(model, model$theta)

  verboser("  Computing density matrices...", level = 2)
  model$f_matrix_nieobs <- generate_conditional_density_matrix(model)
  model$C_matrix_nieobs <- generate_C_matrix(model)

  return(exptilt_estimator_core(
    model = model,
    respondent_mask = respondent_mask,
    on_failure = on_failure,
    ...
  ))
}

#' @keywords internal
exptilt_estimator_core <- function(model, respondent_mask,
                                   on_failure = "return", ...) {
  model$cols_required <- colnames(model$x)
  verboser <- model$verboser

# Optimized solver: Let nleqslv do its own iteration instead of manual loop
# Add early stopping based on score magnitude
  early_stop_threshold <- model$stopping_threshold

  verboser("", level = 1)
  verboser("-- NONLINEAR SOLVER (nleqslv) --", level = 1)
  verboser(sprintf("  Early stopping threshold: %.4f", early_stop_threshold), level = 2)

  target_function <- function(theta) {
# model$theta <<- theta
    O_matrix_nieobs_current <- generate_Odds(model, theta)
    result <- step_func(model, theta, O_matrix_nieobs_current)

# Early stopping: if score is very small, return zero to signal convergence
    if (max(abs(result)) < early_stop_threshold) {
      return(rep(0, length(result)))
    }
    return(result)
  }

# Use nleqslv with user-provided control parameters
  nleqslv_args <- list(
    x = model$theta,
    fn = target_function
  )

# Add method and global if provided in control
  if (!is.null(model$control$method)) {
    nleqslv_args$method <- model$control$method
  }
  if (!is.null(model$control$global)) {
    nleqslv_args$global <- model$control$global
  }

# Add other control parameters if any are provided
  control_params <- model$control[!names(model$control) %in% c("method", "global")]
  if (length(control_params) > 0) {
    nleqslv_args$control <- control_params
  }

  if (!is.null(nleqslv_args$method)) {
    verboser(sprintf("  Method:                   %s", nleqslv_args$method), level = 2)
  }
  if (!is.null(nleqslv_args$global)) {
    verboser(sprintf("  Global strategy:          %s", nleqslv_args$global), level = 2)
  }
  if (length(control_params) > 0) {
    verboser("  Control parameters:", level = 3)
    for (param_name in names(control_params)) {
      verboser(sprintf("    %s = %s", param_name, control_params[[param_name]]), level = 3)
    }
  }

  verboser("  Solving...", level = 1)
  solution <- do.call(nleqslv, nleqslv_args)

  model$theta <- solution$x
  model$loss_value <- solution$fvec
  model$iterations <- solution$iter

# Convergence status messages
  conv_status <- if (solution$termcd <= 2) "OK Converged" else "Warning: Convergence issue"
  verboser("", level = 1)
  verboser(sprintf("  %s", conv_status), level = 1)
  verboser(sprintf("  Iterations:               %d", solution$iter), level = 1)
  verboser(sprintf("  Termination code:         %d", solution$termcd), level = 2)
  verboser(sprintf("  Max |score|:              %.6f", max(abs(solution$fvec))), level = 2)

# Show final theta with labels at level 3
  verboser("  Final parameter estimates (scaled):", level = 3)
  for (i in seq_along(model$theta)) {
    param_name <- names(model$theta)[i]
    if (param_name == "(Intercept)") {
      label <- "[response intercept]"
    } else if (param_name == model$col_y) {
      label <- "[outcome -> response]"
    } else {
      label <- "[missingness pred.]"
    }
    verboser(sprintf("    %-20s = %9.6f  %s", param_name, model$theta[i], label), level = 3)
  }

# Check convergence
  if (solution$termcd > 2) {
    cat("Warning: nleqslv termcd =", solution$termcd, "| max|score| =", max(abs(solution$fvec)), "\n")
  }

  if (model$standardize) {
    verboser("  Unscaling coefficients to original scale...", level = 2)
    unscale <- unscale_coefficients(model$theta, matrix(0, length(model$theta), length(model$theta)), model$nmar_scaling_recipe)
    model$theta <- unscale$coefficients

    verboser("  Final parameter estimates (original scale):", level = 3)
    for (i in seq_along(model$theta)) {
      param_name <- names(model$theta)[i]
      if (param_name == "(Intercept)") {
        label <- "[response intercept]"
      } else if (param_name == model$col_y) {
        label <- "[outcome -> response]"
      } else {
        label <- "[missingness pred.]"
      }
      verboser(sprintf("    %-20s = %9.6f  %s", param_name, model$theta[i], label), level = 3)
    }
  }

  model$x_1 <- model$x[respondent_mask, , drop = FALSE]
  model$x_0 <- model$x[!respondent_mask, , drop = FALSE]
  model$y_1 <- if (nrow(model$x_1)) model$x_1[, model$col_y, drop = TRUE] else numeric(0)

  model$x_for_y_obs <- model$x_1[, model$cols_y_observed, drop = FALSE]
  model$x_for_y_unobs <- model$x_0[, model$cols_y_observed, drop = FALSE]
# From this point the feature matrices are on the original (unscaled) space
# density helpers will internally re-apply the scaling recipe captured during
# model fitting so that gamma_hat remains consistent
  model$features_are_scaled <- FALSE

# model$O_matrix_nieobs <- generate_Odds(model, model$theta)
# model$f_matrix_nieobs <- generate_conditional_density_matrix(model)
# model$C_matrix_nieobs <- generate_C_matrix(model)

## VARIANCE LOGIC - Cleaned version

  default_vcov <- matrix(NA_real_, nrow = length(model$theta), ncol = length(model$theta))
  var_results <- list(var_est = NA_real_, vcov = default_vcov)
  se_final <- NaN

# Initialize bootstrap flag if not present
  if (is.null(model$is_bootstrap_running)) {
    model$is_bootstrap_running <- FALSE
  }

# Determine which variance method to use
  use_delta_method <- function(model) {
# If we're already in a bootstrap run (exptilt_fit_model, which is current function is run to estimate bootstrap variance), never use delta
    if (isTRUE(model$is_bootstrap_running)) {
      return(FALSE)
    }

# Check if delta method is requested and applicable
    if (!identical(model$variance_method, "delta")) {
      return(FALSE)
    }

    if (!identical(model$chosen_y_dens, "normal")) {
      warning("Delta variance unavailable for y_dens='", model$chosen_y_dens,
              "'; falling back to bootstrap.", call. = FALSE)
      return(FALSE)
    }

# Check sample size and weight conditions

    n_resp <- as.numeric(length(model$design_weights))
    design_varies <- (n_resp > 1) && (max(abs(model$design_weights - model$design_weights[1])) > 1e-6)
    few_resp <- n_resp < 40
    if (design_varies || few_resp) {
      warning("Delta variance may be unreliable with the current sample; using bootstrap instead.",
              call. = FALSE)
      return(FALSE)
    }

    if (isTRUE(model$is_survey)) {
      warning("Delta variance is not supported for survey designs; using bootstrap instead.",
              call. = FALSE)
      return(FALSE)
    }

    return(TRUE)
  }

# Try delta method if applicable
  if (use_delta_method(model)) {
    verboser("", level = 1)
    verboser("-- VARIANCE ESTIMATION (Delta Method) --", level = 1)
    delta_attempt <- try(estim_var(model), silent = TRUE)

    if (inherits(delta_attempt, "try-error")) {
      error_msg <- conditionMessage(attr(delta_attempt, "condition"))
      verboser("  Warning: Delta method failed, switching to bootstrap", level = 2)
      warning("Delta variance failed to evaluate; using bootstrap instead.\n",
              call. = FALSE)
    } else {
      var_results <- delta_attempt
      se_final <- sqrt(var_results$var_est)
      verboser("  OK Delta method complete", level = 1)
      verboser(sprintf("  Standard error:           %.6f", se_final), level = 2)
    }
  }

# Use bootstrap if delta wasn't used or failed
  use_bootstrap <- is.nan(se_final) && !isTRUE(model$is_bootstrap_running)

  if (use_bootstrap) {
    verboser("", level = 1)
    verboser("-- VARIANCE ESTIMATION (Bootstrap) --", level = 1)
    verboser(sprintf("  Bootstrap replications:   %d", model$bootstrap_reps), level = 1)
    bootstrap_runner <- function(data, ...) {
# Create a copy for bootstrap with delta method and bootstrap flag
      bootstrap_model <- model$original_params
      bootstrap_model$variance_method <- "delta"
      bootstrap_model$is_bootstrap_running <- TRUE
      bootstrap_model$verboser <- function(...) invisible(NULL)
# bootstrap_model$verbose <- FALSE
      bootstrap_model$trace_level = 0

# If 'data' is a survey.design, extract variables and weights to reuse the
# unified data.frame path (avoids recursion and enforces consistent internals)
      if (inherits(data, "survey.design")) {
        design_vars <- data$variables
        design_weights <- as.numeric(stats::weights(data))
        bootstrap_model$is_survey <- TRUE
        bootstrap_model$design <- data
        bootstrap_model$design_weights <- design_weights
        data <- design_vars
      } else {
        bootstrap_model$is_survey <- FALSE
      }

# Call without passing on_failure explicitly to avoid duplicate argument issues
      call_args <- list(
        data = data,
        model = bootstrap_model
      )

# Add any additional arguments except on_failure to avoid conflicts
      dots <- list(...)
      dots$on_failure <- NULL
      call_args <- c(call_args, dots)

      do.call(exptilt_fit_model, call_args)
    }

    base_args <- list(
      data = if (model$is_survey) model$design else model$data,
      estimator_func = bootstrap_runner,
      point_estimate = estim_mean(model),
      bootstrap_reps = model$bootstrap_reps
    )

# Use the default survey_na_policy from bootstrap_variance() (strict).

    if (!model$is_survey) {
      respondent_mask_guard <- respondent_mask
      base_args$resample_guard <- function(indices, data) {
        any(respondent_mask_guard[indices])
      }
    }

    bootstrap_results <- do.call(bootstrap_variance, base_args)
    se_final <- bootstrap_results$se
    verboser("  OK Bootstrap complete", level = 1)
    verboser(sprintf("  Standard error:           %.6f", se_final), level = 1)
  } else if (isTRUE(model$is_bootstrap_running)) {
# If we're in bootstrap mode but not actually running bootstrap (delta succeeded)
    se_final <- sqrt(var_results$var_est)
  }

# Final result assembly
  mean_estimate <- estim_mean(model)

  verboser("", level = 1)
  verboser("============================================================", level = 1, type = "result")
  verboser("  ESTIMATION COMPLETE", level = 1, type = "result")
  verboser("============================================================", level = 1, type = "result")
  verboser(sprintf("  Mean estimate:            %.6f", mean_estimate), level = 1, type = "result")
  verboser(sprintf("  Standard error:           %.6f", se_final), level = 1, type = "result")
  verboser(sprintf("  95%% CI:                   [%.6f, %.6f]",
                   mean_estimate - 1.96 * se_final,
                   mean_estimate + 1.96 * se_final), level = 1, type = "result")

  verboser("", level = 2)
  verboser("  Response model coefficients:", level = 2)
  for (i in seq_along(model$theta)) {
    param_name <- names(model$theta)[i]
    if (param_name == "(Intercept)") {
      desc <- "Intercept"
    } else if (param_name == model$col_y) {
      desc <- sprintf("Effect of %s on response prob.", model$col_y)
    } else {
      desc <- sprintf("Effect of %s on response prob.", param_name)
    }
    verboser(sprintf("    %-20s: %9.6f  (%s)", param_name, model$theta[i], desc), level = 2)
  }
  verboser("============================================================", level = 1, type = "result")

  result <- new_nmar_result_exptilt(
    estimate = mean_estimate,
    se = se_final,
    coefficients = model$theta,
    vcov = if (use_bootstrap) default_vcov else var_results$vcov,
    model = model,
    converged = TRUE,
    weights = model$design_weights,
    variance_message = NA_character_
  )

  return(validate_nmar_result(result, "nmar_result_exptilt"))
}
