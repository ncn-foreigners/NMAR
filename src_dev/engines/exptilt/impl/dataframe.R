#' @importFrom nleqslv nleqslv
#' @importFrom stats as.formula coef dnorm dgamma sd setNames
#' @exportS3Method exptilt data.frame
exptilt.data.frame <- function(data, formula,
                               standardize = TRUE,
                               prob_model_type = c("logit", "probit"),
                               y_dens = c("normal", "lognormal", "exponential", "binomial"),
                               variance_method = c("bootstrap", 'none'),
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
  is_survey <- !is.null(survey_design)


  res = et_extract_formula(format(formula), data)
  Y = res$Y
  X = res$X
  Z = res$Z
  outcome_var <- as.vector(colnames(Y)[1])
  if (is.null(Z)) {
    Z = matrix(1, nrow = nrow(X), ncol = 1)
    colnames(Z) = "(Intercept)"
  }

# browser()
  if (length(colnames(Y)) > 1) {
    stop("Exptilt supports only single outcome variable.")
  }

  X = X[, !colnames(X) %in% c("(Intercept)"), drop = FALSE]

  if (outcome_var %in% colnames(Z)) {
    warning(
      sprintf(
        "Outcome variable (%s) found in missingness predictors; Performance with / without %s on the right side is the same \n",
        outcome_var,
        outcome_var
      ),
      call. = FALSE
    )
  }

  Z = Z[, !colnames(Z) %in% c("(Intercept)", outcome_var), drop = FALSE]

  cols_y_observed <- as.vector(colnames(X))
  cols_delta <- as.vector(colnames(Z))

  et_validate_df(X, Y, Z)


  data <- cbind(
    as.data.frame(Y),
    as.data.frame(X),
    as.data.frame(Z)
  )

  n_total <- nrow(data)
  respondent_idx <- which(!is.na(data[, outcome_var]))
  nonrespondent_idx <- which(is.na(data[, outcome_var]))
  n_resp <- length(respondent_idx)
  n_nonresp <- length(nonrespondent_idx)

  sampling_performed <- FALSE
# sampled_idx initially entire dataset
  sampled_idx <- seq_len(n_total)
  original_n_total <- n_total
  original_n_resp <- n_resp
  original_n_nonresp <- n_nonresp
  if (!is_survey) {
    design_weights <- rep(1, nrow(data))
  }

  if (n_total > sample_size) {
# if(!is_survey){
    resp_ratio <- n_resp / n_total
    nonresp_ratio <- n_nonresp / n_total

    n_resp_sample <- round(sample_size * resp_ratio)
    n_nonresp_sample <- sample_size - n_resp_sample

    n_resp_sample <- max(1, min(n_resp_sample, n_resp))
    n_nonresp_sample <- max(1, min(n_nonresp_sample, n_nonresp))

    sampled_resp_idx <- sample(respondent_idx, n_resp_sample, replace = FALSE)
    sampled_nonresp_idx <- sample(nonrespondent_idx, n_nonresp_sample, replace = FALSE)

    sampled_idx <- c(sampled_resp_idx, sampled_nonresp_idx)

    data <- data[sampled_idx, , drop = FALSE]
    design_weights <- design_weights[sampled_idx]

    sampling_performed <- TRUE
# }
  }


  respondent_mask <- !is.na(data[, outcome_var])

  model <- list(
    data = data,
# required_cols = required_cols,
    col_y = outcome_var,
    cols_y_observed = cols_y_observed,
    cols_delta = cols_delta,
    prob_model_type = prob_model_type,
    y_dens = y_dens,
    stopping_threshold = stopping_threshold,
# auxiliary_means = auxiliary_means,
    standardize = standardize,
    control = control,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    is_survey = is_survey,
    supress_warnings = supress_warnings,
    design_weights = design_weights,
    design = survey_design,
    formula = formula,
    trace_level = trace_level,
    on_failure = on_failure,
    call = match.call(),
# respondent_mask = respondent_mask,
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
  if (model$is_survey) {
    tmp <- survey_design[sampled_idx, ]
    tmp$variables <- data
    data <- tmp

  }


  model$original_params <- model # for bootstrap purposes, to re-run exptilt_fit_model

# browser()
  exptilt_fit_model(data, model, ...)
}

exptilt_fit_model <- function(data, model, ...) {
# browser()

  model$data <- if (model$is_survey) data$variables else data
  model$respondent_mask <- !is.na(model$data[, model$col_y])
  model$design_weights <- if (model$is_survey) {
    as.numeric(stats::weights(data))
  } else {
    model$design_weights # 1 1 1 1 1
  }

  verboser <- model$verboser
# browser()
  verboser("============================================================", level = 1, type = "step")
  verboser("  EXPTILT ESTIMATION STARTED", level = 1, type = "step")
  verboser("============================================================", level = 1, type = "step")


  trace_msg <- sprintf("Running with trace_level = %d", model$trace_level)
  if (model$trace_level < 3) {
    trace_msg <- paste0(trace_msg, sprintf(" | For more detailed output, use trace_level = %d", model$trace_level + 1), ". Available trace_level = c(1,2,3)")
  }
  verboser(trace_msg, level = 1)


  formula_str <- deparse(model$formula)
  verboser(sprintf("Formula: %s", formula_str), level = 1)

  pct_nonresp <- 100 * model$original_n_resp / model$original_n_total

  verboser("", level = 1)
  verboser("-- DATA SUMMARY --", level = 1)
  verboser(sprintf("  Total observations:      %d", model$original_n_total), level = 1)
  verboser(sprintf("  Respondents:             %d (%.1f%%)", model$original_n_resp, 100 - pct_nonresp), level = 1)
  verboser(sprintf("  Non-respondents:         %d (%.1f%%)", model$original_n_nonresp, pct_nonresp), level = 1)

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
    verboser(sprintf("  Sampled size:            %d (sample_size=%d)", model$original_n_total, model$sample_size), level = 1)
    verboser(sprintf("    Respondents:           %d (%.1f%%)", model$original_n_resp, 100 - pct_nonresp), level = 1)
    verboser(sprintf("    Non-respondents:       %d (%.1f%%)", model$original_n_nonresp, pct_nonresp), level = 1)
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

# When scaling we only want respondent rows to contribute; pass a mask so the
# helper can zero-out nonrespondents without reallocating weights

  scaling_result <- validate_and_apply_nmar_scaling(
    standardize = model$standardize,
    has_aux = F,
    response_model_matrix_unscaled = model$data[, c(model$col_y, model$cols_delta), drop = FALSE],
    aux_matrix_unscaled = model$data[, model$cols_y_observed, drop = FALSE],
    mu_x_unscaled = NULL,
    weights = model$design_weights,
    weight_mask = model$respondent_mask
  )

  if (model$standardize) {
    verboser("  OK Variables standardized", level = 2)
  }

  model$nmar_scaling_recipe <- scaling_result$nmar_scaling_recipe
  response_model_matrix_scaled <- scaling_result$response_model_matrix_scaled
  auxiliary_matrix_scaled <- scaling_result$auxiliary_matrix_scaled
  mu_x_scaled <- scaling_result$mu_x_scaled

  model$data_1 <- response_model_matrix_scaled[model$respondent_mask, , drop = FALSE] # observed
  model$data_0 <- response_model_matrix_scaled[!model$respondent_mask, , drop = FALSE] # unobserved
  model$y_1 <- if (nrow(model$data_1)) model$data_1[, model$col_y, drop = TRUE] else numeric(0) # observed y
  model$data_for_y_obs <- auxiliary_matrix_scaled[model$respondent_mask, , drop = FALSE]
  model$data_for_y_unobs <- auxiliary_matrix_scaled[!model$respondent_mask, , drop = FALSE]


  model$features_are_scaled <- TRUE


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
    data = data,
    model = model,
    ...
  ))
}

#' @keywords internal
exptilt_estimator_core <- function(data, model, ...) {
  model$cols_required <- colnames(model$data)
  verboser <- model$verboser

# Optimized solver: Let nleqslv do its own iteration instead of manual loop
# Add early stopping based on score magnitude
  early_stop_threshold <- model$stopping_threshold

  verboser("", level = 1)
  verboser("-- NONLINEAR SOLVER (nleqslv) --", level = 1)
  verboser(sprintf("  Early stopping threshold: %.4f", early_stop_threshold), level = 2)
# browser()
  target_function <- function(theta) {
# model$theta <<- theta
    O_matrix_nieobs_current <- generate_Odds(model, theta)
    result <- step_func(model, theta, O_matrix_nieobs_current)

# Early stopping: if score is very small, return zero to signal convergence
# browser()
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
# browser()
  solution <- do.call(nleqslv, nleqslv_args)
# browser()

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

  model$data_1 <- model$data[model$respondent_mask, , drop = FALSE]
  model$data_0 <- model$data[!model$respondent_mask, , drop = FALSE]
  model$y_1 <- if (nrow(model$data_1)) model$data_1[, model$col_y, drop = TRUE] else numeric(0)

  model$data_for_y_obs <- model$data_1[, model$cols_y_observed, drop = FALSE]
  model$data_for_y_unobs <- model$data_0[, model$cols_y_observed, drop = FALSE]
# From this point the feature matrices are on the original (unscaled) space
# density helpers will internally re-apply the scaling recipe captured during
# model fitting so that gamma_hat remains consistent
  model$features_are_scaled <- FALSE


# default_vcov <- matrix(NA_real_, nrow = length(model$theta), ncol = length(model$theta))
  var_results <- list(var_est = NA_real_
# , vcov = default_vcov
                      )
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
  use_bootstrap <- is.nan(se_final) && !isTRUE(model$is_bootstrap_running) && model$variance_method != "none"
# browser()
  if (use_bootstrap) {
    verboser("", level = 1)
    verboser("-- VARIANCE ESTIMATION (Bootstrap) --", level = 1)
    verboser(sprintf("  Bootstrap replications:   %d", model$bootstrap_reps), level = 1)
    bootstrap_runner <- function(data, ...) {
# Create a copy for bootstrap with delta method and bootstrap flag
      bootstrap_model <- model$original_params
      bootstrap_model$variance_method <- "none"
      bootstrap_model$is_bootstrap_running <- TRUE
      bootstrap_model$verboser <- function(...) invisible(NULL)
# bootstrap_model$verbose <- FALSE
      bootstrap_model$trace_level = 0

# If 'data' is a survey.design, extract variables and weights to reuse the
# unified data.frame path (avoids recursion and enforces consistent internals)
# if(model$is_survey) {
#   # design_vars <- data$variables
#   # design_weights <- as.numeric(stats::weights(data))
#   # bootstrap_model$is_survey <- TRUE
#   # bootstrap_model$design <- data
#   # bootstrap_model$design_weights <- design_weights
#   # data <- design_vars
# } else {
#   bootstrap_model$is_survey <- FALSE
# }

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
# data = if (model$is_survey) model$design else model$data,
# data=data,
      data = data,
      estimator_func = bootstrap_runner,
      point_estimate = estim_mean(model),
      bootstrap_reps = model$bootstrap_reps
    )

# Use the default survey_na_policy from bootstrap_variance() (strict).

    if (!model$is_survey) {
      base_args$resample_guard <- function(indices, data) {
        any(model$respondent_mask[indices])
      }
    }


    bootstrap_results <- do.call(bootstrap_variance, base_args)
    se_final <- bootstrap_results$se
# browser()
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
  if (model$variance_method != 'none') {
  verboser(sprintf("  Standard error:           %.6f", se_final), level = 1, type = "result")
  verboser(sprintf("  95%% CI:                   [%.6f, %.6f]",
                   mean_estimate - 1.96 * se_final,
                   mean_estimate + 1.96 * se_final), level = 1, type = "result")
  }
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
# vcov = if (use_bootstrap) default_vcov else var_results$vcov,
    model = model,
    converged = TRUE,
    weights = model$design_weights,
    variance_message = NA_character_
  )

  return(validate_nmar_result(result, "nmar_result_exptilt"))
}
