#' @importFrom nleqslv nleqslv
#' @importFrom stats as.formula coef dnorm dgamma sd setNames
#' @exportS3Method exptilt data.frame
exptilt.data.frame <- function(data, formula, response_predictors = NULL,
                               auxiliary_means = NULL,
                               standardize = TRUE,
                               prob_model_type = c("logit", "probit"),
                               y_dens = c("auto", "normal", "lognormal", "exponential"),
                               variance_method = c("delta", "bootstrap"),
                               bootstrap_reps = 10,
                               min_iter = 10,
                               max_iter = 100,
                               tol_value = 1e-5,
                               optim_method = c("Newton", "Broyden"),
                               on_failure = c("return", "error"),
                               supress_warnings = FALSE,
                               design_weights = NULL,
                               survey_design = NULL,
                               ...) {
  prob_model_type <- match.arg(prob_model_type)
  y_dens <- match.arg(y_dens)
  variance_method <- match.arg(variance_method)
  optim_method <- match.arg(optim_method)
  on_failure <- match.arg(on_failure)

  outcome_var <- all.vars(formula[[2]])[1]
  aux_vars <- unique(all.vars(formula[[3]]))
  if (length(aux_vars) == 0) aux_vars <- character()

  if (is.null(response_predictors)) response_predictors <- character()
  response_predictors <- unique(response_predictors)

  required_cols <- unique(c(outcome_var, aux_vars, response_predictors))
  data_subset <- data[, required_cols, drop = FALSE]

  model <- list(
    data = data,
    col_y = outcome_var,
    cols_y_observed = aux_vars,
    cols_delta = response_predictors,
    prob_model_type = prob_model_type,
    y_dens = y_dens,
    tol_value = tol_value,
    min_iter = min_iter,
    auxiliary_means = auxiliary_means,
    standardize = standardize,
    max_iter = max_iter,
    optim_method = optim_method,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    supress_warnings = supress_warnings,
    design_weights = design_weights,
    design = survey_design,
    formula = formula,
    call = match.call()
  )

  class(model) <- "nmar_exptilt"

  model$is_survey <- !is.null(survey_design)
  model$family <- if (prob_model_type == "logit") {
    logit_family()
  } else {
    probit_family()
  }

  model$original_params <- unserialize(serialize(model, NULL))

  exptilt_fit_model(data_subset, model, on_failure = on_failure, ...)
}

exptilt_fit_model <- function(data, model, on_failure = c("return", "error"), ...) {
  on_failure <- match.arg(on_failure)
  model$x <- data
  model$is_survey <- isTRUE(model$is_survey)
  if (is.null(model$standardize)) {
    model$standardize <- TRUE
  }
# Cache the pre-fit template so bootstrap replicates can start from the same
# initial state (rather than inheriting mutated fields from the first fit)
  bootstrap_template <- unserialize(serialize(model, NULL))
  model$cols_required <- colnames(model$x)
  bootstrap_template$cols_required <- model$cols_required
# model$x_1 <- model$x[!is.na(model$x[,model$col_y]),,drop=FALSE] #observed
# model$x_0 <- model$x[is.na(model$x[,model$col_y]),,drop=FALSE] #unobserved
# model$y_1 <- model$x_1[,model$col_y,drop=TRUE] #observed y
#
# model$x_for_y_obs <- model$x_1[,model$cols_y_observed,drop=FALSE]
# model$x_for_y_unobs <- model$x_0[,model$cols_y_observed,drop=FALSE]

  has_aux <- length(model$cols_y_observed) > 0 && !is.null(model$auxiliary_means)
# Retain auxiliary means only for RHS predictors present in this fit, this
# keeps scaling consistent when users supply a superset of moments
  filtered_aux_means <- if (has_aux) {
    aux_names <- intersect(names(model$auxiliary_means), model$cols_y_observed)
    model$auxiliary_means[aux_names]
  } else {
    NULL
  }

# Store design weights so downstream components (score, variance, C-matrix)
# all see the same vector. For plain data frames, fall back to weight 1 if
# the caller has not already supplied model$design_weights (the survey
# method pre-populates this vector)
# if (is.null(model$design_weights) || length(model$design_weights) != nrow(model$x)) {
#   model$design_weights <- rep(1, nrow(model$x))
# }

  respondent_mask <- !is.na(model$x[, model$col_y])
  model$respondent_mask <- respondent_mask
# browser()

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
  model$nmar_scaling_recipe <- scaling_result$nmar_scaling_recipe
  response_model_matrix_scaled <- scaling_result$response_model_matrix_scaled
  auxiliary_matrix_scaled <- scaling_result$auxiliary_matrix_scaled
  mu_x_scaled <- scaling_result$mu_x_scaled

  model$x_1 <- response_model_matrix_scaled[respondent_mask, , drop = FALSE] # observed
  model$x_0 <- response_model_matrix_scaled[!respondent_mask, , drop = FALSE] # unobserved
  model$y_1 <- if (nrow(model$x_1)) model$x_1[, model$col_y, drop = TRUE] else numeric(0) # observed y
  model$x_for_y_obs <- auxiliary_matrix_scaled[respondent_mask, , drop = FALSE]
  model$x_for_y_unobs <- auxiliary_matrix_scaled[!respondent_mask, , drop = FALSE]

  if (model$is_survey == F) model$respondent_weights <- rep(1, nrow(model$x_1)) else model$respondent_weights <- weights(model$x_1)


# Track the current scale of feature matrices. We fit f1(.) on the scaled
# space and compute EM steps there. After unscaling coefficients for
# presentation we flip this flag to FALSE so downstream density evaluations
# can re-apply the same scaling recipe to inputs as needed
  model$features_are_scaled <- TRUE


# We will derive respondent/nonrespondent weight slices on the fly from
# design_weights using the stored masks to avoid duplicating state

  # model$theta <- stats::runif(length(model$cols_delta) + 2, 0, 1)
  model$theta <- c(0.66,-0.36)
  # reg_model <- glm('Y ~ .', data = model$x[!is.na(model$x$Y), ], family = gaussian)
  # # summary(reg_model)
  # model$theta <- coef(reg_model)
  #everything apart of 1st index should be 0
  # model$theta[-c(1, length(model$theta))] <- 0
  # model$theta[1]<- 0
# Name the parameter vector to align with the design vector used throughout
# the ET implementation: [ (Intercept), x1 (cols_delta...), y ].  These names
# are required by the shared scaling unscaler so coefficients are returned on
# the original (unscaled) feature space after solving
  names(model$theta) <- c("(Intercept)", model$cols_delta, model$col_y)

  dens_response <- generate_conditional_density(model)

  model$density_fun <- dens_response$density_function
  model$density_fun_gradient <- dens_response$density_function_grad
  model$density_fun_hess <- dens_response$density_function_hess
  model$density_num_of_coefs <- dens_response$num_of_coefs
  model$chosen_y_dens <- dens_response$chosen_distribution
  model$O_matrix_nieobs <- generate_Odds(model,model$theta)

# const
  model$f_matrix_nieobs <- generate_conditional_density_matrix(model)
  model$C_matrix_nieobs <- generate_C_matrix(model)

  return(exptilt_estimator_core(
    model = model,
    bootstrap_template = bootstrap_template,
    respondent_mask = respondent_mask,
    on_failure = on_failure,
    ...
  ))
}

#' @keywords internal
exptilt_estimator_core <- function(model, bootstrap_template, respondent_mask,
                                   on_failure = "return", ...) {
  model$cols_required <- colnames(model$x)
  bootstrap_template$cols_required <- model$cols_required

  target_function <- function(theta) {
    # theta=c(0.67,-0.32)
    # theta=c(0.8,-0.1)
    model$theta <<- theta
    O_matrix_nieobs_current <- generate_Odds(model,theta)
    step_func(model, theta, O_matrix_nieobs_current)
  }

  solution <- nleqslv(
    x = model$theta,
    fn = target_function,
    method = "Newton",
    # jacobian = T,
    control = list(maxit = 1)
  )

  theta_prev <- model$theta
  model$theta <- solution$x
  model$loss_value <- solution$fvec
  iter <- 1

  while (sum((model$theta - theta_prev)^2) > model$tol_value &&
    (model$min_iter <=iter && iter <=model$max_iter)) {
    solution <- nleqslv(
      x = model$theta,
      fn = target_function,
      method = "Newton",
      # jacobian = T,
      control = list(maxit = 1)
    )

    theta_prev <- model$theta
    model$theta <- solution$x
    model$loss_value <- solution$fvec
    model$O_matrix_nieobs <- generate_Odds(model,model$theta)
    iter <- iter + 1
  }
  model$iterations <- iter

  if (model$standardize) {
    unscale <- unscale_coefficients(model$theta, matrix(0, length(model$theta), length(model$theta)), model$nmar_scaling_recipe)
    model$theta <- unscale$coefficients
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

  model$O_matrix_nieobs <- generate_Odds(model,model$theta)
  model$f_matrix_nieobs <- generate_conditional_density_matrix(model)
  model$C_matrix_nieobs <- generate_C_matrix(model)

  default_vcov <- matrix(NA_real_, nrow = length(model$theta), ncol = length(model$theta))
  var_results <- list(var_est = NA_real_, vcov = default_vcov)
  se_final <- NaN

  want_delta <- identical(model$variance_method, "delta") && identical(model$chosen_y_dens, "normal")
  use_bootstrap <- !want_delta

  if (want_delta) {
    delta_attempt <- try(estim_var(model), silent = TRUE)
    if (inherits(delta_attempt, "try-error")) {
      warning("Delta variance failed to evaluate; using bootstrap instead.", call. = FALSE)
      use_bootstrap <- TRUE
    } else {
      var_results <- delta_attempt
    }
  }

# If heuristics later steer us back to the bootstrap path, discard the delta
# variance so downstream code reports NA rather than an unreliable number
  if (isTRUE(model$is_survey) && !use_bootstrap) {
# The current implementation of analytic delta variance assumes IID sampling: it plugs
# respondent/nonrespondent scores into the Fisher blocks (F11, F21, F22) and
# combines them with the linearized estimating system S2(phi) to get var(tau_hat).
# Under complex designs the weights themselves change across replicates and
# the design induces correlation within PSUs/strata, so the IID "closed-form
# Sigma" (i.e., B = Var(S) where S = sum_i s_i is computed as a simple crossproduct/covariance)
# is no longer valid.
#
# To obtain a design-consistent analytic variance one of two routes is needed:
#  1) Replicate-weight linearization: for each set of replicate weights,
#     recompute the ET quantities that depend on weights (S2(phi), F11, F21,
#     F22, fractional weights w_ij, and any density-derivative sums) and then
#     pass replicate estimates (or score totals) to survey::withReplicates()
#     and survey::svrVar() to aggregate the design variance with the proper
#     scaling. The "score contrasts" here are the per-replicate
#     deviations of the totals/estimates from the full-sample target used by
#     svrVar to form Sigma.
#  2) Analytic influence-function path: derive per-unit linearized variables
#     (analytic influence functions) for phi_hat and tau_hat under ET, evaluated at the
#     full-sample fit, and hand them to survey::svyrecvar() (or the design's
#     replicate engine) to compute the design-based variance without refitting
#     replicates. This still requires careful propagation of respondent and
#     fractional-imputation weights through the ET linearization.
# Either approach would touch estim_var() (replace the IID "closed-form Sigma"
# with replicate or linearized aggregation) and likely add a small adapter akin
# to our bootstrap helper to drive replication. Until such a derivation is
# implemented, we always use bootstrap variance for survey designs to avoid
# under- or over-stating uncertainty. See Remark 2 form the Riddles paper for the
# weighting subtleties that must be preserved under replication.
    use_bootstrap <- TRUE
  }
  if (identical(model$variance_method, "delta") && !identical(model$chosen_y_dens, "normal")) {
    warning("Delta variance unavailable for y_dens='", model$chosen_y_dens, "'; falling back to bootstrap.", call. = FALSE)
  }
  if (identical(model$variance_method, "delta") && identical(model$chosen_y_dens, "normal")) {
# Heuristic guard: the linearization in Riddles et al. assumes a reasonably
# large respondent pool and light weighting. If the sample has few
# respondents or highly uneven weights, steering to the bootstrap tends to
# give more stable uncertainty estimates

    resp_w_local <- model$design_weights[respondent_mask]
    n_resp <- length(resp_w_local)
    design_varies <- n_resp > 1 && max(abs(resp_w_local - resp_w_local[1])) > 1e-6
    few_resp <- n_resp < 40
    if (design_varies || few_resp) {
      warning("Delta variance may be unreliable with the current sample; using bootstrap instead.", call. = FALSE)
      use_bootstrap <- TRUE
    }
  }

  if (use_bootstrap) {
    var_results <- list(var_est = NA_real_, vcov = default_vcov)
  }

  if (use_bootstrap) {
    bootstrap_runner <- function(data, ...) {
      template <- unserialize(serialize(bootstrap_template, NULL))
      if (inherits(data, "survey.design")) {
        data_df <- data$variables[, template$cols_required, drop = FALSE]
        template$design <- data
        template$is_survey <- TRUE
        template$design_weights <- as.numeric(stats::weights(data))
        exptilt_fit_model(data_df, template, on_failure = on_failure, ...)
      } else {
        data_df <- data[, template$cols_required, drop = FALSE]
        template$design <- NULL
        template$is_survey <- FALSE
        template$design_weights <- if (is.null(template$design_weights) ||
          length(template$design_weights) != nrow(data_df)) {
          rep(1, nrow(data_df))
        } else {
          template$design_weights
        }
        exptilt_fit_model(data_df, template, on_failure = on_failure, ...)
      }
    }

    base_args <- list(
      data = if (model$is_survey) model$design else model$x,
      estimator_func = bootstrap_runner,
      point_estimate = estim_mean(model),
      bootstrap_reps = model$bootstrap_reps
    )
    if (!model$is_survey) {
# Guard: resample until a respondent is present so EM/variance code can
# evaluate each replicate. Fallback to NA (with warning) if this fails
      respondent_mask_guard <- respondent_mask
      base_args$resample_guard <- function(indices, data) {
        any(respondent_mask_guard[indices])
      }
    }

    bootstrap_results <- do.call(bootstrap_variance, base_args)
    se_final <- bootstrap_results$se
  } else {
    se_final <- sqrt(var_results$var_est)
  }

  result <- new_nmar_result_exptilt(
    estimate = estim_mean(model),
    se = se_final,
    coefficients = model$theta,
    vcov = var_results$vcov,
    model = model,
    converged = TRUE,
    weights = model$design_weights[respondent_mask],
    variance_message = NA_character_
  )

  return(validate_nmar_result(result, "nmar_result_exptilt"))
}
