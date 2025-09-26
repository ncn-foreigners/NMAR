#' @importFrom nleqslv nleqslv
#' @importFrom stats as.formula coef dnorm dgamma sd setNames
#' @export
exptilt.data.frame <- function(data,model,on_failure=c('return'), ...){ #' #todo on_failure logic
  model$x=data
  model$is_survey <- isTRUE(model$is_survey)
  if (is.null(model$standardize)) {
    model$standardize <- TRUE
  }
  # Cache the pre-fit template so bootstrap replicates can start from the same
  # initial state (rather than inheriting mutated fields from the first fit)
  bootstrap_template <- unserialize(serialize(model, NULL))
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
  if (is.null(model$design_weights) || length(model$design_weights) != nrow(model$x)) {
    model$design_weights <- rep(1, nrow(model$x))
  }

  respondent_mask <- !is.na(model$x[, model$col_y])
  model$respondent_mask <- respondent_mask
  scaling_weights <- model$design_weights
  # When scaling we only want respondent rows to contribute; pass a mask so the
  # helper can zero-out nonrespondents without reallocating weights
  weight_mask <- if (length(respondent_mask) == nrow(model$x)) respondent_mask else NULL

  scaling_result <- validate_and_apply_nmar_scaling(
    standardize = model$standardize,
    has_aux = has_aux,
    response_model_matrix_unscaled = model$x[,c(model$col_y,model$cols_delta),drop=FALSE],
    auxiliary_matrix_unscaled = model$x[,model$cols_y_observed,drop=FALSE],
    mu_x_unscaled = filtered_aux_means,
    weights = scaling_weights,
    weight_mask = weight_mask
  )
  model$nmar_scaling_recipe <- scaling_result$nmar_scaling_recipe
  response_model_matrix_scaled <- scaling_result$response_model_matrix_scaled
  auxiliary_matrix_scaled <- scaling_result$auxiliary_matrix_scaled
  mu_x_scaled <- scaling_result$mu_x_scaled

  model$x_1 <- response_model_matrix_scaled[respondent_mask,,drop=FALSE] #observed
  model$x_0 <- response_model_matrix_scaled[!respondent_mask,,drop=FALSE] #unobserved
  model$y_1 <- if (nrow(model$x_1)) model$x_1[, model$col_y, drop = TRUE] else numeric(0) #observed y
  model$x_for_y_obs <- auxiliary_matrix_scaled[respondent_mask,,drop=FALSE]
  model$x_for_y_unobs <- auxiliary_matrix_scaled[!respondent_mask,,drop=FALSE]
  # Track the current scale of feature matrices. We fit f1(.) on the scaled
  # space and compute EM steps there. After unscaling coefficients for
  # presentation we flip this flag to FALSE so downstream density evaluations
  # can re-apply the same scaling recipe to inputs as needed
  model$features_are_scaled <- TRUE


  # model$respondent_weights <- weights(model$x_1)
  # Keep weight vectors even when one of the groups is empty so downstream code
  # (score, variance, bootstrap guards) can safely rely on zero-length vectors
  model$respondent_weights <- if (any(respondent_mask)) model$design_weights[respondent_mask] else numeric(0)
  model$nonrespondent_weights <- if (any(!respondent_mask)) model$design_weights[!respondent_mask] else numeric(0)

  model$theta = stats::runif(length(model$cols_delta) + 2, 0, 0.1)
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
  model$O_matrix_nieobs <- generate_Odds(model)

  #const
  model$f_matrix_nieobs <- generate_conditional_density_matrix(model)
  model$C_matrix_nieobs <- generate_C_matrix(model)



  target_function <- function(theta) {
    model$theta <<- theta
    step_func(model,theta, model$O_matrix_nieobs)
  }


  solution <- nleqslv(
    x = model$theta,
    fn = target_function,
    method = "Newton",
    control = list(
      maxit = 1

    )
  )

  theta_prev=model$theta
  model$theta <- solution$x
  model$loss_value <- solution$fvec
  iter=0

  while(sum((model$theta-theta_prev)^2) > model$tol_value || (iter < model$min_iter && iter < model$max_iter)) {

    solution <- nleqslv(
      x = model$theta,
      fn = target_function,
      method = "Newton",
      control = list(
        maxit = 1

      )

    )

    theta_prev <- model$theta
    model$theta <- solution$x
    model$loss_value <- solution$fvec
    # O_matrix_nieobs <- generate_Odds(theta_prev,x_0[,cols_delta],y_1)
    #
    model$O_matrix_nieobs <- generate_Odds(model)
    iter<-iter+1



  }
  model$iterations <- iter

  if (model$standardize) {
    unscale <- unscale_coefficients(model$theta, matrix(0, length(model$theta), length(model$theta)), model$nmar_scaling_recipe)
    model$theta <- unscale$coefficients
  }

  model$x_1 <- model$x[respondent_mask,,drop=FALSE]
  model$x_0 <- model$x[!respondent_mask,,drop=FALSE]
  model$y_1 <- if (nrow(model$x_1)) model$x_1[, model$col_y, drop = TRUE] else numeric(0)

  model$x_for_y_obs <- model$x_1[, model$cols_y_observed, drop = FALSE]
  model$x_for_y_unobs <- model$x_0[, model$cols_y_observed, drop = FALSE]
  # From this point the feature matrices are on the original (unscaled) space
  # density helpers will internally re-apply the scaling recipe captured during
  # model fitting so that gamma_hat remains consistent
  model$features_are_scaled <- FALSE
  model$respondent_weights <- if (any(respondent_mask)) model$design_weights[respondent_mask] else numeric(0)
  model$nonrespondent_weights <- if (any(!respondent_mask)) model$design_weights[!respondent_mask] else numeric(0)

  model$O_matrix_nieobs <- generate_Odds(model)
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
    # The analytic delta variance we inherited assumes IID sampling: it plugs
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
    n_resp <- length(model$respondent_weights)
    design_varies <- n_resp > 1 && max(abs(model$respondent_weights - model$respondent_weights[1])) > 1e-6
    few_resp <- n_resp < 40
    if (design_varies || few_resp) {
      warning("Delta variance may be unreliable with the current sample; using bootstrap instead.", call. = FALSE)
      use_bootstrap <- TRUE
    }
  }

  if (use_bootstrap) {
    var_results <- list(var_est = NA_real_, vcov = default_vcov)
  }

  if(use_bootstrap){

    bootstrap_runner <- function(data, ...) {
      template <- unserialize(serialize(bootstrap_template, NULL))
      if (is.null(template$standardize)) template$standardize <- TRUE
      exptilt(data, template, ...)
    }

    base_args <- list(
      data = if (model$is_survey) model$design else model$x,
      estimator_func = bootstrap_runner,
      point_estimate = estim_mean(model),
      bootstrap_reps = model$bootstrap_reps
    )
    if (!model$is_survey) {
      # Guard: resample until at least one respondent is present so the EM score
      # and variance code can evaluate the replicate. The bootstrap helper will
      # fall back to NA (with a warning) if this fails repeatedly
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
    std_error = se_final,
    coefficients = model$theta,
    vcov = var_results$vcov,
    model = model,
    converged = TRUE,
    weights = model$respondent_weights,
    variance_message = NA_character_
  )

  return(validate_nmar_result(result, "nmar_result_exptilt"))
}
