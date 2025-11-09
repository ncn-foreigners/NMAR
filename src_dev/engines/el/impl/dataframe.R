#' Empirical likelihood for data frames (NMAR)
#' @description Internal method dispatched by `el()` when `data` is a `data.frame`.
#'   Returns `c('nmar_result_el','nmar_result')` with the point estimate, optional
#'   bootstrap SE, weights, coefficients, diagnostics, and metadata.
#' @param data A `data.frame` where the outcome column contains `NA` for nonrespondents.
#' @param formula Two-sided formula `Y_miss ~ auxiliaries`.
#' @param auxiliary_means Named numeric vector of population means for auxiliary variables (names must match RHS of outcome formula).
#' @param standardize Logical; whether to standardize predictors prior to estimation.
#' @param design_matrices Internal list of precomputed design matrices produced by
#'   the NMAR input pipeline. Users should leave this as `NULL`; engines call the
#'   method with populated matrices to avoid recomputing model matrices.
#' @param trim_cap Numeric; cap for EL weights (`Inf` = no trimming).
#' @param control List; optional solver control parameters for `nleqslv(control=...)`.
#' @param on_failure Character; one of `"return"` or `"error"` on solver failure.
#' @param variance_method Character; one of `"delta"`, `"bootstrap"`, or `"none"`.
#' @param bootstrap_reps Integer; number of bootstrap reps if `variance_method = "bootstrap"`.
#' @param ... Additional arguments passed to the solver.
#' @details Implements the empirical likelihood estimator for IID data with
#'   optional auxiliary moment constraints. The response-model score is the
#'   Bernoulli derivative with respect to the linear predictor, supporting logit
#'   and probit links. When respondents-only data is supplied (no NA in the
#'   outcome), set \code{n_total} to the total number of sampled units; otherwise
#'   the total is taken as \code{nrow(data)}. If respondents-only data is used
#'   and auxiliaries are requested, you must also provide population auxiliary
#'   means via \code{auxiliary_means}. Result weights are the unnormalized EL
#'   masses \code{d_i/D_i(theta)} on this analysis scale.
#'
#' @keywords internal
el.data.frame <- function(data, formula,
                          user_formula = formula,
                          auxiliary_means = NULL, standardize = TRUE,
                          design_matrices,
                          outcome_label,
                          trim_cap = Inf, control = list(),
                          on_failure = c("return", "error"),
                          variance_method = c("delta", "bootstrap", "none"),
                          bootstrap_reps = 500,
                          n_total = NULL, start = NULL, trace_level = 0,
                          design_info = NULL, ...) {
  cl <- match.call()
  on_failure <- match.arg(on_failure)
  if (is.null(variance_method)) variance_method <- "none"
  variance_method <- match.arg(variance_method)
# Coerce unsupported mode to 'none' locally (engine-level warning already issued)
  if (identical(variance_method, "delta")) variance_method <- "none"


  if (is.null(design_matrices)) {
    stop("`design_matrices` must be supplied by run_engine().", call. = FALSE)
  }
  if (is.null(design_info)) {
    stop("`design_info` must be supplied by run_engine().", call. = FALSE)
  }

  runtime_inputs <- el_build_runtime_inputs(
    data = data,
    design_info = design_info,
    auxiliary_means = auxiliary_means,
    n_total = n_total,
    require_na = is.null(n_total),
    context = "data.frame"
  )
  estimation_data <- runtime_inputs$data
  internal_formula <- runtime_inputs$internal_formula
  response_var <- runtime_inputs$response_var
  observed_indices <- runtime_inputs$observed_indices
  outcome_name <- runtime_inputs$outcome_name
  precomputed_design <- el_build_precomputed_design(
    design_matrices = design_matrices,
    estimation_data = estimation_data,
    outcome_var = outcome_name,
    respondent_indices = observed_indices
  )

# (Guard handled above, before any n_total enforcement.)

  respondent_weights <- rep(1, length(observed_indices))
# For data frames: if n_total not specified, use total sample size
# This represents the "population" we're inferring about
  N_pop <- n_total %||% nrow(estimation_data)

  user_args <- list(
    formula = formula,
    auxiliary_means = auxiliary_means, standardize = standardize,
    trim_cap = trim_cap, control = control, ...,
    n_total = n_total
  )

  core_results <- el_estimator_core(
    full_data = estimation_data,
    respondent_data = estimation_data[observed_indices, ],
    respondent_weights = respondent_weights, N_pop = N_pop,
    internal_formula = internal_formula, auxiliary_means = auxiliary_means,
    standardize = standardize, trim_cap = trim_cap, control = control,
    on_failure = on_failure,
    variance_method = variance_method, bootstrap_reps = bootstrap_reps,
    user_args = user_args, start = start, trace_level = trace_level,
    precomputed_design = precomputed_design, ...
  )

  sample_info <- list(
    outcome_var = outcome_name,
    outcome_label = outcome_label,
    response_var = response_var,
    formula = user_formula,
    nobs = nrow(estimation_data),
    nobs_resp = length(observed_indices),
    n_total = N_pop, # Store N_pop for weights() method
    is_survey = FALSE,
    design = NULL,
    variance_method = variance_method
  )

  if (!core_results$converged) {
    diag_list <- core_results$diagnostics
    if (is.null(diag_list)) diag_list <- list()
    msg <- diag_list$message
    if (is.null(msg)) msg <- NA_character_
    result <- new_nmar_result(
      estimate = NA_real_,
      estimate_name = sample_info$outcome_var,
      se = NA_real_,
      converged = FALSE,
      model = list(coefficients = NULL, vcov = NULL),
      weights_info = list(values = numeric(0), trimmed_fraction = NA_real_),
      sample = list(n_total = N_pop, n_respondents = sample_info$nobs_resp, is_survey = FALSE, design = NULL),
      inference = list(variance_method = variance_method, df = NA_real_, message = msg),
      diagnostics = diag_list,
      meta = list(engine_name = "empirical_likelihood", call = cl, formula = user_formula),
      extra = list(nmar_scaling_recipe = core_results$nmar_scaling_recipe),
      class = "nmar_result_el"
    )
    return(validate_nmar_result(result, "nmar_result_el"))
  }

  result <- new_nmar_result_el(
    y_hat = core_results$y_hat, se = core_results$se, weights = core_results$weights,
    coefficients = core_results$coefficients, vcov = core_results$vcov,
    converged = core_results$converged, diagnostics = core_results$diagnostics,
    data_info = sample_info,
    nmar_scaling_recipe = core_results$nmar_scaling_recipe,
    fitted_values = core_results$fitted_values, call = cl
  )

  validate_nmar_result(result, "nmar_result_el")
}
