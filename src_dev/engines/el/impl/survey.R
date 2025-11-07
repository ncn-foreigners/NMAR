#' Empirical likelihood for survey designs (NMAR)
#' @description Internal method dispatched by `el()` when `data` is a `survey.design`.
#'   Variance via bootstrap is supported. Analytical delta variance for EL is
#'   not implemented and returns NA when requested.
#' @param data A `survey.design` created with [survey::svydesign()].
#' @param formula Two-sided formula: NA-valued outcome on LHS; auxiliaries on RHS.
#' @param auxiliary_means Named numeric vector of population means for auxiliaries.
#' @param standardize Logical; standardize predictors.
#' @param trim_cap Numeric; cap for EL weights (Inf = no trimming).
#' @param control List; solver control for `nleqslv(control=...)`.
#' @param on_failure Character; "return" or "error" on solver failure.
#' @param variance_method Character; "delta" or "bootstrap".
#' @param bootstrap_reps Integer; reps when `variance_method = "bootstrap"`.
#' @param ... Passed to solver.
#' @details Implements the empirical likelihood estimator with design weights.
#'   If \code{n_total} is supplied, design weights are rescaled internally to
#'   ensure \code{sum(weights(design))} and \code{n_total} are on the same scale;
#'   this guarantees the response-multiplier formula uses consistent totals. If
#'   \code{n_total} is not supplied, \code{sum(weights(design))} is used as the
#'   population total \code{N_pop}. When respondents-only designs are used (no
#'   NA in the outcome), \code{n_total} must be provided; if auxiliaries are
#'   requested you must also provide population auxiliary means via
#'   \code{auxiliary_means}. Result weights are the unnormalized EL masses
#'   \code{d_i/D_i(theta)} on this design scale;\code{weights(result, scale = "population")} sums to \code{N_pop}.
#' @return `c('nmar_result_el','nmar_result')`.
#'
#' @keywords internal
el.survey.design <- function(data, formula,
                             auxiliary_means = NULL, standardize = TRUE,
                             trim_cap = Inf, control = list(),
                             on_failure = c("return", "error"),
                             variance_method = c("delta", "bootstrap", "none"),
                             bootstrap_reps = 500,
                             n_total = NULL, start = NULL, trace_level = 0, ...) {
  cl <- match.call()
  on_failure <- match.arg(on_failure)
  if (is.null(variance_method)) variance_method <- "none"
  variance_method <- match.arg(variance_method)
# Coerce unsupported mode to 'none' locally (engine-level warning already issued when called via nmar)
  if (identical(variance_method, "delta")) variance_method <- "none"


  design <- data

# If respondents-only design is supplied, enforce auxiliary_means first (if auxiliaries requested),
# then require n_total for scale coherence.
  outcome_var_check <- all.vars(formula[[2]])
  rhs <- formula[[3]]
  aux_expr <- rhs
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) aux_expr <- rhs[[2L]]
  has_aux_rhs <- length(all.vars(aux_expr)) > 0
  respondents_only_0 <- length(outcome_var_check) == 1 && !anyNA(design$variables[[outcome_var_check]])
  if (respondents_only_0 && has_aux_rhs && is.null(auxiliary_means)) {
    stop(
      paste0(
        "Respondents-only survey design (no NAs in outcome) and auxiliary constraints were requested, ",
        "but 'auxiliary_means' was not provided. Provide population auxiliary means via auxiliary_means=."
      ),
      call. = FALSE
    )
  }
  if (respondents_only_0 && is.null(n_total)) {
    stop("Respondents-only survey design detected (no NAs in outcome), but 'n_total' was not provided. Set el_engine(n_total = <total design weight or population total>).", call. = FALSE)
  }

  parsed_inputs <- prepare_el_inputs(formula, design$variables,
                                     require_na = is.null(n_total))
  design$variables <- parsed_inputs$data
  internal_formula <- parsed_inputs$formula_list
  response_var <- all.vars(internal_formula$response)[1]
  observed_mask <- design$variables[[response_var]] == 1
  observed_indices <- which(observed_mask)
  resp_design <- subset(design, observed_mask)

# (Guard handled earlier.)

# Scale coherence: ensure N_pop and design weights are on the same scale
  design_weight_sum <- sum(weights(design))

  if (!is.null(n_total)) {
    N_pop <- n_total
    scale_factor <- N_pop / design_weight_sum
    scale_mismatch_pct <- abs(scale_factor - 1) * 100

# Graduated warnings based on severity
    if (scale_mismatch_pct > 10) {
# Large mismatch (>10%): likely user error
      warning(sprintf(
        paste0(
          "Large scale mismatch detected (%.1f%%):\n",
          "  User-supplied n_total: %g\n",
          "  Design sum(weights):    %g\n",
          "  Scale factor:           %.4f\n\n",
          "This likely indicates a data preparation error.\n",
          "Verify that n_total and design weights are on the same scale.\n",
          "For example, if weights were rescaled to mean=1 for other analyses,\n",
          "provide n_total on that same rescaled scale, not the original population size."
        ),
        scale_mismatch_pct, n_total, design_weight_sum, scale_factor
      ), call. = FALSE)
      scale_mismatch_detected <- TRUE

    } else if (scale_mismatch_pct > 1) {
# Moderate rescaling (1-10%): upgrade to warning for visibility
      warning(sprintf(
        paste0(
          "Scale mismatch detected (%.1f%%):\n",
          "  User-supplied n_total: %g\n",
          "  Design sum(weights):    %g\n",
          "  Scale factor:           %.4f\n\n",
          "Design weights will be rescaled automatically to ensure internal coherence.\n",
          "This ensures the Lagrange multiplier formula lambda_W = (N_pop/sum(d_i) - 1)/(1 - W)\n",
          "uses consistent scales. Estimates remain unbiased.\n\n",
          "If this is unexpected, verify that n_total and design weights are on the same scale."
        ),
        scale_mismatch_pct, n_total, design_weight_sum, scale_factor
      ), call. = FALSE)
      scale_mismatch_detected <- TRUE

    } else {
# Negligible (<1%) - likely rounding, no message
      scale_mismatch_detected <- FALSE
    }

# Apply scaling to respondent weights
    respondent_weights <- weights(resp_design) * scale_factor

  } else {
# No n_total supplied: use design total as population size
    N_pop <- design_weight_sum
    respondent_weights <- weights(resp_design)
    scale_factor <- 1.0
    scale_mismatch_detected <- FALSE
    scale_mismatch_pct <- 0
  }

  user_args <- list(
    formula = formula,
    auxiliary_means = auxiliary_means, standardize = standardize,
    trim_cap = trim_cap, control = control, ...
  )

  core_results <- el_estimator_core(
    full_data = design,
    respondent_data = resp_design$variables,
    respondent_weights = respondent_weights, N_pop = N_pop,
    internal_formula = internal_formula, auxiliary_means = auxiliary_means,
    standardize = standardize, trim_cap = trim_cap, control = control,
    on_failure = on_failure,
    variance_method = variance_method, bootstrap_reps = bootstrap_reps,
    user_args = user_args, start = start, trace_level = trace_level, ...
  )

  sample_info <- list(
    outcome_var = all.vars(internal_formula$outcome)[1],
    response_var = response_var,
    formula = formula,
    nobs = nrow(design$variables),
    nobs_resp = length(observed_indices),
    n_total = N_pop, # Store N_pop for weights() method
    is_survey = TRUE,
    design = design,
    variance_method = variance_method,
    scale_factor = scale_factor, # Store scale diagnostics
    scale_mismatch_detected = scale_mismatch_detected,
    scale_mismatch_pct = scale_mismatch_pct
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
      sample = list(n_total = N_pop, n_respondents = sample_info$nobs_resp, is_survey = TRUE, design = design),
      inference = list(variance_method = variance_method, df = NA_real_, message = msg),
      diagnostics = diag_list,
      meta = list(engine_name = "empirical_likelihood", call = cl, formula = formula),
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
