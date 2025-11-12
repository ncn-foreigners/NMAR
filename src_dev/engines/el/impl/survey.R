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
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193-200.
#'
#' Wu, C., and Sitter, R. R. (2001). A model-calibration approach to using complete
#' auxiliary information from survey data. Journal of the American Statistical Association,
#' 96(453), 185-193.
#' @return `c('nmar_result_el','nmar_result')`.
#'
#' @name el_survey
#' @keywords internal
el.survey.design <- function(data, formula,
                             auxiliary_means = NULL, standardize = TRUE,
                             trim_cap = Inf, control = list(),
                             on_failure = c("return", "error"),
                             variance_method = c("delta", "bootstrap", "none"),
                             bootstrap_reps = 500,
                             n_total = NULL, start = NULL, trace_level = 0,
                             family = logit_family(), ...) {
  cl <- match.call()
  on_failure <- match.arg(on_failure)
  if (is.null(variance_method)) variance_method <- "none"
  variance_method <- match.arg(variance_method)
# Coerce unsupported mode to 'none' locally (engine-level warning already issued when called via nmar)
  if (identical(variance_method, "delta")) variance_method <- "none"


  design <- data

  respondents_only <- el_validate_respondents_only(formula, design$variables, auxiliary_means, context_label = "survey design")
  if (respondents_only && is.null(n_total)) {
    stop("Respondents-only survey design detected (no NAs in outcome), but 'n_total' was not provided. Set el_engine(n_total = <total design weight or population total>).", call. = FALSE)
  }

# Prepare inputs, appending survey design context to any validation error
  el_get_design_context <- function(design) {
    ctx <- list(ids = "<unspecified>", strata = "<unspecified>")
    dc <- try(getCall(design), silent = TRUE)
    if (!inherits(dc, "try-error") && !is.null(dc)) {
      args <- as.list(dc)[-1]
      get_arg <- function(nm) if (!is.null(args[[nm]])) args[[nm]] else NULL
      ids_val <- get_arg("ids"); id_val <- get_arg("id")
      ids_obj <- if (!is.null(ids_val)) ids_val else id_val
      strata_obj <- get_arg("strata")
      deparse1 <- function(x) paste(deparse(x, width.cutoff = 200L), collapse = " ")
      if (!is.null(ids_obj)) ctx$ids <- deparse1(ids_obj)
      if (!is.null(strata_obj)) ctx$strata <- deparse1(strata_obj)
    }
    ctx
  }

  survey_ctx <- el_get_design_context(design)

  spec <- tryCatch(
    el_prepare_inputs(
      formula = formula,
      data = design$variables,
      require_na = is.null(n_total),
      auxiliary_means = auxiliary_means
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      msg2 <- paste0(msg, sprintf("\nSurvey design info: ids = %s, strata = %s", survey_ctx$ids, survey_ctx$strata))
      stop(msg2, call. = FALSE)
    }
  )
  design$variables <- spec$data

# Scale coherence: ensure N_pop and design weights are on the same scale
  weights_all <- as.numeric(weights(design))
  design_weight_sum <- sum(weights_all)

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

    respondent_weights_full <- weights_all * scale_factor

  } else {
# No n_total supplied: use design total as population size
    N_pop <- design_weight_sum
    respondent_weights_full <- weights_all
    scale_factor <- 1.0
    scale_mismatch_detected <- FALSE
    scale_mismatch_pct <- 0
  }

  context <- el_build_context(
    prepared_inputs = spec,
    full_data = design,
    formula = formula,
    respondent_weights_full = respondent_weights_full,
    N_pop = N_pop,
    is_survey = TRUE,
    design = design,
    variance_method = variance_method
  )

  aux_summary <- el_resolve_auxiliaries(
    spec$auxiliary_matrix,
    spec$auxiliary_matrix_full,
    auxiliary_means,
    weights_full = weights_all
  )

  core_results <- el_estimator_core(
    full_data = context$full_data,
    respondent_data = context$respondent_data,
    respondent_weights = context$respondent_weights,
    N_pop = context$N_pop,
    response_matrix = spec$response_matrix,
    auxiliary_matrix = aux_summary$matrix,
    mu_x = aux_summary$means,
    standardize = standardize,
    trim_cap = trim_cap,
    control = control,
    on_failure = on_failure,
    family = family,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    user_args = list(
      formula = formula,
      auxiliary_means = auxiliary_means,
      standardize = standardize,
      trim_cap = trim_cap,
      control = control,
      ...
    ),
    start = start,
    trace_level = trace_level,
    outcome_var = spec$outcome_var,
    has_aux = aux_summary$has_aux,
    auxiliary_means = auxiliary_means
  )

  el_build_result(core_results, context$analysis_info, cl, formula)
}
