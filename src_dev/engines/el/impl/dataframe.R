#' Empirical likelihood for data frames (NMAR)
#' @description Internal method dispatched by `el()` when `data` is a `data.frame`.
#'   Returns `c('nmar_result_el','nmar_result')` with the point estimate, optional
#'   bootstrap SE, weights, coefficients, diagnostics, and metadata.
#' @param data A `data.frame` where the outcome column contains `NA` for nonrespondents.
#' @param formula Two-sided formula `Y_miss ~ auxiliaries`.
#' @param auxiliary_means Named numeric vector of population means for auxiliary variables (names must match RHS of outcome formula).
#' @param standardize Logical; whether to standardize predictors prior to estimation.
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
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193-200.
#'
#' Wu, C., and Sitter, R. R. (2001). A model-calibration approach to using complete
#' auxiliary information from survey data. Journal of the American Statistical Association,
#' 96(453), 185-193.
#'
#' @name el_dataframe
#' @keywords internal
el.data.frame <- function(data, formula,
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
# Coerce unsupported mode to 'none' locally (engine-level warning already issued)
  if (identical(variance_method, "delta")) variance_method <- "none"


  respondents_only <- el_validate_respondents_only(formula, data, auxiliary_means, context_label = "data frame")
# el_prepare_inputs(require_na = is.null(n_total)) still enforces the NA requirement.

  spec <- el_prepare_inputs(
    formula = formula,
    data = data,
    require_na = is.null(n_total),
    auxiliary_means = auxiliary_means
  )

  context <- el_build_context(
    prepared_inputs = spec,
    full_data = spec$data,
    formula = formula,
    respondent_weights_full = NULL,
    N_pop = n_total,
    is_survey = FALSE,
    design = NULL,
    variance_method = variance_method
  )

  aux_summary <- el_resolve_auxiliaries(
    spec$auxiliary_matrix,
    spec$auxiliary_matrix_full,
    auxiliary_means,
    weights_full = NULL
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
