#' Empirical likelihood for data frames (NMAR)
#' @description Internal method dispatched by `el()` when `data` is a `data.frame`.
#'   Returns `c('nmar_result_el','nmar_result')` with the point estimate, optional
#'   bootstrap SE, weights, coefficients, diagnostics, and metadata.
#' @param data A `data.frame` where the outcome column contains `NA` for nonrespondents.
#' @param formula Two-sided formula `Y_miss ~ auxiliaries`.
#' @param auxiliary_means Named numeric vector of population means for auxiliary
#'   design columns. Names must match the materialized `model.matrix` columns on
#'   the first RHS (after formula expansion), including factor indicators and
#'   transformed terms. The intercept is always excluded.
#' @param standardize Logical; whether to standardize predictors prior to estimation.
#' @param trim_cap Numeric; cap for EL weights (`Inf` = no trimming).
#' @param control List; optional solver control parameters for `nleqslv(control=...)`.
#' @param on_failure Character; one of `"return"` or `"error"` on solver failure.
#' @param variance_method Character; one of `"delta"`, `"bootstrap"`, or `"none"`.
#' @param bootstrap_reps Integer; number of bootstrap reps if `variance_method = "bootstrap"`.
#' @param n_total Optional integer population size. When the outcome contains
#'   at least one `NA`, `n_total` defaults to `nrow(data)`; when respondents-only
#'   data are supplied (no `NA` in the outcome), `n_total` must be provided.
#' @param start Optional list of starting values passed to the solver helpers.
#' @param trace_level Integer 0-3 controlling estimator logging detail.
#' @param family Missingness (response) model family specification (defaults to the logit bundle).
#' @param ... Additional arguments passed to the solver.
#' @details Implements the empirical likelihood estimator for IID data with
#'   optional auxiliary moment constraints. The missingness-model score is the
#'   Bernoulli derivative with respect to the linear predictor, supporting logit
#'   and probit links. When respondents-only data are supplied (no `NA` in the
#'   outcome), `n_total` is required so the response-rate equation targets the
#'   full sample size. When missingness is observed (`NA` present), the default
#'   population total is `nrow(data)`. If respondents-only data are used and
#'   auxiliaries are requested, you must also provide population auxiliary
#'   means via `auxiliary_means`. Result weights are the unnormalized EL
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
  if (identical(variance_method, "delta")) variance_method <- "none"


  design <- el_prepare_design(
    formula = formula,
    data = data,
    require_na = FALSE
  )

  extra_args <- list(...)

  el_run_core_analysis(
    call = cl,
    formula = formula,
    raw_data = data,
    design = design,
    weights_full = NULL,
    n_total = n_total,
    variance_method = variance_method,
    is_survey = FALSE,
    design_object = NULL,
    auxiliary_means = auxiliary_means,
    standardize = standardize,
    trim_cap = trim_cap,
    control = control,
    on_failure = on_failure,
    family = family,
    bootstrap_reps = bootstrap_reps,
    start = start,
    trace_level = trace_level,
    extra_user_args = extra_args
  )
}
