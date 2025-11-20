#' Extract a strata factor from a survey.design object
#'
#' Uses the original svydesign() call stored in the object to recreate the
#' stratum labels as a single factor. When multiple stratification variables
#' are supplied, their interaction is used.
#' @keywords internal
el_extract_strata_factor <- function(design) {
  if (!inherits(design, "survey.design")) return(NULL)
  dc <- try(getCall(design), silent = TRUE)
  if (inherits(dc, "try-error") || is.null(dc)) return(NULL)
  args <- as.list(dc)[-1L]
  get_arg <- function(nm) if (!is.null(args[[nm]])) args[[nm]] else NULL
  strata_expr <- get_arg("strata")
  if (is.null(strata_expr)) strata_expr <- get_arg("strata")
  if (is.null(strata_expr)) return(NULL)
  vars <- design$variables
# Coerce to formula if needed
  strata_formula <- if (inherits(strata_expr, "formula")) {
    strata_expr
  } else {
    tryCatch(stats::as.formula(strata_expr),
      error = function(e) return(NULL)
    )
  }
  if (is.null(strata_formula)) return(NULL)
  mf <- tryCatch(
    stats::model.frame(strata_formula, data = vars, na.action = stats::na.pass),
    error = function(e) NULL
  )
  if (is.null(mf) || nrow(mf) != nrow(vars)) return(NULL)
  if (ncol(mf) == 1L) {
    as.factor(mf[[1L]])
  } else {
    interaction(mf, drop = TRUE)
  }
}

#' Empirical likelihood for survey designs (NMAR)
#' @description Internal method dispatched by `el()` when `data` is a `survey.design`.
#'   Variance via bootstrap is supported. Analytical delta variance for EL is
#'   not implemented and returns NA when requested.
#' @param data A `survey.design` created with [survey::svydesign()].
#' @param formula Two-sided formula: NA-valued outcome on LHS; auxiliaries on RHS.
#' @param auxiliary_means Named numeric vector of population means for auxiliary
#'   design columns. Names must match the materialized `model.matrix` columns on
#'   the first RHS (after formula expansion), including factor indicators and
#'   transformed terms. The intercept is always excluded.
#' @param standardize Logical; standardize predictors.
#' @param strata_augmentation Logical; when \code{TRUE} (default), augment the
#'   auxiliary design with stratum indicators and stratum shares when a strata
#'   structure is present in the survey design.
#' @param trim_cap Numeric; cap for EL weights (Inf = no trimming).
#' @param control List; solver control for `nleqslv(control=...)`.
#' @param on_failure Character; "return" or "error" on solver failure.
#' @param variance_method Character; "delta", "bootstrap", or "none".
#' @param bootstrap_reps Integer; reps when `variance_method = "bootstrap"`.
#' @param n_total Optional analysis-scale population size \code{N_pop}; required for respondents-only designs.
#' @param start Optional list of starting values passed to solver helpers.
#' @param trace_level Integer 0-3 controlling estimator logging detail.
#' @param family Missingness (response) model family specification (defaults to logit).
#' @param ... Passed to solver.
#' @details Implements the empirical likelihood estimator with design weights.
#'   If \code{n_total} is supplied, it is treated as the analysis-scale population
#'   size \code{N_pop} used in the design-weighted QLS system. If \code{n_total}
#'   is not supplied, \code{sum(weights(design))} is used as \code{N_pop}. Design
#'   weights are not rescaled internally; the EL equations use respondent weights
#'   and \code{N_pop} via \code{T0 = N_pop - sum(d_i)} in the linkage equation.
#'   When respondents-only designs are used (no NA in the outcome), \code{n_total}
#'   must be provided; if auxiliaries are requested you must also provide
#'   population auxiliary means via \code{auxiliary_means}. Result weights are the
#'   unnormalized EL masses \code{d_i/D_i(theta)} on this analysis scale;
#'   \code{weights(result, scale = "population")} sums to \code{N_pop}.
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193-200.
#'
#' @return `c('nmar_result_el','nmar_result')`.
#'
#' @name el_survey
#' @keywords internal
el.survey.design <- function(data, formula,
                             auxiliary_means = NULL, standardize = TRUE,
                             strata_augmentation = TRUE,
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
# Coerce unsupported variance modes to "none"; nmar() already warned at dispatch time.
  if (identical(variance_method, "delta")) variance_method <- "none"


  design <- data
# Capture a readable summary so validation errors can report the original survey call.
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

  weights_initial <- as.numeric(weights(design))
  design_weight_sum <- sum(weights_initial)

  strata_factor <- el_extract_strata_factor(design)
  if (!is.null(strata_factor)) {
    design$variables[["..nmar_strata_factor.."]] <- strata_factor
  }

  if (!is.null(n_total)) {
    N_pop <- n_total
    respondent_weights_full <- weights_initial
  } else {
    N_pop <- design_weight_sum
    respondent_weights_full <- weights_initial
  }

  extra_args <- list(...)

  input_spec <- tryCatch(
    el_build_input_spec(
      formula = formula,
      data = design$variables,
      weights_full = respondent_weights_full,
      population_total = N_pop,
      population_total_supplied = !is.null(n_total),
      is_survey = TRUE,
      design_object = design,
      auxiliary_means = auxiliary_means
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      msg2 <- paste0(msg, sprintf("\nSurvey design info: ids = %s, strata = %s", survey_ctx$ids, survey_ctx$strata))
      stop(msg2, call. = FALSE)
    }
  )

  el_run_core_analysis(
    call = cl,
    formula = formula,
    input_spec = input_spec,
    variance_method = variance_method,
    auxiliary_means = auxiliary_means,
    standardize = standardize,
    strata_augmentation = strata_augmentation,
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
