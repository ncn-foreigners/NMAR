#' Rescale survey design weights in place
#' @keywords internal
el_rescale_survey_design_weights <- function(design, scale_factor) {
  if (!inherits(design, "survey.design")) {
    stop("Internal error: expected a survey.design object for weight rescaling.", call. = FALSE)
  }
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    stop("scale_factor must be a positive, finite numeric value.", call. = FALSE)
  }
  design[["prob"]] <- design[["prob"]] / scale_factor
  if (!is.null(design[["allprob"]])) {
    design[["allprob"]] <- design[["allprob"]] / scale_factor
  }
  design
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
#' @param trim_cap Numeric; cap for EL weights (Inf = no trimming).
#' @param control List; solver control for `nleqslv(control=...)`.
#' @param on_failure Character; "return" or "error" on solver failure.
#' @param variance_method Character; "delta" or "bootstrap".
#' @param bootstrap_reps Integer; reps when `variance_method = "bootstrap"`.
#' @param n_total Optional population size used to rescale design weights; required for respondents-only designs.
#' @param start Optional list of starting values passed to solver helpers.
#' @param trace_level Integer 0-3 controlling estimator logging detail.
#' @param family Missingness (response) model family specification (defaults to logit).
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

# Scale coherence: ensure `N_pop` and design weights remain on the same scale.
  weights_initial <- as.numeric(weights(design))
  design_weight_sum <- sum(weights_initial)

  if (!is.null(n_total)) {
    N_pop <- n_total
    scale_factor <- N_pop / design_weight_sum
    scale_mismatch_pct <- abs(scale_factor - 1) * 100

# Emit warnings with severity that reflects the amount of rescaling.
    if (scale_mismatch_pct > 10) {
# Large mismatch (>10%) is most often a data-preparation error.
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
    } else if (scale_mismatch_pct > 1) {
# Moderate rescaling (1-10%) still deserves a warning so users can confirm intent.
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
    } else {
# Negligible (<1%) differences usually stem from rounding, so we stay silent.
    }

    if (!isTRUE(all.equal(scale_factor, 1))) {
      design <- el_rescale_survey_design_weights(design, scale_factor)
    }
    respondent_weights_full <- as.numeric(weights(design))

  } else {
# If no population size is supplied, default to the design-weight total.
    N_pop <- design_weight_sum
    respondent_weights_full <- weights_initial
    scale_factor <- 1.0
    scale_mismatch_pct <- 0
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
