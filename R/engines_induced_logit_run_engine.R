#' Run method for induced-logit engine
#'
#' @param engine An object of class `nmar_engine_induced_logit`.
#' @param formula A two-sided formula passed through by `nmar()`.
#' @param data A `data.frame` or `survey.design`.
#' @param trace_level Unused (reserved for future verbosity controls).
#'
#' @return An object of class `nmar_result_induced_logit`.
#' @keywords internal
#' @noRd
#' @exportS3Method run_engine nmar_engine_induced_logit
run_engine.nmar_engine_induced_logit <- function(engine, formula, data, trace_level = 0) {
  args <- list(
    data = data,
    formula = formula,
    variance_method = engine$variance_method,
    bootstrap_reps = engine$bootstrap_reps %||% 500,
    standardize = engine$standardize,
    control = engine$control %||% list(),
    on_failure = engine$on_failure,
    survey_design_policy = engine$survey_design_policy %||% "strict",
    keep_fits = engine$keep_fits %||% FALSE
  )

  res <- do.call(induced_logit, args)

  if (!inherits(res, "nmar_result_induced_logit")) {
    stop("Induced-logit engine did not return an 'nmar_result_induced_logit' object.")
  }
  if (is.list(res$meta)) {
    engine_call <- res$meta$call %||% NULL
    res$meta$engine_call <- engine_call
  }

  res
}
