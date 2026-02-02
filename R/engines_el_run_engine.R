#' Run method for EL engine
#'
#' @param engine An object of class \code{nmar_engine_el}.
#' @param formula A two-sided formula passed through by \code{nmar()}.
#' @param data A \code{data.frame} or \code{survey.design}.
#' @param trace_level Integer 0-3 controlling verbosity.
#'
#' @return An object of class \code{nmar_result_el}.
#' @keywords internal
#' @exportS3Method run_engine nmar_engine_el
run_engine.nmar_engine_el <- function(engine, formula, data, trace_level = 0) {
  args <- list(
    data = data,
    formula = formula,
    auxiliary_means = engine$auxiliary_means,
    standardize = engine$standardize,
    n_total = engine$n_total,
    start = engine$start,
    trim_cap = engine$trim_cap,
    control = engine$control,
    on_failure = engine$on_failure,
    variance_method = engine$variance_method,
    bootstrap_reps = engine$bootstrap_reps,
    family = engine$family,
    trace_level = trace_level
  )

  if (inherits(data, "survey.design")) {
    args$strata_augmentation <- engine$strata_augmentation %||% TRUE
  }

  res <- do.call(el, args)

  if (!inherits(res, "nmar_result_el")) {
    stop("EL engine did not return an 'nmar_result_el' object.")
  }
  if (is.list(res$meta)) {
    engine_call <- res$meta$call %||% NULL
    res$meta$engine_call <- engine_call
  }

  res
}
