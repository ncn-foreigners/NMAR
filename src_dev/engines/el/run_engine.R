#' Run method for EL engine
#' @keywords internal
#' @exportS3Method run_engine nmar_engine_el
run_engine.nmar_engine_el <- function(engine, formula, data, trace_level = 0) {
# Build argument list for EL implementation directly from engine config
  args <- list(
    data = data,
    formula = formula,
    auxiliary_means = engine$auxiliary_means,
    standardize = engine$standardize,
    strata_augmentation = engine$strata_augmentation %||% TRUE,
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

# Dispatch to EL implementation (data.frame or survey.design)
  res <- do.call(el, args)

# Ensure class includes the NMAR parent for downstream compatibility
  if (!inherits(res, "nmar_result_el")) {
    stop("EL engine did not return an 'nmar_result_el' object.")
  }
# Patch result metadata: expose outer nmar() call as meta$call and engine call as meta$engine_call
  if (is.list(res$meta)) {
    engine_call <- res$meta$call %||% NULL
    res$meta$engine_call <- engine_call
  }
  res
}
