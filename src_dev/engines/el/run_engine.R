#' Run method for EL engine
#' @keywords internal
#' @exportS3Method run_engine nmar_engine_el
run_engine.nmar_engine_el <- function(engine, task) {
# Reuse the shared design preparation so EL mirrors the ET workflow for
# survey designs, scaling, and auxiliary moment injection
  design_info <- prepare_nmar_design(
    task,
    standardize = engine$standardize,
    auxiliary_means = engine$auxiliary_means,
    include_response = TRUE,
    include_auxiliary = TRUE
  )

# Use the original user formula (possibly partitioned with `|`) so that
# factors and transformations are honored by model.matrix().
  f_use <- task$formula_original

  args <- list(
    data = design_info$survey_design %||% design_info$data,
    formula = f_use,
    auxiliary_means = design_info$auxiliary_means,
    standardize = design_info$standardize,
    n_total = engine$n_total,
    start = engine$start,
    trim_cap = engine$trim_cap,
    control = engine$control,
    on_failure = engine$on_failure,
    variance_method = engine$variance_method,
    bootstrap_reps = engine$bootstrap_reps,
    family = engine$family,
    trace_level = task$trace_level
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
# Use outer nmar() call, but ensure the formula component is the actual formula object
    outer_call <- task$nmar_call %||% res$meta$call
    if (!is.null(outer_call)) {
      outer_call$formula <- task$formula_original %||% outer_call$formula
    }
    res$meta$call <- outer_call
  }
  res
}
