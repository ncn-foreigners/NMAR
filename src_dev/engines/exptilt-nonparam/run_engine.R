#' @exportS3Method NULL
run_engine.nmar_engine_exptilt_nonparam <- function(engine, formula, data, trace_level = 0) {


args <- list(
  data = data,
  formula = formula,
  trace_level = trace_level,
  refusal_col = engine$refusal_col,
  max_iter = engine$max_iter,
  tol = engine$tol_value
)
res <- do.call(exptilt_nonparam, args)
if (!inherits(res, "nmar_result_exptilt_nonparam")) {
  stop("Exptilt engine did not return an 'nmar_result_exptilt_nonparam' object.")
}
res
}
