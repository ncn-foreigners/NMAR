#' @exportS3Method run_engine nmar_engine_exptilt_nonparam
run_engine.nmar_engine_exptilt_nonparam <- function(engine, formula, data, trace_level) {


args <- list(
  data = data,
  formula = formula,
  trace_level = trace_level,
  refusal_col = engine$refusal_col,
  max_iter = engine$max_iter,
  tol_value = engine$tol_value,
  design_weights = NULL
)
# print(exptilt_nonparam)
# browser()
res <- do.call(exptilt_nonparam, args)
if (!inherits(res, "nmar_result_exptilt_nonparam")) {
  stop("Exptilt engine did not return an 'nmar_result_exptilt_nonparam' object.")
}
res
}
