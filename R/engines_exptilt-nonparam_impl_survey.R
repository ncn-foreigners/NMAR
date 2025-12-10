#' @importFrom stats aggregate
#' @exportS3Method exptilt_nonparam data.frame
exptilt_nonparam.survey <- function(
    data,
    formula,
    trace_level = 0,
    refusal_col,
    max_iter = 100,
    tol_value = 1e-6,
    design_weights = NULL, # Added to match S3 pattern
    ...
) {
  design_vars <- data$variables

  design_weights <- as.numeric(stats::weights(data))

  exptilt_nonparam.data.frame(
    data = design_vars,
    formula = formula,
    trace_level = trace_level,
    refusal_col = refusal_col,
    max_iter = max_iter,
    tol_value = tol_value,
    design_weights = design_weights
  )
}
