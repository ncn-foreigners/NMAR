#' @importFrom stats weights
#' @exportS3Method exptilt_nonparam survey.design
exptilt_nonparam.survey.design <- function(
    data,
    formula,
    trace_level = 0,
    refusal_col,
    max_iter = 100,
    tol = 1e-6,
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
    tol = tol,
    design_weights = design_weights,
    ...
  )
}
