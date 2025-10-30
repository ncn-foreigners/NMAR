#' @importFrom stats weights
#' @exportS3Method exptilt survey.design
exptilt.survey.design <- function(data, formula,
                                  auxiliary_means = NULL,
                                  standardize = TRUE,
                                  prob_model_type = c("logit", "probit"),
                                  y_dens = c("auto", "normal", "lognormal", "exponential"),
                                  variance_method = c("delta", "bootstrap"),
                                  bootstrap_reps = 10,
                                  control = list(),
                                  stopping_threshold = 1,
                                  on_failure = c("return", "error"),
                                  supress_warnings = FALSE,
                                  trace_level = 0,
                                  ...) {
  design_vars <- data$variables
  design_weights <- as.numeric(stats::weights(data))
# The survey method is a thin adapter that harvests the analysis weights and
# then reuses the data.frame implementation so that both flows share identical
# preprocessing, scaling, EM fitting and bootstrap variance logic

# Delegate to the data.frame method after harvesting weights so both flows
# share identical preprocessing, scaling, and bootstrap behavior
  exptilt.data.frame(
    design_vars,
    formula = formula,
    auxiliary_means = auxiliary_means,
    standardize = standardize,
    prob_model_type = prob_model_type,
    y_dens = y_dens,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    control = control,
    stopping_threshold = stopping_threshold,
    on_failure = on_failure,
    supress_warnings = supress_warnings,
    design_weights = design_weights,
    survey_design = data,
    trace_level = trace_level,
    ...
  )
}
