#' @importFrom stats weights
#' @export
exptilt.survey.design <- function(data, model, on_failure = c("return"), ...) {
  on_failure <- match.arg(on_failure)

  # Keep the survey design on the model so diagnostics/variance can reuse it
  model$design <- data
  model$is_survey <- TRUE

  # Trim the design variables to the columns selected upstream (run_engine)
  design_vars <- data$variables

  # Extract the analysis weights once so every downstream component sees the
  # same vector. Use the survey package's helper to honor calibration etc.
  # stats::weights() dispatches to survey::weights for survey designs, which
  # honors calibration / analysis weight choices
  survey_weights <- stats::weights(data)
  model$design_weights <- as.numeric(survey_weights)

  # Reuse the data.frame implementation now that we have the design weights
  # and the same column structure, so survey and IID flows share one code path
  exptilt.data.frame(design_vars, model, on_failure = on_failure)
}
