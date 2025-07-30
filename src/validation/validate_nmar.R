validate_nmar <- function(data,outcome_variable,covariates_for_outcome,covariates_for_missingness=c()){
  browser()
  stopifnot(
    "At least one observed variable must be specified in 'covariates_for_outcome'" =
      length(covariates_for_outcome) > 0

,"'outcome_variable' must be a single variable" =
  length(outcome_variable) == 1,

"Column 'outcome_variable' must be present in the data frame" =
  outcome_variable %in% colnames(data),

"Columns in 'covariates_for_outcome' must be present in the data frame" =
  all(covariates_for_outcome %in% colnames(data)),

"Column 'outcome_variable' must not be in 'covariates_for_outcome'" =
  !(outcome_variable %in% covariates_for_outcome),

"Columns in covariates_for_missingness must be present in the data frame" =
  all(covariates_for_missingness %in% colnames(data)),

"Columns in covariates_for_missingness must not be in covariates_for_outcome" =
  all(!(covariates_for_missingness %in% covariates_for_outcome)),

"Outcome variable must not be in covariates_for_missingness" =
  !(outcome_variable %in% covariates_for_missingness),

"Outcome variable must not be in covariaties_for_outcome" =
  !(outcome_variable %in% covariates_for_outcome),

"Data must be a data frame" =
  is.data.frame(data),

"Data must not be empty" =
  nrow(data) > 0,

"All data must be numeric" =
  all(sapply(data, is.numeric) | sapply(data, is.integer) | sapply(data, is.double) | sapply(data, is.logical)),

"NA's should be present in outcome_variable" =
  any(is.na(data[[outcome_variable]]))

  )
}
