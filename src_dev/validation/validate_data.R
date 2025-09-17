#' @title Validate Data for NMAR Analysis
#' @description A robust function to validate a data frame or survey object before performing NMAR analysis.
#' The function checks for common errors like missing variables, data type inconsistencies,
#' and inappropriate variable overlaps. It provides detailed, actionable error messages to facilitate debugging.
#'
#' @param data A data frame or a survey object.
#' @param outcome_variable A string specifying the outcome variable, which is expected to contain NA values.
#' @param covariates_for_outcome A character vector of covariates explaining the outcome.
#' @param covariates_for_missingness A character vector of covariates explaining missingness,
#'   which must be distinct from `covariates_for_outcome`.
#' @return Returns `invisible(NULL)` on success, stopping with a descriptive error on failure.
validate_data <- function(data, outcome_variable, covariates_for_outcome, covariates_for_missingness = c()) {

  # Validate data object type
  if (!inherits(data, c("data.frame", "survey.design"))) {
    stop("'data' must be a data.frame or survey.design object. Received: ", class(data)[1])
  }

  # Extract variables from survey object
  if (inherits(data, "survey.design")) {
    data <- data$variables
  }

  # Check for empty data
  if (nrow(data) == 0) {
    stop("Input dataset is empty (0 rows).")
  }

  # Combine all required variables
  all_vars <- c(outcome_variable, covariates_for_outcome, covariates_for_missingness)

  # Check for duplicate variable specifications
  duplicates <- all_vars[duplicated(all_vars)]
  if (length(duplicates) > 0) {
    stop("Duplicate variables found in input lists: ", paste(duplicates, collapse = ", "))
  }

  # Check variable existence in data
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Validate outcome variable
  if (!is.numeric(data[[outcome_variable]])) {
    bad_val <- data[[outcome_variable]][which(!is.numeric(data[[outcome_variable]]))[1]]
    stop(
      "Outcome variable '", outcome_variable, "' must be numeric.\n",
      "First invalid value: '", bad_val, "' at row ", which(!is.numeric(data[[outcome_variable]]))[1]
    )
  }

  # Check for required NAs in outcome
  if (!anyNA(data[[outcome_variable]])) {
    stop("Outcome variable '", outcome_variable, "' must contain NA values.")
  }

  # Check for at least one non-NA value in outcome
  if (all(is.na(data[[outcome_variable]]))) {
    stop("Outcome variable '", outcome_variable, "' cannot be entirely NA.")
  }

  # Validate covariates
  covariate_vars <- c(covariates_for_outcome, covariates_for_missingness)

  for (var in covariate_vars) {
    # Check type
    if (!is.numeric(data[[var]]) && !is.logical(data[[var]])) {
      bad_val <- data[[var]][which(!is.numeric(data[[var]]) & !is.logical(data[[var]]))[1]]
      stop(
        "Covariate '", var, "' must be numeric or logical.\n",
        "First invalid value: '", bad_val, "' at row ", which(!is.numeric(data[[var]]) & !is.logical(data[[var]]))[1]
      )
    }

    # Check for NAs
    if (anyNA(data[[var]])) {
      stop(
        "Covariate '", var, "' contains NA values.\n",
        "First NA at row ", which(is.na(data[[var]]))[1]
      )
    }
  }

  # Check for overlaps between covariate sets
  if (length(intersect(covariates_for_outcome, covariates_for_missingness)) > 0) {
    common <- intersect(covariates_for_outcome, covariates_for_missingness)
    stop(
      "Covariate sets must be mutually exclusive. Overlapping variables: ",
      paste(common, collapse = ", ")
    )
  }

  invisible(NULL)
}
