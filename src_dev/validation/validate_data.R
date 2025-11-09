#' @title Validate Data for NMAR Analysis
#' @description A robust function to validate a data frame or survey object before performing NMAR analysis.
#' The function checks for common errors like missing variables, data type inconsistencies,
#' and inappropriate variable overlaps. It provides detailed, actionable error messages to facilitate debugging.
#'
#' @param data A data frame or a survey object.
#' @param outcome_variables Character vector listing outcome variables, expected to contain NA values unless respondents-only is allowed.
#' @param covariates_for_outcome A character vector of covariates explaining the outcome.
#' @param covariates_for_missingness A character vector of covariates explaining missingness.
#' @param allow_outcome_in_missingness Logical; allow the outcome to also appear in the
#'   response-model covariates (default `FALSE`).
#' @param allow_covariate_overlap Logical; allow overlap between outcome and response
#'   covariate sets (default `FALSE`).
#' @param allow_respondents_only Logical; allow datasets with no missing outcome
#'   values (respondents-only). When TRUE, the caller is expected to provide any
#'   additional information required by the engine (e.g., total sample size).
#' @return Returns `invisible(NULL)` on success, stopping with a descriptive error on failure.
#' @keywords internal
validate_data <- function(data,
                          outcome_variables,
                          covariates_for_outcome,
                          covariates_for_missingness = character(),
                          allow_outcome_in_missingness = FALSE,
                          allow_covariate_overlap = FALSE,
                          allow_respondents_only = FALSE) {
# Validate data object type
  if (!inherits(data, c("data.frame", "survey.design"))) {
    stop("'data' must be a data.frame or survey.design object. Received: ", class(data)[1], call. = FALSE)
  }

# Extract variables from survey object
  if (inherits(data, "survey.design")) {
    data <- data$variables
  }

# Check for empty data
  if (nrow(data) == 0) {
    stop("Input dataset is empty (0 rows).", call. = FALSE)
  }

  outcome_variables <- unique(outcome_variables)
  if (!is.character(outcome_variables) || !length(outcome_variables)) {
    stop("`outcome_variable` must be a character vector with at least one entry.", call. = FALSE)
  }
  if (!is.character(covariates_for_outcome)) {
    stop("`covariates_for_outcome` must be a character vector.", call. = FALSE)
  }
  if (!is.character(covariates_for_missingness)) {
    stop("`covariates_for_missingness` must be a character vector.", call. = FALSE)
  }

  if (anyDuplicated(covariates_for_outcome)) {
    dup <- unique(covariates_for_outcome[duplicated(covariates_for_outcome)])
    stop("Duplicate variables found in covariates_for_outcome: ", paste(dup, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(covariates_for_missingness)) {
    dup <- unique(covariates_for_missingness[duplicated(covariates_for_missingness)])
    stop("Duplicate variables found in covariates_for_missingness: ", paste(dup, collapse = ", "), call. = FALSE)
  }

  if (!allow_outcome_in_missingness && any(outcome_variables %in% covariates_for_missingness)) {
    stop("Outcome variable cannot be reused as a response covariate unless explicitly allowed.", call. = FALSE)
  }

  overlap <- intersect(covariates_for_outcome, covariates_for_missingness)
  if (length(overlap) > 0 && !allow_covariate_overlap) {
    stop(
      "Covariate sets must be mutually exclusive. Overlapping variables: ",
      paste(overlap, collapse = ", ")
    )
  }

# Combine all required variables
  all_vars <- unique(c(outcome_variables, covariates_for_outcome, covariates_for_missingness))

# Check variable existence in data
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

# Check for non-finite values (Inf, -Inf, NaN) in all variables
  for (var in all_vars) {
    if (is.numeric(data[[var]])) {
      bad_indices <- which(!is.na(data[[var]]) & !is.finite(data[[var]]))
      if (length(bad_indices) > 0) {
        first_bad_idx <- bad_indices[1]
        bad_val <- data[[var]][first_bad_idx]
        var_type <- if (var %in% outcome_variables) "Outcome variable" else "Covariate"
        stop(
          var_type, " '", var, "' contains non-finite values (Inf, -Inf, or NaN).\n",
          "First non-finite value: ", bad_val, " at row ", first_bad_idx
        )
      }
    }
  }

# Validate outcome variables
  for (outcome_var in outcome_variables) {
    if (!is.numeric(data[[outcome_var]])) {
      bad_val <- data[[outcome_var]][which(!is.numeric(data[[outcome_var]]))[1]]
      stop(
        "Outcome variable '", outcome_var, "' must be numeric.\n",
        "First invalid value: '", bad_val, "' at row ", which(!is.numeric(data[[outcome_var]]))[1]
      )
    }
    if (all(is.na(data[[outcome_var]]))) {
      stop(
        "Outcome variable '", outcome_var, "' cannot be entirely NA.",
        call. = FALSE
      )
    }
  }

# Check for required NAs unless respondents-only is allowed
  if (!isTRUE(allow_respondents_only)) {
    has_missing <- vapply(outcome_variables, function(var) anyNA(data[[var]]), logical(1))
    if (!any(has_missing)) {
      stop("At least one outcome variable must contain NA values.", call. = FALSE)
    }
  }

# Validate covariates
  covariates_for_missingness_checked <- if (allow_outcome_in_missingness) {
    setdiff(covariates_for_missingness, outcome_variables)
  } else {
    covariates_for_missingness
  }

  nmar_validate_covariates(data, covariates_for_outcome, block_label = "auxiliary")
  nmar_validate_covariates(data, covariates_for_missingness_checked, block_label = "response")

  invisible(NULL)
}
