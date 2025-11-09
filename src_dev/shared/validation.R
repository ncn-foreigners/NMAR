#' Validation utilities for NMAR input
#'
#' These functions enforce basic structural contracts between outcomes,
#' auxiliary predictors, and explicit response-only predictors. Importantly,
#' the outcome is always included implicitly in the response model for NMAR;
#' validation only concerns explicit RHS variables supplied by the user.
#'
#' @keywords internal

validate_nmar_args <- function(spec, traits = list()) {
  if (!inherits(spec, "nmar_input_spec")) {
    stop("`spec` must be created by `parse_nmar_spec()`.", call. = FALSE)
  }
  traits <- utils::modifyList(NMAR_DEFAULT_TRAITS, traits)

  primary_outcome <- spec$outcome_primary %||% spec$outcome[[1]]
  if (traits$requires_single_outcome && (is.null(primary_outcome) || is.na(primary_outcome))) {
    stop("The formula must have exactly one outcome variable on the left-hand side.", call. = FALSE)
  }
  if (traits$requires_single_outcome && length(spec$outcome) != 1L) {
    stop("The formula must have exactly one outcome variable on the left-hand side.", call. = FALSE)
  }

  aux_vars_canonical <- unique(spec$auxiliary_vars %||% character())
  aux_vars_raw <- unique(spec$auxiliary_vars_raw %||% character())
  response_vars_canonical <- unique(spec$response_predictors %||% character())
  response_vars_raw <- unique(spec$response_predictors_raw %||% character())

  aux_vars_check <- if (length(aux_vars_canonical)) aux_vars_canonical else aux_vars_raw
# Enforce engine policy using the design-expanded response predictors; fall back to
# raw symbols only if blueprint information is unavailable.
  response_vars_check <- if (length(response_vars_canonical)) response_vars_canonical else response_vars_raw

  outcomes_for_relationships <- if (length(spec$outcome)) spec$outcome else primary_outcome
  validate_predictor_relationships(
    outcomes = outcomes_for_relationships,
    auxiliary_vars = aux_vars_check,
    response_vars = response_vars_check,
    allow_outcome_in_missingness = traits$allow_outcome_in_missingness,
    allow_covariate_overlap = traits$allow_covariate_overlap
  )

  if (length(spec$outcome) == 1L) {
    validate_data(
      data = spec$original_data,
      outcome_variable = spec$outcome[[1]],
      covariates_for_outcome = aux_vars_check,
      covariates_for_missingness = response_vars_check,
      allow_outcome_in_missingness = traits$allow_outcome_in_missingness,
      allow_covariate_overlap = traits$allow_covariate_overlap,
      allow_respondents_only = isTRUE(traits$allow_respondents_only)
    )
  } else {
    validate_multi_outcome_data(spec$data, spec$outcome)
    nmar_validate_covariates(spec$data, aux_vars_check, block_label = "auxiliary")
    nmar_validate_covariates(spec$data, response_vars_check, block_label = "response")
  }

  invisible(spec)
}

#' @keywords internal
validate_predictor_relationships <- function(outcomes,
                                             auxiliary_vars,
                                             response_vars,
                                             allow_outcome_in_missingness,
                                             allow_covariate_overlap) {
  auxiliary_vars <- unique(auxiliary_vars %||% character())
  response_vars <- unique(response_vars %||% character())
  if (!allow_outcome_in_missingness) {
    outcome_on_rhs <- intersect(outcomes, response_vars)
    if (length(outcome_on_rhs) > 0) {
      stop(
        "Outcome variable cannot appear on the response-model side unless the engine allows it.\n",
        "Variables: ", paste(outcome_on_rhs, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!allow_covariate_overlap) {
    overlap <- intersect(auxiliary_vars, response_vars)
    if (length(overlap) > 0) {
      stop(
        "Covariate sets must be mutually exclusive. Overlapping variables: ",
        paste(overlap, collapse = ", "),
        call. = FALSE
      )
    }
  }
}

#' @keywords internal
nmar_validate_covariates <- function(data, vars, block_label) {
  vars <- unique(vars)
  if (!length(vars)) return(invisible(NULL))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      block_label, " covariates not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  for (var in vars) {
    column <- data[[var]]
    valid_type <- is.numeric(column) || is.logical(column) || is.factor(column) || is.character(column)
    if (!valid_type) {
      stop(
        "Covariate '", var, "' in ", block_label, " predictors must be numeric, logical, factor, or character.\n",
        "Detected type: ", paste(class(column), collapse = "/"),
        call. = FALSE
      )
    }
    if (is.numeric(column)) {
      bad_indices <- which(!is.na(column) & !is.finite(column))
      if (length(bad_indices) > 0) {
        first_bad_idx <- bad_indices[1]
        bad_val <- column[first_bad_idx]
        stop(
          block_label, " covariate '", var, "' contains non-finite values (Inf, -Inf, or NaN).\n",
          "First non-finite value: ", bad_val, " at row ", first_bad_idx,
          call. = FALSE
        )
      }
    }
    if (anyNA(column)) {
      stop(
        block_label, " covariate '", var, "' contains NA values.\n",
        "First NA at row ", which(is.na(column))[1],
        call. = FALSE
      )
    }
  }
}

#' @keywords internal
validate_multi_outcome_data <- function(data, outcome_vars) {
  missing_outcomes <- setdiff(outcome_vars, names(data))
  if (length(missing_outcomes) > 0) {
    stop("Outcome variables not found in data: ", paste(missing_outcomes, collapse = ", "))
  }
  for (outcome_var in outcome_vars) {
    col <- data[[outcome_var]]
    if (!is.numeric(col)) {
      bad_val <- col[which(!is.numeric(col))[1]]
      stop(
        "Outcome variable '", outcome_var, "' must be numeric.\n",
        "First invalid value: '", bad_val, "' at row ", which(!is.numeric(col))[1]
      )
    }
    if (all(is.na(col))) {
      stop(
        "Outcome variable '", outcome_var, "' cannot be entirely NA.",
        call. = FALSE
      )
    }
  }
}
