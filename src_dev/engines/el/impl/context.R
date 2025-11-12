el_validate_respondents_only <- function(formula, data, auxiliary_means, context_label = "data") {
  outcome_var <- all.vars(formula[[2L]])
  respondents_only <- length(outcome_var) == 1 && !anyNA(data[[outcome_var]])

  has_aux <- {
    rhs <- formula[[3L]]
    aux_expr <- if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) rhs[[2L]] else rhs
    length(all.vars(aux_expr)) > 0
  }

  if (respondents_only && has_aux && is.null(auxiliary_means)) {
    stop(
      paste0(
        "Respondents-only ", context_label, " detected (no NAs in outcome) and auxiliary constraints were requested, ",
        "but 'auxiliary_means' was not provided. Provide population auxiliary means via auxiliary_means=."
      ),
      call. = FALSE
    )
  }

  respondents_only
}

#' Prepare respondent-level inputs shared by IID and survey entry points
#'
#' Attaches the NMAR delta indicator, slices respondent rows/weights, and
#' prepares the metadata needed by downstream estimators and result builders.
#'
#' @param data Data frame (or survey design variables) containing the outcome.
#' @param outcome_var Character outcome column name.
#' @param mask Logical vector identifying respondents (non-missing outcomes).
#' @param weights_full Optional weight vector aligned with `data` (survey designs supply this).
#' @param N_pop Optional population size; defaults to the augmented row count when missing.
#' @param variance_method Character variance label to store in result metadata.
#' @param is_survey Logical; `TRUE` when called from the survey method.
#' @param design Survey design object when `is_survey = TRUE`.
#'
#' @return A list with fields `data_aug`, `respondent_weights`, `respondent_indices`,
#'   `N_pop`, and `data_info`.
#' @keywords internal
el_prepare_analysis_inputs <- function(data,
                                       outcome_var,
                                       mask,
                                       weights_full = NULL,
                                       N_pop = NULL,
                                       variance_method,
                                       is_survey = FALSE,
                                       design = NULL) {
  if (length(mask) != nrow(data)) {
    stop("Internal error: respondent mask must have the same length as data.", call. = FALSE)
  }

  delta <- el_make_delta_column(data, outcome_var, mask)
  data_aug <- delta$data
  respondent_indices <- which(mask)
  if (length(respondent_indices) == 0) {
    stop("No respondents detected in data after preprocessing.", call. = FALSE)
  }

  if (!is.null(weights_full)) {
    if (length(weights_full) != nrow(data)) {
      stop("`weights_full` must align with the number of rows in the data.", call. = FALSE)
    }
    respondent_weights <- weights_full[mask]
  } else {
    respondent_weights <- rep(1, length(respondent_indices))
  }

  N_pop_val <- N_pop %||% nrow(data_aug)

  data_info <- list(
    outcome_var = outcome_var,
    n_total = N_pop_val,
    nobs_resp = length(respondent_indices),
    is_survey = is_survey,
    design = if (is_survey) design else NULL,
    variance_method = variance_method
  )

  list(
    data_aug = data_aug,
    respondent_weights = respondent_weights,
    respondent_indices = respondent_indices,
    N_pop = N_pop_val,
    data_info = data_info
  )
}
