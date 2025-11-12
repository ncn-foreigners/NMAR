#' Build EL analysis context (respondent data, weights, metadata)
#'
#' Stores the augmented data, respondent slices, weights, population size, and
#' minimal metadata needed by downstream EL routines (solver, variance, result
#' builder). Both IID and survey entry points call this helper.
#'
#' @param data_aug Data with the NMAR delta column appended.
#' @param respondent_mask Logical vector identifying rows with observed outcomes.
#' @param outcome_var Character outcome name.
#' @param formula Original model formula.
#' @param full_data Underlying data object (survey designs supply the design).
#' @param respondent_weights_full Optional weight vector aligned with `data_aug`.
#' @param N_pop Optional population size; defaults to `nrow(data_aug)` if `NULL`.
#' @param is_survey Logical; TRUE when called from the survey method.
#' @param design Survey design object when `is_survey = TRUE`.
#' @param variance_method Character variance method label to carry through logging.
#' @keywords internal
el_build_context <- function(data_aug,
                             respondent_mask,
                             outcome_var,
                             formula,
                             full_data = NULL,
                             respondent_weights_full = NULL,
                             N_pop = NULL,
                             is_survey = FALSE,
                             design = NULL,
                             variance_method) {
  respondent_indices <- which(respondent_mask)
  if (length(respondent_indices) == 0) {
    stop("No respondents detected in data after preprocessing.", call. = FALSE)
  }

  if (is.null(respondent_weights_full)) {
    respondent_weights <- rep(1, length(respondent_indices))
  } else {
    if (length(respondent_weights_full) != nrow(data_aug)) {
      stop("`respondent_weights` must have the same length as the number of rows in the data.", call. = FALSE)
    }
    respondent_weights <- respondent_weights_full[respondent_indices]
  }

  N_pop_val <- N_pop %||% nrow(data_aug)

  context <- list(
    data = data_aug,
    full_data = full_data %||% design %||% data_aug,
    respondent_data = data_aug[respondent_indices, , drop = FALSE],
    respondent_weights = respondent_weights,
    respondent_indices = respondent_indices,
    N_pop = N_pop_val,
    outcome_var = outcome_var,
    formula = formula,
    variance_method = variance_method,
    is_survey = is_survey,
    design = if (is_survey) design else NULL
  )
  class(context) <- "el_context"
  context
}
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
