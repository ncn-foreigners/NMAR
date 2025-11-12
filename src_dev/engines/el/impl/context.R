#' Build EL analysis context (respondent data, weights, metadata)
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
    N_pop = N_pop_val,
    meta = list(
      outcome_var = outcome_var,
      nobs_resp = length(respondent_indices),
      n_total = N_pop_val,
      is_survey = is_survey,
      design = if (is_survey) design else NULL,
      variance_method = variance_method,
      call = NULL,
      formula = formula
    )
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
