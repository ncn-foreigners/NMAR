#' Build EL analysis context (weights, population size, metadata)
#' @keywords internal
el_build_context <- function(prepared_inputs,
                             full_data,
                             formula,
                             respondent_weights_full = NULL,
                             N_pop = NULL,
                             is_survey = FALSE,
                             design = NULL,
                             variance_method) {
  data_aug <- prepared_inputs$data
  respondent_indices <- which(prepared_inputs$respondent_mask)
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

  analysis_info <- list(
    outcome_var = prepared_inputs$outcome_var,
    nobs_resp = length(respondent_indices),
    n_total = N_pop_val,
    is_survey = is_survey,
    design = if (is_survey) design else NULL,
    variance_method = variance_method
  )

  list(
    full_data = full_data,
    respondent_data = data_aug[respondent_indices, , drop = FALSE],
    respondent_weights = respondent_weights,
    N_pop = N_pop_val,
    analysis_info = analysis_info,
    respondent_indices = respondent_indices,
    formula = formula
  )
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
