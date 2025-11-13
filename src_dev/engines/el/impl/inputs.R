#' Validate respondents-only workflows
#'
#' Ensures that respondents-only data (no NA values in the outcome) satisfy the
#' ancillary requirements of the empirical likelihood estimator: population
#' totals (`n_total`) must be supplied elsewhere and auxiliary constraints
#' require `auxiliary_means`.
#'
#' @param formula Two-sided formula supplied to `nmar()`.
#' @param data Data frame or survey variables referenced by the formula.
#' @param auxiliary_means Optional named vector of auxiliary population means.
#' @param context_label Character string inserted into error messages to clarify the source (e.g., "data frame").
#' @return Logical; `TRUE` when the outcome contains no NA values.
#' @keywords internal
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
#' @param data Data frame (or survey design variables) containing the outcome.
#' @param outcome_var Character outcome column name.
#' @param mask Logical vector identifying respondents (non-missing outcomes).
#' @param weights_full Optional weight vector aligned with `data` (survey designs supply this).
#' @param N_pop Optional population size; defaults to the augmented row count when missing.
#' @param variance_method Character variance label to store in result metadata.
#' @param is_survey Logical; `TRUE` when called from the survey method.
#' @param design Survey design object when `is_survey = TRUE`.
#' @return A list with fields `data_aug`, `respondent_weights`, `respondent_indices`, `N_pop`, and `data_info`.
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

#' Create the NMAR delta indicator column
#' @keywords internal
el_make_delta_column <- function(data, outcome_var, respondent_mask = NULL) {
  if (is.null(respondent_mask)) {
    respondent_mask <- !is.na(data[[outcome_var]])
  }
  if (length(respondent_mask) != nrow(data)) {
    stop("Internal error: respondent mask must align with data.", call. = FALSE)
  }

  delta_name <- "..nmar_delta.."
  if (delta_name %in% names(data)) {
    i <- 1L
    while (paste0(delta_name, i) %in% names(data)) i <- i + 1L
    delta_name <- paste0(delta_name, i)
  }

  data[[delta_name]] <- as.integer(respondent_mask)
  list(data = data, delta_name = delta_name)
}

#' Launch EL estimation once design matrices are parsed
#'
#' Wraps the shared respondent-level preparation, auxiliary resolution,
#' and downstream estimator call so that IID and survey entry points stay in sync.
#' @keywords internal
el_run_core_analysis <- function(call,
                                 formula,
                                 raw_data,
                                 design_inputs,
                                 weights_full,
                                 n_total,
                                 variance_method,
                                 is_survey,
                                 design_object,
                                 auxiliary_means,
                                 standardize,
                                 trim_cap,
                                 control,
                                 on_failure,
                                 family,
                                 bootstrap_reps,
                                 start,
                                 trace_level,
                                 extra_user_args = list()) {
  prep <- el_prepare_analysis_inputs(
    data = raw_data,
    outcome_var = design_inputs$outcome_var,
    mask = design_inputs$respondent_mask,
    weights_full = weights_full,
    N_pop = n_total,
    variance_method = variance_method,
    is_survey = is_survey,
    design = if (is_survey) design_object else NULL
  )

  analysis_object <- if (is_survey) {
    design_object$variables <- prep$data_aug
    design_object
  } else {
    prep$data_aug
  }

  aux_summary <- el_resolve_auxiliaries(
    design_inputs$aux_full[design_inputs$respondent_mask, , drop = FALSE],
    design_inputs$aux_full,
    auxiliary_means,
    weights_full = weights_full
  )

  user_args <- c(
    list(
      formula = formula,
      auxiliary_means = auxiliary_means,
      standardize = standardize,
      trim_cap = trim_cap,
      control = control,
      n_total = n_total
    ),
    extra_user_args
  )

  core_results <- el_estimator_core(
    response_matrix = design_inputs$response_matrix,
    response_outcome = design_inputs$y_obs,
    auxiliary_matrix = aux_summary$matrix,
    mu_x = aux_summary$means,
    respondent_weights = prep$respondent_weights,
    full_data = analysis_object,
    outcome_var = design_inputs$outcome_var,
    N_pop = prep$N_pop,
    standardize = standardize,
    trim_cap = trim_cap,
    control = control,
    on_failure = on_failure,
    family = family,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    user_args = user_args,
    start = start,
    trace_level = trace_level,
    auxiliary_means = auxiliary_means
  )

  el_build_result(core_results, prep$data_info, call, formula)
}
