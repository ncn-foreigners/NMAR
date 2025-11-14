#' Enforce respondents-only requirements for a given design
#' @keywords internal
el_check_respondents_only_requirements <- function(design, n_total, auxiliary_means, context_label) {
  if (!isTRUE(all(design$respondent_mask))) return(invisible(NULL))
  noun <- if (identical(context_label, "data frame")) "data" else context_label
  message_prefix <- sprintf("Respondents-only %s detected (no NAs in outcome)", noun)
  if (is.null(n_total)) {
    stop(
      sprintf("%s, but 'n_total' was not provided.", message_prefix),
      call. = FALSE
    )
  }
  has_aux <- is.matrix(design$auxiliary_design_full) && ncol(design$auxiliary_design_full) > 0
  if (has_aux && is.null(auxiliary_means)) {
    stop(
      sprintf(
        "%s and auxiliary constraints were requested, but 'auxiliary_means' was not provided. Provide population auxiliary means via auxiliary_means=.",
        message_prefix
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Prepare respondent-level inputs shared by IID and survey entry points
#'
#' Builds the analysis object and respondent metadata used by the estimator.
#' Adds a delta indicator column to the data and returns respondent weights,
#' indices, and the analysis-scale population size.
#'
#' @param data Data frame (or survey design variables) containing the outcome.
#' @param design_inputs Design list returned by `el_prepare_design()`.
#' @param weights_full Optional weight vector aligned with `data` (survey designs supply this).
#' @param N_pop Optional population size; defaults to the augmented row count when missing.
#' @param is_survey Logical; `TRUE` when called from the survey method.
#' @param design_object Survey design object when `is_survey = TRUE`.
#' @return A list with fields `analysis_object`, `respondent_weights`, `respondent_indices`, and `N_pop`.
#' @keywords internal
el_prepare_analysis_context <- function(data,
                                        design_inputs,
                                        weights_full = NULL,
                                        N_pop = NULL,
                                        is_survey = FALSE,
                                        design_object = NULL) {
  outcome_var <- design_inputs$outcome_var
  mask <- design_inputs$respondent_mask
  if (length(mask) != nrow(data)) {
    stop("Internal error: respondent mask must have the same length as data.", call. = FALSE)
  }

  delta <- el_make_delta_column_name(data, outcome_var, mask)
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

  analysis_object <- if (isTRUE(is_survey)) {
    design_object$variables <- data_aug
    design_object
  } else {
    data_aug
  }

  list(
    analysis_object = analysis_object,
    respondent_weights = respondent_weights,
    respondent_indices = respondent_indices,
    N_pop = N_pop_val
  )
}

#' Create the NMAR delta indicator column
#' @keywords internal
el_make_delta_column_name <- function(data, outcome_var, respondent_mask = NULL) {
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
  list(data = data, delta_column_name = delta_name)
}

#' Launch EL estimation once design matrices are parsed
#'
#' Runs the shared respondent prep, auxiliary resolution, and solver call so IID
#' and survey methods stay aligned.
#' @keywords internal
el_run_core_analysis <- function(call,
                                 formula,
                                 raw_data,
                                 design,
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
  el_validate_design_spec(
    design = design,
    data_nrow = nrow(raw_data),
    context_label = if (isTRUE(is_survey)) "survey design" else "data"
  )
  context_label <- if (isTRUE(is_survey)) "survey design" else "data"
  el_check_respondents_only_requirements(
    design = design,
    n_total = n_total,
    auxiliary_means = auxiliary_means,
    context_label = context_label
  )
  prep <- el_prepare_analysis_context(
    data = raw_data,
    design_inputs = design,
    weights_full = weights_full,
    N_pop = n_total,
    is_survey = is_survey,
    design_object = if (is_survey) design_object else NULL
  )
  analysis_object <- prep$analysis_object

  aux_summary <- el_resolve_auxiliaries(
    auxiliary_design_full = design$auxiliary_design_full,
    respondent_mask = design$respondent_mask,
    auxiliary_means = auxiliary_means,
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
    missingness_design = design$missingness_design,
    auxiliary_matrix = aux_summary$auxiliary_design,
    mu_x = aux_summary$means,
    respondent_weights = prep$respondent_weights,
    full_data = analysis_object,
    outcome_var = design$outcome_var,
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

  data_info <- list(
    outcome_var = design$outcome_var,
    n_total = prep$N_pop,
    nobs_resp = length(prep$respondent_indices),
    is_survey = is_survey,
    design = if (is_survey) analysis_object else NULL,
    variance_method = variance_method
  )
  el_build_result(core_results, data_info, call, formula)
}
