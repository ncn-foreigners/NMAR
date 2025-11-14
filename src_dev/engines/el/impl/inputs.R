#' Build the combined EL input specification
#'
#' Constructs the parsed design matrices and augments the data with the
#' respondent indicator so iid and survey entry points can share the same
#' downstream workflow.
#' @keywords internal
el_build_input_spec <- function(formula,
                                data,
                                weights_full = NULL,
                                population_total = NULL,
                                n_total_arg = population_total,
                                is_survey = FALSE,
                                design_object = NULL,
                                auxiliary_means = NULL) {
  context_label <- if (isTRUE(is_survey)) "survey design" else "data frame"
  design <- el_prepare_design(
    formula = formula,
    data = data,
    require_na = FALSE
  )
  el_validate_design_spec(design, data_nrow = nrow(data), context_label = context_label)
  el_require_population_inputs(design, n_total_arg, auxiliary_means, context_label)

  outcome_var <- design$outcome_source %||% design$outcome_var
  mask <- design$respondent_mask
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

  analysis_object <- if (isTRUE(is_survey)) {
    if (is.null(design_object)) {
      stop("Internal error: design object must be supplied for survey workflows.", call. = FALSE)
    }
    design_object$variables <- data_aug
    design_object
  } else {
    data_aug
  }

  N_pop_val <- population_total %||% nrow(data_aug)

  list(
    missingness_design = design$missingness_design,
    auxiliary_design_full = design$auxiliary_design_full,
    respondent_mask = mask,
    outcome_var = design$outcome_var,
    analysis_object = analysis_object,
    respondent_weights = respondent_weights,
    respondent_indices = respondent_indices,
    N_pop = N_pop_val,
    weights_full = weights_full,
    is_survey = isTRUE(is_survey)
  )
}

#' Enforce respondents-only requirements for a given design
#' @keywords internal
el_require_population_inputs <- function(design, n_total, auxiliary_means, context_label) {
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
                                 input_spec,
                                 variance_method,
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
  aux_summary <- el_resolve_auxiliaries(
    auxiliary_design_full = input_spec$auxiliary_design_full,
    respondent_mask = input_spec$respondent_mask,
    auxiliary_means = auxiliary_means,
    weights_full = input_spec$weights_full
  )

  user_args <- c(
    list(
      formula = formula,
      auxiliary_means = auxiliary_means,
      standardize = standardize,
      trim_cap = trim_cap,
      control = control,
      n_total = input_spec$N_pop
    ),
    extra_user_args
  )

  core_results <- el_estimator_core(
    missingness_design = input_spec$missingness_design,
    auxiliary_matrix = aux_summary$auxiliary_design,
    mu_x = aux_summary$means,
    respondent_weights = input_spec$respondent_weights,
    full_data = input_spec$analysis_object,
    outcome_var = input_spec$outcome_var,
    N_pop = input_spec$N_pop,
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

  input_spec$variance_method <- variance_method
  el_build_result(core_results, input_spec, call, formula)
}
