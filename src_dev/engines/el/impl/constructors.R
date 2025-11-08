#' Construct EL Result Object
#' @keywords internal
new_nmar_result_el <- function(y_hat, se, weights, coefficients, vcov,
                               converged, diagnostics, data_info,
                               nmar_scaling_recipe, fitted_values, call) {
  diagnostics <- diagnostics %||% list()
  if (is.null(data_info$method)) data_info$method <- "Empirical Likelihood (EL)"
  outcome_name <- data_info$outcome_label %||% data_info$outcome_var %||% NA_character_
  trim_fraction <- diagnostics$trimmed_fraction %||% NA_real_

  sample <- list(
# Preserve the engine-supplied analysis-scale population total (N_pop).
# Use data_info$n_total when provided (set by dataframe/survey methods), not
# data_info$nobs. This ensures weights(object, scale = "population") sums
# to N_pop and survey-scale reporting is correct.
    n_total = data_info$n_total %||% data_info$nobs %||% NA_integer_,
    n_respondents = data_info$nobs_resp %||% NA_integer_,
    is_survey = isTRUE(data_info$is_survey),
    design = if (isTRUE(data_info$is_survey)) data_info$design else NULL
  )

  df_val <- NA_real_
  if (isTRUE(sample$is_survey) && !is.null(sample$design) && requireNamespace("survey", quietly = TRUE)) {
    df_val <- tryCatch(survey::degf(sample$design), error = function(e) NA_real_)
  }

  if (!is.list(coefficients)) {
    response_coeffs <- NULL
    nuisance_terms <- list()
  } else {
    response_coeffs <- coefficients$response_model %||% NULL
    nuisance_terms <- coefficients$nuisance %||% list()
  }

  inference <- list(
    variance_method = data_info$variance_method %||% NA_character_,
    df = df_val,
    message = diagnostics$vcov_message %||% NA_character_
  )

  meta <- list(
    engine_name = "empirical_likelihood",
    call = call,
    formula = data_info$formula
  )

  result <- new_nmar_result(
    estimate = y_hat,
    estimate_name = outcome_name,
    se = se,
    converged = converged,
    model = list(
      coefficients = response_coeffs,
      vcov = vcov,
      nuisance = nuisance_terms
    ),
    weights_info = list(values = weights, trimmed_fraction = trim_fraction),
    sample = sample,
    inference = inference,
    diagnostics = diagnostics,
    meta = meta,
    extra = list(
      nuisance = nuisance_terms,
      nmar_scaling_recipe = nmar_scaling_recipe,
      fitted_values = fitted_values
    ),
    class = "nmar_result_el"
  )

  result
}

el_build_runtime_inputs <- function(data,
                                    design_info,
                                    auxiliary_means = NULL,
                                    n_total = NULL,
                                    require_na = TRUE,
                                    context = c("data.frame", "survey.design")) {
  context <- match.arg(context)
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame containing the NMAR variables.", call. = FALSE)
  }
  if (is.null(design_info$outcome_columns)) {
    stop("`design_info` must include outcome metadata from `prepare_nmar_design()`.", call. = FALSE)
  }
  outcome_col <- design_info$outcome_column %||% design_info$outcome_columns[[1]]
  if (!outcome_col %in% names(data)) {
    stop("Outcome column '", outcome_col, "' not found in supplied data.", call. = FALSE)
  }
  outcome_label <- design_info$outcome_label %||% outcome_col
  respondents_only <- !anyNA(data[[outcome_col]])
  has_aux <- length(design_info$auxiliary_vars) > 0

  if (isTRUE(require_na) && respondents_only) {
    stop(
      "Outcome variable '", outcome_label, "' must contain NA values to indicate nonresponse.",
      call. = FALSE
    )
  }
  if (respondents_only && has_aux && is.null(auxiliary_means)) {
    msg <- if (identical(context, "survey.design")) {
      "Respondents-only survey design (no NAs in outcome) and auxiliary constraints were requested, but 'auxiliary_means' was not provided. Provide population auxiliary means via auxiliary_means=."
    } else {
      "Respondents-only data detected (no NAs in outcome) and auxiliary constraints were requested, but 'auxiliary_means' was not provided. Provide population auxiliary means via auxiliary_means=."
    }
    stop(msg, call. = FALSE)
  }
  if (respondents_only && is.null(n_total)) {
    msg <- if (identical(context, "survey.design")) {
      "Respondents-only survey design detected (no NAs in outcome), but 'n_total' was not provided. Set el_engine(n_total = <total design weight or population total>)."
    } else {
      "Respondents-only data detected (no NAs in outcome), but 'n_total' was not provided. Set el_engine(n_total = <total sample size>)."
    }
    stop(msg, call. = FALSE)
  }

  delta_name <- nmar_make_unique_colname("..nmar_delta..", names(data))
  data[[delta_name]] <- as.integer(!is.na(data[[outcome_col]]))

  aux_lang_use <- if (has_aux) design_info$aux_rhs_lang else NULL
  env <- design_info$environment %||% parent.frame()
  forms <- nmar_build_internal_formulas(
    delta_name = delta_name,
    outcome_var = outcome_col,
    aux_rhs_lang = aux_lang_use,
    response_rhs_lang = design_info$response_rhs_lang,
    env = env
  )

  response_var <- all.vars(forms$response)[1]
  observed_mask <- data[[response_var]] == 1

  list(
    data = data,
    internal_formula = forms,
    response_var = response_var,
    observed_indices = which(observed_mask),
    observed_mask = observed_mask,
    outcome_name = all.vars(forms$outcome)[1],
    respondents_only = respondents_only
  )
}
