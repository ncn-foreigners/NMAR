#' Construct EL Result Object
#' @keywords internal
new_nmar_result_el <- function(y_hat, se, weights, coefficients, vcov,
                               converged, diagnostics, data_info,
                               nmar_scaling_recipe, fitted_values, call) {
  diagnostics <- diagnostics %||% list()
  if (is.null(data_info$method)) data_info$method <- "Empirical Likelihood (EL)"
  outcome_name <- data_info$outcome_var %||% NA_character_
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

#' Prepare inputs for EL estimation
#' @details Validates the two-sided outcome formula and constructs three
#'   internal formulas: outcome (~ outcome_var), response (for the missingness
#'   model), and auxiliary (RHS only, no intercept). Response-only predictors
#'   may include variables not on the outcome RHS; such variables enter only the
#'   response model (no auxiliary moment constraint). Only variables on the
#'   outcome RHS are treated as auxiliaries and, when provided, must match the
#'   names in `auxiliary_means`. See Qin, Leung and Shao (2002) for the EL
#'   formulation.
#' @keywords internal
prepare_el_inputs <- function(formula, data, require_na = TRUE) {
  if (!inherits(formula, "formula") || length(formula) != 3 || length(all.vars(formula[[2]])) != 1) {
    stop("`formula` must be a two-sided formula with a single variable on the LHS, e.g., y ~ x1 + x2.", call. = FALSE)
  }
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  outcome_var <- all.vars(formula[[2]])

# Split RHS by `|` if present: auxiliaries on the left, response-only on the right
  rhs <- formula[[3]]
  aux_expr <- rhs
  resp_expr <- NULL
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    aux_expr <- rhs[[2L]]
    resp_expr <- rhs[[3L]]
  }
  rhs_vars <- all.vars(aux_expr)
  response_predictors <- if (is.null(resp_expr)) character() else unique(all.vars(resp_expr))
  if (!outcome_var %in% names(data)) stop(sprintf("Outcome variable '%s' not found in the data.", outcome_var), call. = FALSE)
  if (isTRUE(require_na) && !anyNA(data[[outcome_var]])) stop(sprintf("Outcome variable '%s' must contain NA values to indicate nonresponse.", outcome_var), call. = FALSE)
  response_predictors_full <- c(outcome_var, response_predictors)
  missing_in_data <- setdiff(response_predictors, names(data))
  if (length(missing_in_data) > 0) stop(sprintf("Variables in response part not found in data: %s", paste(missing_in_data, collapse = ", ")), call. = FALSE)
  delta_name <- "..nmar_delta.."
  if (delta_name %in% names(data)) {
    i <- 1L
    while (paste0(delta_name, i) %in% names(data)) i <- i + 1L
    delta_name <- paste0(delta_name, i)
  }
  data2 <- data
  data2[[delta_name]] <- as.integer(!is.na(data2[[outcome_var]]))
  outcome_fml <- stats::as.formula(paste(outcome_var, "~ 1"))
  environment(outcome_fml) <- env
  response_rhs <- if (length(response_predictors_full) > 0) paste(response_predictors_full, collapse = " + ") else "1"
  response_fml <- stats::as.formula(paste(delta_name, "~", response_rhs))
  environment(response_fml) <- env
  auxiliary_fml <- if (length(rhs_vars) > 0) {
    f <- stats::as.formula(paste("~", paste(rhs_vars, collapse = " + "), "- 1"))
    environment(f) <- env
    f
  } else {
    NULL
  }
  list(data = data2, formula_list = list(outcome = outcome_fml, response = response_fml, auxiliary = auxiliary_fml))
}
