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
#' Build internal EL formulas and data
#' @keywords internal
el_prepare_inputs <- function(formula, data, require_na = TRUE) {
  if (!inherits(formula, "formula") || length(formula) != 3 || length(all.vars(formula[[2]])) != 1) {
    stop("`formula` must be a two-sided formula with a single variable on the LHS, e.g., y ~ x1 + x2.", call. = FALSE)
  }
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()

# LHS outcome
  outcome_sym <- formula[[2L]]
  outcome_var <- all.vars(outcome_sym)[1]

# Split RHS by `|` (preserve language objects)
  rhs <- formula[[3L]]
  aux_expr <- rhs
  resp_expr <- NULL
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    aux_expr <- rhs[[2L]]
    resp_expr <- rhs[[3L]]
  }

# Validation of referenced variables (use symbol sets)
  aux_vars <- unique(all.vars(aux_expr))
  resp_vars <- if (is.null(resp_expr)) character() else unique(all.vars(resp_expr))
  if (!outcome_var %in% names(data)) {
    stop(sprintf("Variables not found in data: %s", outcome_var), call. = FALSE)
  }
  missing_in_aux <- setdiff(aux_vars, names(data))
  if (length(missing_in_aux) > 0) {
    stop(sprintf("Variables not found in data: %s", paste(missing_in_aux, collapse = ", ")), call. = FALSE)
  }
  if (isTRUE(require_na) && !anyNA(data[[outcome_var]])) stop(sprintf("Outcome variable '%s' must contain NA values to indicate nonresponse.", outcome_var), call. = FALSE)
  missing_in_resp <- setdiff(resp_vars, names(data))
  if (length(missing_in_resp) > 0) stop(sprintf("Variables in response part not found in data: %s", paste(missing_in_resp, collapse = ", ")), call. = FALSE)

# Create response indicator column name (avoid collision)
  delta_name <- "..nmar_delta.."
  if (delta_name %in% names(data)) {
    i <- 1L
    while (paste0(delta_name, i) %in% names(data)) i <- i + 1L
    delta_name <- paste0(delta_name, i)
  }

# Build internal outcome and response formulas from language objects
  data2 <- data
  data2[[delta_name]] <- as.integer(!is.na(data2[[outcome_var]]))

# outcome: y ~ 1
  outcome_fml <- as.formula(call("~", outcome_sym, 1L))
  environment(outcome_fml) <- env

# response: delta ~ (outcome + resp_expr?)
  rhs_resp <- if (is.null(resp_expr)) {
    outcome_sym
  } else {
    call("+", outcome_sym, resp_expr)
  }
  response_fml <- as.formula(call("~", as.name(delta_name), rhs_resp))
  environment(response_fml) <- env

# auxiliary: ~ 0 + aux_expr (no intercept) without using update(. ~ . - 1)
  auxiliary_fml <- NULL
  if (!is.null(aux_expr) && length(all.vars(aux_expr)) > 0) {
    rhs <- call("+", 0, aux_expr)
    auxiliary_fml <- as.formula(call("~", rhs))
    environment(auxiliary_fml) <- env
  }

  list(data = data2, formula_list = list(outcome = outcome_fml, response = response_fml, auxiliary = auxiliary_fml))
}

# Backward-compatible alias (kept silent to avoid test warnings). Prefer el_prepare_inputs().
prepare_el_inputs <- function(formula, data, require_na = TRUE) {
  el_prepare_inputs(formula = formula, data = data, require_na = require_na)
}
