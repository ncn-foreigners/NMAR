#' Construct EL Result Object
#' @keywords internal
# new_nmar_result_el <- function(y_hat, se, weights, coefficients, vcov,
#                                converged, diagnostics, data_info,
#                                nmar_scaling_recipe, fitted_values, call) {
#   if (is.null(data_info$method)) data_info$method <- "Empirical Likelihood (EL)"
#   # Standardize diagnostics and data_info bundles
#   diagnostics <- new_nmar_diagnostics(diagnostics)
#   data_info <- new_nmar_data_info(data_info)
#   structure(
#     list(
#       y_hat = y_hat, se = se, weights = weights, coefficients = coefficients,
#       vcov = vcov, converged = converged, diagnostics = diagnostics,
#       data_info = data_info, nmar_scaling_recipe = nmar_scaling_recipe,
#       fitted_values = fitted_values, call = call
#     ),
#     class = c("nmar_result_el", "nmar_result")
#   )
# }
#' Construct EL Result Object
#' @keywords internal
new_nmar_result_el <- function(y_hat, se, weights, coefficients, vcov,
                               converged, diagnostics, data_info,
                               nmar_scaling_recipe, fitted_values, call) {


  if (is.null(data_info$method)) {
    data_info$method <- "Empirical Likelihood (EL)"
  }



  result <- new_nmar_result(
    y_hat = y_hat,
    se = se,
    weights = weights,
    coefficients = coefficients,
    vcov = vcov,
    converged = converged,
    class = "nmar_result_el"
  )


  result$diagnostics <- diagnostics
  result$data_info <- data_info
  result$fitted_values <- fitted_values
  result$call <- call
  result$nmar_scaling_recipe <- nmar_scaling_recipe

  result$diagnostics <- new_nmar_diagnostics(diagnostics)
  result$data_info <- new_nmar_data_info(data_info)

  return(result)
}

#' Prepare inputs for EL estimation
#' @details Validates the two-sided outcome formula and constructs three
#'   internal formulas: outcome (~ outcome_var), response (for the missingness
#'   model), and auxiliary (RHS only, no intercept). The `response_predictors`
#'   may include variables not on the outcome RHS; such variables enter only the
#'   response model (no auxiliary moment constraint). Only variables on the
#'   outcome RHS are treated as auxiliaries and, when provided, must match the
#'   names in `auxiliary_means`. See Qin, Leung and Shao (2002) for the EL
#'   formulation.
#' @keywords internal
prepare_el_inputs <- function(formula, data, response_predictors) {
  if (!inherits(formula, "formula") || length(formula) != 3 || length(all.vars(formula[[2]])) != 1) {
    stop("`formula` must be a two-sided formula with a single variable on the LHS, e.g., y ~ x1 + x2.", call. = FALSE)
  }
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  outcome_var <- all.vars(formula[[2]])
  rhs_vars <- all.vars(formula[[3]])
  if (!outcome_var %in% names(data)) stop(sprintf("Outcome variable '%s' not found in the data.", outcome_var), call. = FALSE)
  if (!anyNA(data[[outcome_var]])) stop(sprintf("Outcome variable '%s' must contain NA values to indicate nonresponse.", outcome_var), call. = FALSE)
  response_predictors_full <- c(outcome_var, response_predictors)
  if (!is.null(response_predictors)) {
    if (!is.character(response_predictors)) stop("`response_predictors` must be a character vector of variable names.", call. = FALSE)
    missing_in_data <- setdiff(response_predictors, names(data))
    if (length(missing_in_data) > 0) stop(sprintf("Variables in `response_predictors` not found in data: %s", paste(missing_in_data, collapse = ", ")), call. = FALSE)
  }
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

#' Validator for EL result
#' @keywords internal
# validate_nmar_result_el <- function(x) {
#   stopifnot(is.list(x), inherits(x, "nmar_result_el"))
#   if (isTRUE(x$converged)) {
#     stopifnot(is.finite(x$y_hat), is.numeric(x$se))
#   }
#   x
# }
