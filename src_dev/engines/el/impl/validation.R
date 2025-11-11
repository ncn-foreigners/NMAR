#' EL input validation and parsing (engine-local)
#'
#' Provides EL-specific input validation and constructs internal formulas for
#' outcome, response, and auxiliary models. Centralizing these checks clarifies
#' responsibilities and keeps constructors lean.
#' @keywords internal
NULL

## ---- Helpers ----

#' @keywords internal
el_split_rhs <- function(formula) {
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  rhs <- formula[[3L]]
  aux_expr <- rhs
  resp_expr <- NULL
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    aux_expr <- rhs[[2L]]
    resp_expr <- rhs[[3L]]
  }
  list(aux_expr = aux_expr, resp_expr = resp_expr, env = env)
}

#' @keywords internal
el_validate_outcome <- function(data, outcome_var, require_na) {
  if (!outcome_var %in% names(data)) {
    stop(sprintf("Variables not found in data: %s", outcome_var), call. = FALSE)
  }
  if (!is.numeric(data[[outcome_var]])) {
    stop(sprintf("Outcome variable '%s' must be numeric.", outcome_var), call. = FALSE)
  }
  if (isTRUE(require_na) && !anyNA(data[[outcome_var]])) {
    stop(sprintf("Outcome variable '%s' must contain NA values to indicate nonresponse.", outcome_var), call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
el_validate_vars_present <- function(data, aux_vars, resp_vars) {
  missing_in_aux <- setdiff(aux_vars, names(data))
  if (length(missing_in_aux) > 0) {
    stop(sprintf("Variables not found in data: %s", paste(missing_in_aux, collapse = ", ")), call. = FALSE)
  }
  missing_in_resp <- setdiff(resp_vars, names(data))
  if (length(missing_in_resp) > 0) {
    stop(sprintf("Variables not found in data: %s", paste(missing_in_resp, collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
el_validate_aux_no_na <- function(data, aux_expr, respondent_mask, auxiliary_means) {
  if (is.null(aux_expr) || length(all.vars(aux_expr)) == 0) return(invisible(NULL))
  aux_mf <- tryCatch(
    model.frame(
      as.formula(call("~", call("+", 0, aux_expr))),
      data = data,
      na.action = stats::na.pass
    ),
    error = function(e) NULL
  )
  if (!is.null(aux_mf) && ncol(aux_mf) > 0) {
# Always forbid NA among respondents (aux used in constraints for respondents)
    if (!is.null(respondent_mask) && anyNA(respondent_mask)) {
      stop("Internal error: respondent mask contains NA.", call. = FALSE)
    }
    if (!is.null(respondent_mask) && length(respondent_mask) == nrow(aux_mf)) {
      aux_resp <- aux_mf[respondent_mask, , drop = FALSE]
      if (anyNA(aux_resp)) {
        na_by_col <- vapply(seq_len(ncol(aux_resp)), function(j) anyNA(aux_resp[[j]]), logical(1))
        bad_col <- colnames(aux_resp)[which(na_by_col)[1]]
        bad_row_local <- which(is.na(aux_resp[[bad_col]]))[1]
# Map back to global row index for clarity
        bad_row <- which(respondent_mask)[bad_row_local]
        msg <- sprintf("Covariate '%s' contains NA values among respondents.%s",
                       bad_col,
                       if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else "")
        stop(msg, call. = FALSE)
      }
    }
# If auxiliary means are not provided, we must compute μ_x from full data => forbid NA anywhere
    if (is.null(auxiliary_means) && anyNA(aux_mf)) {
      na_by_col <- vapply(seq_len(ncol(aux_mf)), function(j) anyNA(aux_mf[[j]]), logical(1))
      bad_col <- colnames(aux_mf)[which(na_by_col)[1]]
      bad_row <- which(is.na(aux_mf[[bad_col]]))[1]
      msg <- sprintf("Auxiliary variables contain NA values in full data (needed to compute population means).%s\n",
                     if (is.finite(bad_row)) sprintf("\nFirst NA in '%s' at row %d", bad_col, bad_row) else "")
      msg <- paste0(msg, "Either provide auxiliary_means=... or remove NA values from auxiliary variables.")
      stop(msg, call. = FALSE)
    }
  }
  invisible(NULL)
}

#' @keywords internal
el_make_delta_column <- function(data, outcome_var) {
  delta_name <- "..nmar_delta.."
  if (delta_name %in% names(data)) {
    i <- 1L
    while (paste0(delta_name, i) %in% names(data)) i <- i + 1L
    delta_name <- paste0(delta_name, i)
  }
  data2 <- data
  data2[[delta_name]] <- as.integer(!is.na(data2[[outcome_var]]))
  list(data = data2, delta_name = delta_name)
}

#' @keywords internal
el_build_internal_formulas <- function(outcome_sym, resp_expr, aux_expr, delta_name, env) {
  outcome_fml <- as.formula(call("~", outcome_sym, 1L))
  environment(outcome_fml) <- env
  rhs_resp <- if (is.null(resp_expr)) outcome_sym else call("+", outcome_sym, resp_expr)
  response_fml <- as.formula(call("~", as.name(delta_name), rhs_resp))
  environment(response_fml) <- env
  auxiliary_fml <- NULL
  if (!is.null(aux_expr) && length(all.vars(aux_expr)) > 0) {
    rhs <- call("+", 0, aux_expr)
    auxiliary_fml <- as.formula(call("~", rhs))
    environment(auxiliary_fml) <- env
  }
  list(outcome = outcome_fml, response = response_fml, auxiliary = auxiliary_fml)
}

#' Prepare inputs for EL estimation
#' @details Orchestrates parsing and validation via helpers, then returns the
#'   prepared data and internal formulas (outcome, response, auxiliary).
#'   Error messages and behavior are unchanged.
#' @keywords internal
el_prepare_inputs <- function(formula, data, require_na = TRUE, auxiliary_means = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3 || length(all.vars(formula[[2]])) != 1) {
    stop("`formula` must be a two-sided formula with a single variable on the LHS, e.g., y ~ x1 + x2.", call. = FALSE)
  }
  parts <- el_split_rhs(formula)
  env <- parts$env
  outcome_sym <- formula[[2L]]
  outcome_var <- all.vars(outcome_sym)[1]
  aux_vars <- unique(all.vars(parts$aux_expr))
  resp_vars <- if (is.null(parts$resp_expr)) character() else unique(all.vars(parts$resp_expr))

  el_validate_outcome(data, outcome_var, require_na)
  el_validate_vars_present(data, aux_vars, resp_vars)
  respondent_mask <- !is.na(data[[outcome_var]])
  el_validate_aux_no_na(data, parts$aux_expr, respondent_mask = respondent_mask, auxiliary_means = auxiliary_means)

# Response-model predictors must not contain NA among respondents
# Build model.frame for ~ 0 + (outcome + resp_expr)
  resp_rhs <- if (is.null(parts$resp_expr)) outcome_sym else call("+", outcome_sym, parts$resp_expr)
  resp_mf <- tryCatch(
    model.frame(as.formula(call("~", call("+", 0, resp_rhs))), data = data, na.action = stats::na.pass),
    error = function(e) NULL
  )
  if (!is.null(resp_mf) && ncol(resp_mf) > 0) {
    resp_mf_resp <- resp_mf[respondent_mask, , drop = FALSE]
    if (anyNA(resp_mf_resp)) {
      na_by_col <- vapply(seq_len(ncol(resp_mf_resp)), function(j) anyNA(resp_mf_resp[[j]]), logical(1))
      bad_col <- colnames(resp_mf_resp)[which(na_by_col)[1]]
      bad_row_local <- which(is.na(resp_mf_resp[[bad_col]]))[1]
      bad_row <- which(respondent_mask)[bad_row_local]
      msg <- sprintf("Response-model predictor '%s' contains NA values among respondents.%s",
                     bad_col,
                     if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else "")
      stop(msg, call. = FALSE)
    }
  }

  d <- el_make_delta_column(data, outcome_var)
  forms <- el_build_internal_formulas(outcome_sym, parts$resp_expr, parts$aux_expr, d$delta_name, env)
  list(data = d$data, formula_list = forms)
}
