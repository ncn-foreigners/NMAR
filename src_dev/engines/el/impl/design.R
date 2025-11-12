#' EL design construction and Formula workflow
#'
#' Centralizes Formula parsing and design-matrix construction for the EL engine.
#' Builds a single model.frame from a Formula object and materializes:
#'  - response_matrix: (Intercept) + outcome + RHS2 predictors (RHS2 intercept dropped)
#'  - y_obs: outcome values on respondent rows (non-missing)
#'  - aux_full: RHS1 auxiliary design on full data (no intercept; outcome dropped)
#'  - respondent_mask: logical mask for respondents (non-missing outcomes)
#'  - outcome_var: outcome column name
#'
#' @keywords internal
el_parse_design <- function(formula, data, require_na = TRUE) {
  if (missing(formula)) stop("`formula` must be supplied.", call. = FALSE)

  base_formula <- stats::as.formula(formula)
  if (!inherits(base_formula, "formula")) {
    stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  }
  lhs_part <- tryCatch(base_formula[[2L]], error = function(e) NULL)
  rhs_part <- tryCatch(base_formula[[3L]], error = function(e) NULL)
  if (is.null(lhs_part) || is.null(rhs_part)) {
    stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  }

  outcome_expr <- base_formula[[2L]]
  outcome_vars <- all.vars(outcome_expr)
  if (length(outcome_vars) != 1L) {
    stop("The left-hand side must contain exactly one outcome variable.", call. = FALSE)
  }
  outcome_var <- outcome_vars[1L]
  el_validate_outcome(data, outcome_var, require_na)

  fml <- el_as_formula(base_formula)
  parts <- length(fml)
  n_rhs_parts <- if (length(parts) >= 2L) parts[2L] else 1L
  if (n_rhs_parts > 2L) {
    stop("EL formulas support at most two RHS partitions (auxiliaries | response predictors).", call. = FALSE)
  }

  model_frame <- tryCatch(
    stats::model.frame(fml, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE),
    error = el_rethrow_data_error
  )

  response_vector <- tryCatch(
    Formula::model.part(fml, data = model_frame, lhs = 1, drop = TRUE),
    error = el_rethrow_data_error
  )
  if (!is.numeric(response_vector)) {
    stop("Outcome variable must evaluate to a numeric vector.", call. = FALSE)
  }
  if (length(response_vector) != nrow(model_frame)) {
    stop("Internal error: response extraction did not align with data.", call. = FALSE)
  }

  mask <- !is.na(response_vector)
  respondent_indices <- which(mask)

  aux_full <- el_build_aux_matrix(fml, model_frame, n_rhs_parts, outcome_var)
  if (ncol(aux_full) > 0) {
    el_abort_if_na(aux_full[mask, , drop = FALSE], respondent_indices, "Auxiliary covariate")
  }

  response_matrix <- el_build_response_matrix(
    fml = fml,
    model_frame = model_frame,
    response_vector = response_vector,
    mask = mask,
    n_rhs_parts = n_rhs_parts,
    outcome_var = outcome_var
  )
  el_abort_if_na(response_matrix, respondent_indices, "Missingness-model predictor")

  structure(
    list(
      response_matrix = response_matrix,
      y_obs = response_vector[mask],
      aux_full = aux_full,
      respondent_mask = mask,
      outcome_var = outcome_var
    ),
    class = "el_design"
  )
}

el_as_formula <- function(base_formula) {
# Wrap in as.list() before as.Formula() to avoid edge cases when the caller
# supplies language objects instead of plain formula literals.
  tryCatch(
    Formula::as.Formula(as.list(base_formula)),
    error = function(err) {
      stop(conditionMessage(err), call. = FALSE)
    }
  )
}

el_build_response_matrix <- function(fml, model_frame, response_vector, mask, n_rhs_parts, outcome_var) {
  n_resp <- sum(mask)
  if (n_resp == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  intercept_col <- matrix(1, nrow = n_resp, ncol = 1)
  colnames(intercept_col) <- "(Intercept)"
  outcome_col <- matrix(response_vector[mask], ncol = 1)
  colnames(outcome_col) <- outcome_var

  rhs_terms <- if (n_rhs_parts >= 2L) {
    rhs_formula <- stats::formula(fml, lhs = 0, rhs = 2L)
    rhs_terms_obj <- tryCatch(stats::terms(rhs_formula, data = model_frame), error = el_rethrow_data_error)
    if (isTRUE(attr(rhs_terms_obj, "intercept") == 0)) {
      warning("Missingness-model intercept is required; '-1' or '+0' has no effect.", call. = FALSE)
    }
    rhs_mm <- tryCatch(stats::model.matrix(fml, data = model_frame, rhs = 2L), error = el_rethrow_data_error)
    el_remove_intercept(rhs_mm)[mask, , drop = FALSE]
  } else {
    matrix(nrow = n_resp, ncol = 0)
  }

  cbind(intercept_col, outcome_col, rhs_terms)
}

el_build_aux_matrix <- function(fml, model_frame, n_rhs_parts, outcome_var) {
  if (nrow(model_frame) == 0 || n_rhs_parts < 1L) {
    return(matrix(nrow = nrow(model_frame), ncol = 0))
  }

  aux_formula <- stats::formula(fml, lhs = 0, rhs = 1L)
  aux_vars <- all.vars(aux_formula)
  if (outcome_var %in% aux_vars) {
    stop("The outcome cannot appear in auxiliary constraints.", call. = FALSE)
  }

  aux_terms <- tryCatch(stats::terms(aux_formula, data = model_frame), error = el_rethrow_data_error)
  rhs_expr <- if (length(aux_formula) >= 3L) aux_formula[[3L]] else aux_formula[[2L]]
  intercept_requested <- isTRUE(attr(aux_terms, "intercept") == 1) && el_has_explicit_intercept(rhs_expr)
  if (intercept_requested) {
    warning("Auxiliary intercepts are ignored; + 1 has no effect.", call. = FALSE)
  }

  aux_full <- tryCatch(stats::model.matrix(fml, data = model_frame, rhs = 1L), error = el_rethrow_data_error)
  aux_full <- el_remove_intercept(aux_full)
  if (ncol(aux_full) > 0 && outcome_var %in% colnames(aux_full)) {
    aux_full <- aux_full[, colnames(aux_full) != outcome_var, drop = FALSE]
  }
  aux_full
}

el_remove_intercept <- function(mm) {
  if (is.null(mm) || !is.matrix(mm) || ncol(mm) == 0) return(mm)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  mm
}

el_has_explicit_intercept <- function(node) {
  if (is.null(node)) return(FALSE)
  if (is.numeric(node) && length(node) == 1L && isTRUE(all.equal(node, 1))) return(TRUE)
  if (is.call(node)) {
    op <- as.character(node[[1L]])
    if (op %in% c("+", "-", "~")) {
      for (i in 2L:length(node)) {
        if (Recall(node[[i]])) return(TRUE)
      }
    }
  }
  FALSE
}

el_abort_if_na <- function(mat, respondent_indices, label) {
  if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0 || nrow(mat) == 0) return(invisible(NULL))
  if (!anyNA(mat)) return(invisible(NULL))
  na_loc <- which(is.na(mat), arr.ind = TRUE)[1, , drop = TRUE]
  bad_col <- colnames(mat)[na_loc[2]]
  bad_row <- if (length(respondent_indices) >= na_loc[1]) respondent_indices[na_loc[1]] else NA_integer_
  msg <- sprintf(
    "%s '%s' contains NA values among respondents.%s",
    label,
    bad_col,
    if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else ""
  )
  stop(msg, call. = FALSE)
}

el_rethrow_data_error <- function(err) {
  msg <- conditionMessage(err)
  missing_pattern <- "object '([^']+)' not found"
  if (grepl(missing_pattern, msg, perl = TRUE)) {
    missing_var <- sub(missing_pattern, "\\1", msg, perl = TRUE)
    stop(sprintf("Variables not found in data: %s", missing_var), call. = FALSE)
  }
  stop(msg, call. = FALSE)
}

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
