#' EL design construction and Formula workflow
#'
#' Centralizes Formula parsing and design-matrix construction for the EL engine.
#' Builds a single model.frame from a Formula object and materializes:
#'  - missingness_model_matrix: (Intercept) + outcome + RHS2 predictors (RHS2 intercept dropped)
#'  - y_obs: outcome values on respondent rows (non-missing)
#'  - aux_mm_full: RHS1 auxiliary design on full data (no intercept; outcome dropped)
#'  - respondent_mask: logical mask for respondents (non-missing outcomes)
#'  - outcome_var: outcome column name
#'
#' Also provides narrowly-scoped helpers used only by this workflow.
#'
#' @keywords internal
el_construct_design <- function(formula, data, require_na = TRUE) {
  if (missing(formula)) {
    stop("`formula` must be supplied.", call. = FALSE)
  }

  base_formula <- stats::as.formula(formula)
  if (!inherits(base_formula, "formula") || length(base_formula) != 3L) {
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

# Detect intercept presence on RHS2 (missingness-model predictors) via terms();
# the missingness-model intercept is enforced even if suppressed in the formula.
  rhs2_has_intercept <- if (n_rhs_parts >= 2L) {
    rhs_formula <- stats::formula(fml, lhs = 0, rhs = 2L)
    terms_rhs <- tryCatch(stats::terms(rhs_formula, data = model_frame), error = function(e) NULL)
    if (is.null(terms_rhs)) TRUE else isTRUE(attr(terms_rhs, "intercept") == 1)
  } else TRUE

  missingness_model_matrix <- el_build_missingness_matrix(
    response_vector = response_vector,
    mask = mask,
    fml = fml,
    model_frame = model_frame,
    n_rhs_parts = n_rhs_parts,
    outcome_var = outcome_var
  )
  el_validate_missingness_predictors(missingness_model_matrix, respondent_indices)

  if (n_rhs_parts >= 2L && !rhs2_has_intercept) {
    warning("Missingness-model intercept is required; '-1' or '+0' has no effect.", call. = FALSE)
  }

  aux_design <- el_build_aux_matrix_full(
    fml = fml,
    model_frame = model_frame,
    mask = mask,
    n_rhs_parts = n_rhs_parts,
    outcome_var = outcome_var
  )
# Early NA check among respondents for auxiliaries (if any)
  el_validate_aux_respondent_matrix(aux_design$full[mask, , drop = FALSE], respondent_indices)

  if (aux_design$intercept_requested) {
    warning("Auxiliary intercepts are ignored; + 1 has no effect.", call. = FALSE)
  }

  structure(
    list(
      missingness_model_matrix = missingness_model_matrix,
      y_obs = response_vector[mask],
      aux_mm_full = aux_design$full,
      respondent_mask = mask,
      outcome_var = outcome_var
    ),
    class = "el_design"
  )
}

el_as_formula <- function(base_formula) {
# Prefer as.Formula() for idiomatic conversion; keep as.list() wrapper
# to sidestep upstream parsing edge cases when passing plain formula.
  tryCatch(
    Formula::as.Formula(as.list(base_formula)),
    error = function(err) {
      stop(conditionMessage(err), call. = FALSE)
    }
  )
}

el_build_missingness_matrix <- function(response_vector,
                                        mask,
                                        fml,
                                        model_frame,
                                        n_rhs_parts,
                                        outcome_var) {
  n_resp <- sum(mask)
  if (n_resp == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  intercept_col <- matrix(1, nrow = n_resp, ncol = 1)
  colnames(intercept_col) <- "(Intercept)"

  outcome_col <- matrix(response_vector[mask], ncol = 1)
  colnames(outcome_col) <- outcome_var

  rhs_terms <- if (n_rhs_parts >= 2L) {
    mm <- tryCatch(
      stats::model.matrix(fml, data = model_frame, rhs = 2L),
      error = el_rethrow_data_error
    )
    mm_resp <- mm[mask, , drop = FALSE]
    el_drop_intercept(mm_resp)
  } else {
    matrix(nrow = n_resp, ncol = 0)
  }

  cbind(intercept_col, outcome_col, rhs_terms)
}

el_build_aux_matrix_full <- function(fml,
                                     model_frame,
                                     mask,
                                     n_rhs_parts,
                                     outcome_var) {
  n_total <- nrow(model_frame)
  if (n_total == 0 || n_rhs_parts < 1L) {
    empty <- matrix(nrow = n_total, ncol = 0)
    return(list(full = empty, intercept_requested = FALSE))
  }

  aux_formula <- stats::formula(fml, lhs = 0, rhs = 1L)
  aux_vars <- all.vars(aux_formula)
  if (outcome_var %in% aux_vars) {
    stop("The outcome cannot appear in auxiliary constraints.", call. = FALSE)
  }

  aux_terms <- tryCatch(
    stats::terms(aux_formula, data = model_frame),
    error = el_rethrow_data_error
  )
  rhs_expr <- if (length(aux_formula) >= 3L) aux_formula[[3L]] else aux_formula[[2L]]
  intercept_requested <- attr(aux_terms, "intercept") == 1 &&
    el_has_explicit_intercept(rhs_expr)

  aux_full <- tryCatch(
    stats::model.matrix(fml, data = model_frame, rhs = 1L),
    error = el_rethrow_data_error
  )
  aux_full <- el_drop_intercept(aux_full)

  if (ncol(aux_full) > 0 && outcome_var %in% colnames(aux_full)) {
    aux_full <- aux_full[, colnames(aux_full) != outcome_var, drop = FALSE]
  }

  list(full = aux_full, intercept_requested = intercept_requested)
}

# Helpers used only in the EL Formula/design workflow
el_drop_intercept <- function(mm) {
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

el_validate_missingness_predictors <- function(response_matrix, respondent_indices) {
  if (!is.matrix(response_matrix) || ncol(response_matrix) == 0 || nrow(response_matrix) == 0) {
    return(invisible(NULL))
  }
  if (!anyNA(response_matrix)) return(invisible(NULL))

  na_loc <- which(is.na(response_matrix), arr.ind = TRUE)[1, , drop = TRUE]
  bad_col <- colnames(response_matrix)[na_loc[2]]
  bad_row <- if (length(respondent_indices) >= na_loc[1]) respondent_indices[na_loc[1]] else NA_integer_

  msg <- sprintf(
    "Missingness-model predictor '%s' contains NA values among respondents.%s",
    bad_col,
    if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else ""
  )
  stop(msg, call. = FALSE)
}

el_validate_aux_respondent_matrix <- function(aux_resp, respondent_indices) {
  if (is.null(aux_resp) || !is.matrix(aux_resp) || ncol(aux_resp) == 0) return(invisible(NULL))

  if (anyNA(aux_resp)) {
    na_loc <- which(is.na(aux_resp), arr.ind = TRUE)[1, , drop = TRUE]
    bad_col <- colnames(aux_resp)[na_loc[2]]
    bad_row <- if (length(respondent_indices) >= na_loc[1]) respondent_indices[na_loc[1]] else NA_integer_
    msg <- sprintf(
      "Covariate '%s' contains NA values among respondents.%s",
      bad_col,
      if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else ""
    )
    stop(msg, call. = FALSE)
  }
  invisible(NULL)
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
