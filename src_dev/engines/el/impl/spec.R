#' Prepare EL design matrices via Formula workflows
#'
#' Builds a single `model.frame()` with `Formula::Formula()`, extracts the
#' respondent-only response matrix (intercept, outcome, RHS2 predictors) and the
#' auxiliary design matrices (full sample plus respondent slice), and enforces
#' EL-specific constraints such as outcome exclusion from auxiliaries and
#' automatic intercept removal (explicit `+ 1` yields a warning). The result is
#' an `el_design` list consumed by both IID and survey entry points.
#'
#' @keywords internal
el_prepare_design <- function(formula, data, require_na = TRUE) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  }

  outcome_expr <- formula[[2L]]
  outcome_vars <- all.vars(outcome_expr)
  if (length(outcome_vars) != 1L) {
    stop("The left-hand side must contain exactly one outcome variable.", call. = FALSE)
  }

  outcome_var <- outcome_vars[1]
  el_validate_outcome(data, outcome_var, require_na)

  formula_obj <- Formula::Formula(formula)
  parts <- length(formula_obj)
  rhs_parts <- if (length(parts) >= 2L) parts[2L] else 1L
  if (rhs_parts > 2L) {
    stop("EL formulas support at most two RHS partitions (auxiliaries | response predictors).", call. = FALSE)
  }

  model_frame <- tryCatch(
    stats::model.frame(formula_obj, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE),
    error = el_rethrow_data_error
  )

  response_vector <- tryCatch(stats::model.response(model_frame), error = el_rethrow_data_error)
  if (is.null(response_vector)) {
    stop("Unable to extract a response vector from the supplied formula.", call. = FALSE)
  }
  if (is.matrix(response_vector)) {
    response_vector <- as.numeric(response_vector)
  }
  if (!is.numeric(response_vector)) {
    stop("Outcome variable must evaluate to a numeric vector.", call. = FALSE)
  }

  respondent_mask <- !is.na(response_vector)

  response_matrix <- el_build_response_matrix(
    response_vector = response_vector,
    respondent_mask = respondent_mask,
    formula_obj = formula_obj,
    model_frame = model_frame,
    rhs_parts = rhs_parts,
    outcome_var = outcome_var
  )
  respondent_indices <- which(respondent_mask)
  el_validate_response_matrix(response_matrix, respondent_indices)

  aux_design <- el_build_auxiliary_matrices(
    formula_obj = formula_obj,
    model_frame = model_frame,
    respondent_mask = respondent_mask,
    rhs_parts = rhs_parts,
    outcome_var = outcome_var
  )
  el_validate_auxiliary_resp(aux_design$respondent, respondent_indices)

  structure(
    list(
      response = response_matrix,
      aux_resp = aux_design$respondent,
      aux_full = aux_design$full,
      mask = respondent_mask,
      outcome = outcome_var
    ),
    class = "el_design"
  )
}

el_build_response_matrix <- function(response_vector,
                                     respondent_mask,
                                     formula_obj,
                                     model_frame,
                                     rhs_parts,
                                     outcome_var) {
  n_resp <- sum(respondent_mask)
  if (n_resp == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  intercept_col <- matrix(1, nrow = n_resp, ncol = 1)
  colnames(intercept_col) <- "(Intercept)"

  outcome_col <- matrix(response_vector[respondent_mask], ncol = 1)
  colnames(outcome_col) <- outcome_var

  rhs_terms <- if (rhs_parts >= 2L) {
    mm <- tryCatch(
      stats::model.matrix(formula_obj, data = model_frame, rhs = 2L),
      error = el_rethrow_data_error
    )
    mm_resp <- mm[respondent_mask, , drop = FALSE]
    el_drop_intercept_col(mm_resp)
  } else {
    matrix(nrow = n_resp, ncol = 0)
  }

  cbind(intercept_col, outcome_col, rhs_terms)
}

el_build_auxiliary_matrices <- function(formula_obj,
                                        model_frame,
                                        respondent_mask,
                                        rhs_parts,
                                        outcome_var) {
  if (nrow(model_frame) == 0) {
    empty_resp <- matrix(nrow = 0, ncol = 0)
    empty_full <- matrix(nrow = 0, ncol = 0)
    return(list(full = empty_full, respondent = empty_resp))
  }

  aux_formula <- if (rhs_parts >= 1L) {
    stats::formula(formula_obj, lhs = 0, rhs = 1L)
  } else {
    stats::as.formula("~ 0")
  }

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
    el_formula_has_explicit_intercept(rhs_expr)

  aux_full <- tryCatch(
    stats::model.matrix(aux_terms, data = model_frame),
    error = el_rethrow_data_error
  )
  aux_full <- el_drop_intercept_col(aux_full)

  if (ncol(aux_full) > 0 && outcome_var %in% colnames(aux_full)) {
    aux_full <- aux_full[, colnames(aux_full) != outcome_var, drop = FALSE]
  }

  if (intercept_requested) {
    warning("Auxiliary intercepts are ignored; + 1 has no effect.", call. = FALSE)
  }

  if (ncol(aux_full) == 0) {
    full_zero <- matrix(nrow = nrow(model_frame), ncol = 0)
    resp_zero <- matrix(nrow = sum(respondent_mask), ncol = 0)
    return(list(full = full_zero, respondent = resp_zero))
  }

  aux_resp <- aux_full[respondent_mask, , drop = FALSE]
  list(full = aux_full, respondent = aux_resp)
}

el_drop_intercept_col <- function(mm) {
  if (is.null(mm) || !is.matrix(mm) || ncol(mm) == 0) return(mm)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  mm
}
