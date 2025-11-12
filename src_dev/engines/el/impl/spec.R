#' Construct EL design objects via Formula workflows
#'
#' Builds a single `model.frame()` via `Formula::Formula()` (wrapped in `as.list`
#' to side-step upstream warnings), extracts the outcome, respondent mask, response
#' design matrix (intercept + outcome + RHS2 predictors), and auxiliary matrices
#' (full sample and respondent slice) while enforcing EL-specific constraints.
#' The result is an `el_design` list consumed uniformly by both IID and survey
#' entry points.
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

  formula_obj <- el_as_formula_object(base_formula)
  parts <- length(formula_obj)
  rhs_parts <- if (length(parts) >= 2L) parts[2L] else 1L
  if (rhs_parts > 2L) {
    stop("EL formulas support at most two RHS partitions (auxiliaries | response predictors).", call. = FALSE)
  }

  model_frame <- tryCatch(
    stats::model.frame(formula_obj, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE),
    error = el_rethrow_data_error
  )

  response_vector <- tryCatch(
    Formula::model.part(formula_obj, data = model_frame, lhs = 1, drop = TRUE),
    error = el_rethrow_data_error
  )
  if (!is.numeric(response_vector)) {
    stop("Outcome variable must evaluate to a numeric vector.", call. = FALSE)
  }
  if (length(response_vector) != nrow(model_frame)) {
    stop("Internal error: response extraction did not align with data.", call. = FALSE)
  }

  respondent_mask <- !is.na(response_vector)
  respondent_indices <- which(respondent_mask)

  response_intercept_requested <- el_response_intercept_requested(formula_obj, model_frame, rhs_parts)
  response_design <- el_build_response_design(
    response_vector = response_vector,
    respondent_mask = respondent_mask,
    formula_obj = formula_obj,
    model_frame = model_frame,
    rhs_parts = rhs_parts,
    outcome_var = outcome_var,
    intercept_requested = response_intercept_requested
  )
  el_validate_response_matrix(response_design, respondent_indices)

  if (rhs_parts >= 2L && !response_intercept_requested) {
    warning("Response-model intercept is required; '-1' or '+0' has no effect.", call. = FALSE)
  }

  aux_design <- el_build_auxiliary_design(
    formula_obj = formula_obj,
    model_frame = model_frame,
    respondent_mask = respondent_mask,
    rhs_parts = rhs_parts,
    outcome_var = outcome_var
  )
  el_validate_auxiliary_resp(aux_design$respondent, respondent_indices)

  if (aux_design$intercept_requested) {
    warning("Auxiliary intercepts are ignored; + 1 has no effect.", call. = FALSE)
  }

  structure(
    list(
      response_design = response_design,
      response_outcome = response_vector[respondent_mask],
      aux_resp = aux_design$respondent,
      aux_full = aux_design$full,
      mask = respondent_mask,
      outcome = outcome_var
    ),
    class = "el_design"
  )
}

el_as_formula_object <- function(base_formula) {
  tryCatch(
    withCallingHandlers(
      Formula::Formula(as.list(base_formula)),
      warning = function(w) {
        msg <- conditionMessage(w)
        pattern <- "is.name\\(callee\\).*length\\(x\\) = 2 > 1"
        if (grepl(pattern, msg)) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(err) {
      stop(conditionMessage(err), call. = FALSE)
    }
  )
}

el_build_response_design <- function(response_vector,
                                     respondent_mask,
                                     formula_obj,
                                     model_frame,
                                     rhs_parts,
                                     outcome_var,
                                     intercept_requested) {
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

  out <- cbind(intercept_col, outcome_col, rhs_terms)
  attr(out, "intercept_requested") <- intercept_requested
  out
}

el_response_intercept_requested <- function(formula_obj, model_frame, rhs_parts) {
  if (rhs_parts < 2L) return(TRUE)
  rhs_formula <- stats::formula(formula_obj, lhs = 0, rhs = 2L)
  terms_rhs <- tryCatch(
    stats::terms(rhs_formula, data = model_frame),
    error = function(e) NULL
  )
  if (is.null(terms_rhs)) return(TRUE)
  attr(terms_rhs, "intercept") == 1
}

el_build_auxiliary_design <- function(formula_obj,
                                      model_frame,
                                      respondent_mask,
                                      rhs_parts,
                                      outcome_var) {
  n_total <- nrow(model_frame)
  if (n_total == 0 || rhs_parts < 1L) {
    empty <- matrix(nrow = n_total, ncol = 0)
    return(list(full = empty, respondent = empty[respondent_mask, , drop = FALSE], intercept_requested = FALSE))
  }

  aux_formula <- stats::formula(formula_obj, lhs = 0, rhs = 1L)
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
    stats::model.matrix(formula_obj, data = model_frame, rhs = 1L),
    error = el_rethrow_data_error
  )
  aux_full <- el_drop_intercept_col(aux_full)

  if (ncol(aux_full) > 0 && outcome_var %in% colnames(aux_full)) {
    aux_full <- aux_full[, colnames(aux_full) != outcome_var, drop = FALSE]
  }

  aux_resp <- aux_full[respondent_mask, , drop = FALSE]
  list(full = aux_full, respondent = aux_resp, intercept_requested = intercept_requested)
}
