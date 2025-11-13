#' EL design construction and Formula workflow
#'
#' Builds the shared design matrices for the empirical likelihood engine. The
#' function constructs one `model.frame` from the supplied `Formula`, extracts
#' the respondent indicators, and materializes the auxiliary (RHS1) and
#' missingness (intercept + outcome + RHS2) matrices. Offsets (`offset()`) are
#' not supported in either partition; users should add predictors explicitly.
#'
#' @keywords internal
el_prepare_design <- function(formula, data, require_na = TRUE, context_label = "data") {
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
  el_validate_outcome(data, outcome_var, require_na = require_na)

  fml <- el_as_formula(base_formula)
  parts <- length(fml)
  n_rhs_parts <- if (length(parts) >= 2L) parts[2L] else 1L
  if (n_rhs_parts > 2L) {
    stop("EL formulas support at most two RHS partitions (auxiliaries | response predictors).", call. = FALSE)
  }

  model_frame <- el_eval_formula_call(
    stats::model.frame(fml, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE),
    context_label
  )

  response_vector <- el_eval_formula_call(
    Formula::model.part(fml, data = model_frame, lhs = 1, drop = TRUE),
    context_label
  )
  if (!is.numeric(response_vector)) {
    stop("Outcome variable must evaluate to a numeric vector.", call. = FALSE)
  }
  if (length(response_vector) != nrow(model_frame)) {
    stop("Internal error: response extraction did not align with data.", call. = FALSE)
  }

  mask <- !is.na(response_vector)
  respondents_only <- all(mask)
  if (isTRUE(require_na) && respondents_only) {
    stop(sprintf("Outcome variable '%s' must contain NA values to indicate nonresponse.", outcome_var), call. = FALSE)
  }
  respondent_indices <- which(mask)

  aux_design <- el_build_aux_design(
    fml = fml,
    model_frame = model_frame,
    mask = mask,
    n_rhs_parts = n_rhs_parts,
    outcome_var = outcome_var,
    context_label = context_label,
    respondent_indices = respondent_indices
  )

  response_design <- el_build_missingness_design(
    fml = fml,
    model_frame = model_frame,
    response_vector = response_vector,
    mask = mask,
    n_rhs_parts = n_rhs_parts,
    outcome_var = outcome_var,
    context_label = context_label,
    respondent_indices = respondent_indices
  )

  structure(
    list(
      missingness_design = response_design,
      y_obs = response_vector[mask],
      auxiliary_design_full = aux_design,
      respondent_mask = mask,
      respondents_only = respondents_only,
      outcome_var = outcome_var,
      has_aux = attr(aux_design, "has_aux") %||% (ncol(aux_design) > 0)
    ),
    class = "el_design"
  )
}

el_as_formula <- function(base_formula) {
  tryCatch(
    Formula::as.Formula(base_formula),
    error = function(err) {
      stop(conditionMessage(err), call. = FALSE)
    }
  )
}

el_eval_formula_call <- function(expr, context_label) {
  expr <- substitute(expr)
  eval(
    bquote(tryCatch(.(expr), error = function(e) el_rethrow_data_error(e, context_label))),
    envir = parent.frame()
  )
}

el_build_missingness_design <- function(fml,
                                        model_frame,
                                        response_vector,
                                        mask,
                                        n_rhs_parts,
                                        outcome_var,
                                        context_label,
                                        respondent_indices) {
  n_resp <- sum(mask)
  if (n_resp == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  intercept_col <- matrix(1, nrow = n_resp, ncol = 1)
  colnames(intercept_col) <- "(Intercept)"

  outcome_col <- matrix(response_vector[mask], ncol = 1)
  colnames(outcome_col) <- outcome_var

  rhs_predictors <- if (n_rhs_parts >= 2L) {
    rhs_formula <- stats::formula(fml, lhs = 0, rhs = 2L)
    rhs_terms <- el_eval_formula_call(stats::terms(rhs_formula, data = model_frame), context_label)
    el_assert_no_offset(rhs_terms, context_label, "response predictors")
    if (isTRUE(attr(rhs_terms, "intercept") == 0)) {
      warning("Missingness-model intercept is required; '-1' or '+0' has no effect.", call. = FALSE)
    }
    rhs_matrix <- el_eval_formula_call(stats::model.matrix(fml, data = model_frame, rhs = 2L), context_label)
    if (ncol(rhs_matrix) > 0 && "(Intercept)" %in% colnames(rhs_matrix)) {
      rhs_matrix <- rhs_matrix[, colnames(rhs_matrix) != "(Intercept)", drop = FALSE]
    }
    rhs_sub <- rhs_matrix[mask, , drop = FALSE]
    el_assert_no_na(rhs_sub, respondent_indices, "Missingness-model predictor")
    el_check_constant_columns(rhs_sub, label = "Missingness-model predictor", severity = "warn")
    rhs_sub
  } else {
    matrix(nrow = n_resp, ncol = 0)
  }

  cbind(intercept_col, outcome_col, rhs_predictors)
}

el_build_aux_design <- function(fml,
                                model_frame,
                                mask,
                                n_rhs_parts,
                                outcome_var,
                                context_label,
                                respondent_indices) {
  if (nrow(model_frame) == 0 || n_rhs_parts < 1L) {
    out <- matrix(nrow = nrow(model_frame), ncol = 0)
    attr(out, "has_aux") <- FALSE
    return(out)
  }

  aux_formula <- stats::formula(fml, lhs = 0, rhs = 1L)
  aux_terms <- el_eval_formula_call(stats::terms(aux_formula, data = model_frame), context_label)
  el_assert_no_offset(aux_terms, context_label, "auxiliary predictors")
  aux_expr <- el_rhs_part(fml, part = 1L)
  aux_vars_expr <- if (!is.null(aux_expr)) all.vars(aux_expr) else character(0)
  if (length(aux_vars_expr) > 0 && outcome_var %in% aux_vars_expr) {
    stop("The outcome cannot appear in auxiliary constraints.", call. = FALSE)
  }
  intercept_requested <- isTRUE(attr(aux_terms, "intercept") == 1) && el_has_explicit_intercept(aux_expr)
  if (intercept_requested) {
    warning("Auxiliary intercepts are ignored; + 1 has no effect.", call. = FALSE)
  }

  aux_matrix <- el_eval_formula_call(stats::model.matrix(fml, data = model_frame, rhs = 1L), context_label)
  if (ncol(aux_matrix) > 0 && "(Intercept)" %in% colnames(aux_matrix)) {
    aux_matrix <- aux_matrix[, colnames(aux_matrix) != "(Intercept)", drop = FALSE]
  }
  if (ncol(aux_matrix) > 0 && outcome_var %in% colnames(aux_matrix)) {
    aux_matrix <- aux_matrix[, colnames(aux_matrix) != outcome_var, drop = FALSE]
  }

  if (ncol(aux_matrix) > 0) {
    aux_resp <- aux_matrix[mask, , drop = FALSE]
    el_assert_no_na(aux_resp, respondent_indices, "Auxiliary covariate")
    el_check_constant_columns(aux_resp, label = "Auxiliary covariate", severity = "error")
  }

  attr(aux_matrix, "has_aux") <- ncol(aux_matrix) > 0
  aux_matrix
}

el_rhs_part <- function(fml, part = 1L) {
  rhs_parts <- attr(fml, "rhs")
  if (is.null(rhs_parts) || length(rhs_parts) < part) return(NULL)
  rhs_parts[[part]]
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

el_check_constant_columns <- function(mat, label, severity = c("error", "warn")) {
  severity <- match.arg(severity)
  if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0 || nrow(mat) == 0) return(invisible(NULL))

  is_constant <- vapply(seq_len(ncol(mat)), function(j) {
    col <- mat[, j]
    if (all(is.na(col))) return(FALSE)
    all(col == col[which(!is.na(col))[1]], na.rm = TRUE)
  }, logical(1))

  if (!any(is_constant)) return(invisible(NULL))
  cols <- colnames(mat)[is_constant]
  msg <- sprintf(
    "%s %s zero variance among respondents: %s",
    label,
    if (length(cols) == 1) "has" else "have",
    paste(cols, collapse = ", ")
  )
  if (identical(severity, "warn")) {
    warning(msg, call. = FALSE)
    return(invisible(NULL))
  }
  stop(msg, call. = FALSE)
}

el_assert_no_na <- function(mat, row_map = NULL, label, scope_note = NULL, plural_label = FALSE) {
  if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0 || nrow(mat) == 0) return(invisible(NULL))
  if (!anyNA(mat)) return(invisible(NULL))
  na_loc <- which(is.na(mat), arr.ind = TRUE)[1, , drop = TRUE]
  bad_col <- colnames(mat)[na_loc[2]]
  row_label <- if (!is.null(row_map) && length(row_map) >= na_loc[1]) row_map[na_loc[1]] else na_loc[1]
  note <- scope_note %||% if (!is.null(row_map)) " among respondents" else ""
  msg <- if (isTRUE(plural_label)) {
    sprintf(
      "%s contain NA values%s.%s",
      label,
      note,
      if (isTRUE(is.finite(row_label))) sprintf("\nFirst NA in '%s' at row %d", bad_col, row_label) else ""
    )
  } else {
    sprintf(
      "%s '%s' contains NA values%s.%s",
      label,
      bad_col,
      note,
      if (isTRUE(is.finite(row_label))) sprintf("\nFirst NA at row %d", row_label) else ""
    )
  }
  stop(msg, call. = FALSE)
}

el_rethrow_data_error <- function(err, context_label) {
  msg <- conditionMessage(err)
  stop(sprintf("Formula evaluation failed for %s: %s", context_label, msg), call. = FALSE)
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

el_validate_design_spec <- function(design, data_nrow, context_label) {
  mask <- design$respondent_mask
  if (!is.logical(mask)) {
    stop(sprintf("Internal error: respondent mask must be logical for %s.", context_label), call. = FALSE)
  }
  if (length(mask) != data_nrow) {
    stop(sprintf("Internal error: respondent mask length (%d) must equal %s rows (%d).", length(mask), context_label, data_nrow), call. = FALSE)
  }
  missingness_design <- design$missingness_design
  if (!is.null(missingness_design) && nrow(missingness_design) != sum(mask)) {
    stop("Internal error: missingness design rows must equal respondent count.", call. = FALSE)
  }
  if (!is.null(design$auxiliary_design_full) && nrow(design$auxiliary_design_full) != data_nrow) {
    stop("Internal error: auxiliary_design_full must align with original data rows.", call. = FALSE)
  }
  invisible(design)
}

el_assert_no_offset <- function(terms_obj, context_label, label) {
  if (is.null(terms_obj)) return(invisible(NULL))
  offsets <- attr(terms_obj, "offset")
  if (!is.null(offsets) && length(offsets) > 0) {
    stop(
      sprintf(
        "Offsets (offset()) are not supported in %s for %s. Remove offset() from the formula.",
        context_label,
        label
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}
