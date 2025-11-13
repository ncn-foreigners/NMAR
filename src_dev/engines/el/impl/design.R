#' EL design construction and Formula workflow
#'
#' Builds the shared design matrices for the empirical likelihood estimator.
#' A single `model.frame()` is constructed from the supplied `Formula`, the
#' respondent mask is derived from the numeric outcome, and then each RHS part
#' is materialized via the standard Formula `model.matrix()` helpers. Auxiliary
#' matrices are full-data, L-1 coded, and never contain an intercept; missingness
#' matrices are respondent-only with an explicit intercept and outcome column.
#'
#' @keywords internal
el_prepare_design <- function(formula, data, require_na = TRUE) {
  parsed <- el_parse_formula(formula, data, require_na)

  aux_design <- el_build_aux_design(parsed)
  missingness_design <- el_build_missingness_design(parsed)

  design <- list(
    missingness_design = missingness_design,
    auxiliary_design_full = aux_design,
    respondent_mask = parsed$mask,
    outcome_var = parsed$outcome_var,
    response_vector = parsed$response_vector
  )

  structure(design, class = "el_design")
}

el_parse_formula <- function(formula, data, require_na) {
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
  if (!is.symbol(lhs_part)) {
    stop("The left-hand side must be a variable name. Create a column in `data` for transformed outcomes.", call. = FALSE)
  }

  outcome_var <- as.character(lhs_part)
  el_validate_outcome(data, outcome_var, require_na = require_na)

  fml <- el_as_formula(base_formula)
  parts <- length(fml)
  n_rhs_parts <- if (length(parts) >= 2L) parts[2L] else 1L
  if (n_rhs_parts > 2L) {
    stop("EL formulas support at most two RHS partitions (auxiliaries | response predictors).", call. = FALSE)
  }

  model_frame <- el_with_formula_errors(
    stats::model.frame(fml, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE),
    "data"
  )

  response_vector <- el_with_formula_errors(
    Formula::model.part(fml, data = model_frame, lhs = 1, drop = TRUE),
    "data"
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

  list(
    fml = fml,
    model_frame = model_frame,
    response_vector = response_vector,
    outcome_var = outcome_var,
    mask = mask,
    respondent_indices = which(mask),
    n_rhs_parts = n_rhs_parts
  )
}

el_build_aux_design <- function(parsed) {
  n <- nrow(parsed$model_frame)
  if (n == 0 || parsed$n_rhs_parts < 1L) {
    return(matrix(nrow = n, ncol = 0))
  }

  rhs <- el_materialize_rhs(parsed, part = 1L, label = "auxiliary predictors")
  aux_expr <- rhs$expr
  aux_terms <- rhs$terms

  aux_expr_vars <- if (!is.null(aux_expr)) all.vars(aux_expr) else character(0)
  if (length(aux_expr_vars) > 0 && parsed$outcome_var %in% aux_expr_vars) {
    stop("The outcome cannot appear in auxiliary constraints.", call. = FALSE)
  }

  aux_matrix <- drop_intercept(rhs$matrix)
  if (ncol(aux_matrix) > 0 && parsed$outcome_var %in% colnames(aux_matrix)) {
    aux_matrix <- aux_matrix[, colnames(aux_matrix) != parsed$outcome_var, drop = FALSE]
  }

  intercept_requested <- isTRUE(attr(aux_terms, "intercept") == 1) && el_has_explicit_intercept(aux_expr)
  if (intercept_requested) {
    warning("Auxiliary intercepts are ignored; + 1 has no effect.", call. = FALSE)
  }

  if (ncol(aux_matrix) > 0) {
    aux_resp <- aux_matrix[parsed$mask, , drop = FALSE]
    el_assert_no_na(aux_resp, parsed$respondent_indices, "Auxiliary covariate")
    el_check_constant_columns(aux_resp, label = "Auxiliary covariate", severity = "error")
  }

  aux_matrix
}

el_build_missingness_design <- function(parsed) {
  n_resp <- sum(parsed$mask)
  if (n_resp == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  intercept_col <- matrix(1, nrow = n_resp, ncol = 1)
  colnames(intercept_col) <- "(Intercept)"

  outcome_col <- matrix(parsed$response_vector[parsed$mask], ncol = 1)
  colnames(outcome_col) <- parsed$outcome_var

  rhs_predictors <- if (parsed$n_rhs_parts >= 2L) {
    rhs <- el_materialize_rhs(parsed, part = 2L, label = "response predictors")
    if (isTRUE(attr(rhs$terms, "intercept") == 0)) {
      warning("Missingness-model intercept is required; '-1' or '+0' has no effect.", call. = FALSE)
    }
    rhs_matrix <- drop_intercept(rhs$matrix)
    rhs_sub <- rhs_matrix[parsed$mask, , drop = FALSE]
    el_assert_no_na(rhs_sub, parsed$respondent_indices, "Missingness-model predictor")
    el_check_constant_columns(rhs_sub, label = "Missingness-model predictor", severity = "warn")
    rhs_sub
  } else {
    matrix(nrow = n_resp, ncol = 0)
  }

  cbind(intercept_col, outcome_col, rhs_predictors)
}

el_materialize_rhs <- function(parsed, part, label) {
  rhs_formula <- stats::formula(parsed$fml, lhs = 0, rhs = part)
  rhs_expr <- el_rhs_expression(rhs_formula)
  rhs_terms <- el_terms_no_offset(rhs_formula, parsed$model_frame, "data", label)
  terms_for_mm <- rhs_terms
  attr(terms_for_mm, "intercept") <- 1L
  mm <- el_with_formula_errors(stats::model.matrix(terms_for_mm, data = parsed$model_frame), "data")
  list(matrix = mm, expr = rhs_expr, terms = rhs_terms)
}

el_rhs_expression <- function(rhs_formula) {
  tryCatch(rhs_formula[[2L]], error = function(e) NULL)
}

el_as_formula <- function(base_formula) {
  tryCatch(
    Formula::as.Formula(base_formula),
    error = function(err) {
      stop(conditionMessage(err), call. = FALSE)
    }
  )
}

el_with_formula_errors <- function(expr, context_label) {
  tryCatch(
    expr,
    error = function(e) el_rethrow_data_error(e, context_label)
  )
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

el_terms_no_offset <- function(rhs_formula, model_frame, context_label, label) {
  tr <- el_with_formula_errors(stats::terms(rhs_formula, data = model_frame), context_label)
  el_assert_no_offset(tr, context_label, label)
  tr
}

drop_intercept <- function(mm) {
  if (is.null(mm) || !is.matrix(mm) || ncol(mm) == 0) return(mm)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  mm
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
  if (!is.numeric(design$response_vector) || length(design$response_vector) != length(mask)) {
    stop("Internal error: response_vector must align with respondent mask.", call. = FALSE)
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
