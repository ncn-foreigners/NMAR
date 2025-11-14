#' EL design construction and Formula workflow
#'
#' Parses the user formula once via `Formula::model.frame()` and materializes
#' each RHS block with the standard `model.matrix()` machinery. Auxiliary
#' designs are built on the full dataset (with the intercept and outcome
#' columns removed) while the missingness design uses only respondent rows
#' and always includes an explicit intercept and outcome column.
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
    outcome_var = parsed$outcome_label,
    outcome_source = parsed$outcome_source
  )

  structure(design, class = "el_design_spec")
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
  lhs_vars <- unique(all.vars(lhs_part))
  if (length(lhs_vars) == 0) {
    stop("Left-hand side must reference a variable in `data`.", call. = FALSE)
  }
  if (length(lhs_vars) > 1) {
    stop(
      "Left-hand side may reference only one outcome variable; pre-compute other transforms before modeling.",
      call. = FALSE
    )
  }

  outcome_source <- lhs_vars[[1L]]
  el_validate_outcome(data, outcome_source, require_na = require_na)
  outcome_label <- el_outcome_label(lhs_part)

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

  mask <- !is.na(data[[outcome_source]])
  respondents_only <- all(mask)
  if (isTRUE(require_na) && respondents_only) {
    stop(sprintf("Outcome variable '%s' must contain NA values to indicate nonresponse.", outcome_source), call. = FALSE)
  }

  response_na <- is.na(response_vector)
  transformed_na <- which(mask & response_na)
  if (length(transformed_na) > 0) {
    stop(
      sprintf(
        "LHS expression '%s' produced NA/NaN for observed outcome rows. Ensure the transform is defined for all respondents.",
        outcome_label
      ),
      call. = FALSE
    )
  }

  list(
    fml = fml,
    model_frame = model_frame,
    response_vector = response_vector,
    outcome_label = outcome_label,
    outcome_source = outcome_source,
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

  aux_expr_vars <- if (!is.null(aux_expr)) all.vars(aux_expr) else character(0)
  if (length(aux_expr_vars) > 0 && parsed$outcome_source %in% aux_expr_vars) {
    stop("The outcome cannot appear in auxiliary constraints.", call. = FALSE)
  }

  aux_matrix <- rhs$matrix
  intercept_requested <- isTRUE(attr(rhs$terms, "intercept") == 1) && el_has_explicit_intercept(aux_expr)
  if (ncol(aux_matrix) > 0 && "(Intercept)" %in% colnames(aux_matrix)) {
    aux_matrix <- aux_matrix[, colnames(aux_matrix) != "(Intercept)", drop = FALSE]
  }
  if (intercept_requested) {
    warning("Auxiliary intercepts are ignored; + 1 has no effect.", call. = FALSE)
  }

  drop_names <- c(parsed$outcome_source, parsed$outcome_label)
  if (ncol(aux_matrix) > 0) {
    keep <- !colnames(aux_matrix) %in% drop_names
    aux_matrix <- aux_matrix[, keep, drop = FALSE]
  }

  el_validate_matrix_block(
    mat = aux_matrix,
    mask = parsed$mask,
    row_map = parsed$respondent_indices,
    label = "Auxiliary covariate",
    severity = "error"
  )

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
  colnames(outcome_col) <- parsed$outcome_label

  rhs_predictors <- if (parsed$n_rhs_parts >= 2L) {
    rhs <- el_materialize_rhs(parsed, part = 2L, label = "response predictors")
    if (isTRUE(attr(rhs$terms, "intercept") == 0)) {
      warning("Missingness-model intercept is required; '-1' or '+0' has no effect.", call. = FALSE)
    }
    rhs_matrix <- rhs$matrix
    if ("(Intercept)" %in% colnames(rhs_matrix)) {
      rhs_matrix <- rhs_matrix[, colnames(rhs_matrix) != "(Intercept)", drop = FALSE]
    }
    el_validate_matrix_block(
      mat = rhs_matrix,
      mask = parsed$mask,
      row_map = parsed$respondent_indices,
      label = "Missingness-model predictor",
      severity = "warn"
    )
    rhs_matrix[parsed$mask, , drop = FALSE]
  } else {
    matrix(nrow = n_resp, ncol = 0)
  }

  cbind(intercept_col, outcome_col, rhs_predictors)
}

el_materialize_rhs <- function(parsed, part, label) {
  rhs_formula <- stats::formula(parsed$fml, lhs = 0, rhs = part)
  rhs_expr <- tryCatch(rhs_formula[[2L]], error = function(e) NULL)
  rhs_terms <- el_with_formula_errors(stats::terms(rhs_formula, data = parsed$model_frame), "data")
  el_assert_no_offset(rhs_terms, "data", label)
  mm <- el_with_formula_errors(
    stats::model.matrix(parsed$fml, data = parsed$model_frame, rhs = part),
    "data"
  )
  list(matrix = mm, expr = rhs_expr, terms = rhs_terms)
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


el_validate_matrix_block <- function(mat, mask, row_map, label, severity = c("error", "warn")) {
  if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0 || nrow(mat) == 0) return(invisible(NULL))
  mat_sub <- if (is.null(mask)) mat else mat[mask, , drop = FALSE]
  if (nrow(mat_sub) == 0 || ncol(mat_sub) == 0) return(invisible(NULL))
  el_assert_no_na(mat_sub, row_map, label)
  el_check_constant_columns(mat_sub, label = label, severity = severity)
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

el_outcome_label <- function(expr) {
  paste(deparse(expr), collapse = " ")
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
