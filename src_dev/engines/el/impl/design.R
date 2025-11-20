#' EL design construction and Formula workflow
#'
#' Builds a single model frame and materializes RHS blocks with standard
#' `model.matrix()` calls so IID and survey entry points share the same parsing
#' logic. Auxiliary designs are built on the full dataset (intercept removed)
#' while the missingness design uses only respondent rows and always includes an
#' explicit intercept and outcome column.
#'
#' @keywords internal
el_process_design <- function(formula, data) {
  if (missing(formula)) stop("`formula` must be supplied.", call. = FALSE)

  base_formula <- stats::as.formula(formula)
  if (!inherits(base_formula, "formula")) {
    stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  }

  lhs_expr <- tryCatch(base_formula[[2L]], error = function(e) NULL)
  rhs_expr <- tryCatch(base_formula[[3L]], error = function(e) NULL)
  if (is.null(lhs_expr) || is.null(rhs_expr)) {
    stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  }

  lhs_vars <- unique(all.vars(lhs_expr))
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
  if (!outcome_source %in% names(data)) {
    stop(sprintf("Variables not found in data: %s", outcome_source), call. = FALSE)
  }
  outcome_label <- el_outcome_label(lhs_expr)

  fml <- tryCatch(
    Formula::Formula(base_formula),
    error = function(err) {
      stop(conditionMessage(err), call. = FALSE)
    }
  )
  parts <- length(fml)
  n_rhs_parts <- if (length(parts) >= 2L) parts[2L] else 1L
  if (n_rhs_parts > 2L) {
    stop("EL formulas support at most two RHS partitions (auxiliaries | response predictors).", call. = FALSE)
  }

  model_frame <- tryCatch(
    stats::model.frame(fml, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE),
    error = function(e) {
      stop(sprintf("Formula evaluation failed for data: %s", conditionMessage(e)), call. = FALSE)
    }
  )

  response_vector <- tryCatch(
    Formula::model.part(fml, data = model_frame, lhs = 1, drop = TRUE),
    error = function(e) {
      stop(sprintf("Formula evaluation failed for data: %s", conditionMessage(e)), call. = FALSE)
    }
  )
  if (!is.numeric(response_vector)) {
    stop("Outcome variable must be numeric after evaluating the left-hand side.", call. = FALSE)
  }
  if (length(response_vector) != nrow(model_frame)) {
    stop("Internal error: response extraction did not align with data.", call. = FALSE)
  }

  raw_mask <- !is.na(data[[outcome_source]])
  transformed_na <- which(raw_mask & is.na(response_vector))
  if (length(transformed_na) > 0) {
    stop(
      sprintf(
        "LHS expression '%s' produced NA/NaN for observed outcome rows. Ensure the transform is defined for all respondents.",
        outcome_label
      ),
      call. = FALSE
    )
  }

  mask <- unname(!is.na(response_vector))
  if (!any(mask)) {
    stop("No respondents detected in data after preprocessing.", call. = FALSE)
  }
  respondent_indices <- which(mask)

  aux_design_full <- el_build_auxiliary_matrix(
    fml = fml,
    model_frame = model_frame,
    n_rhs_parts = n_rhs_parts,
    outcome_source = outcome_source,
    outcome_label = outcome_label,
    mask = mask
  )

  missingness_design <- el_build_missingness_matrix(
    fml = fml,
    model_frame = model_frame,
    n_rhs_parts = n_rhs_parts,
    outcome_label = outcome_label,
    mask = mask,
    respondent_indices = respondent_indices
  )

  design <- list(
    missingness_design = missingness_design,
    auxiliary_design_full = aux_design_full,
    respondent_mask = mask,
    outcome_var = outcome_label,
    outcome_source = outcome_source,
    strata_factor = data[["..nmar_strata_factor.."]] %||% NULL
  )

  structure(design, class = "el_design_spec")
}

el_build_auxiliary_matrix <- function(fml,
                                      model_frame,
                                      n_rhs_parts,
                                      outcome_source,
                                      outcome_label,
                                      mask) {
  n <- nrow(model_frame)
  if (n == 0 || n_rhs_parts < 1L) {
    return(matrix(nrow = n, ncol = 0))
  }

  rhs_formula <- stats::formula(fml, lhs = 0, rhs = 1)
  rhs_expr <- tryCatch(rhs_formula[[2L]], error = function(e) NULL)
  rhs_terms <- tryCatch(
    stats::terms(rhs_formula, data = model_frame),
    error = function(e) stop(sprintf("Formula evaluation failed for data: %s", conditionMessage(e)), call. = FALSE)
  )
  el_assert_no_offset(rhs_terms, "data", "auxiliary predictors")
  aux_matrix <- tryCatch(
    stats::model.matrix(fml, data = model_frame, rhs = 1, na.action = stats::na.pass),
    error = function(e) stop(sprintf("Formula evaluation failed for data: %s", conditionMessage(e)), call. = FALSE)
  )
  intercept_requested <- isTRUE(attr(rhs_terms, "intercept") == 1) && el_has_explicit_intercept(rhs_expr)
  aux_matrix <- el_drop_intercept_columns(aux_matrix)
  if (intercept_requested) {
    warning("Auxiliary intercepts are ignored; + 1 has no effect.", call. = FALSE)
  }

  aux_expr_vars <- if (!is.null(rhs_expr)) all.vars(rhs_expr) else character(0)
  if (length(aux_expr_vars) > 0 && outcome_source %in% aux_expr_vars) {
    stop("The outcome cannot appear in auxiliary constraints.", call. = FALSE)
  }

  drop_names <- c(outcome_source, outcome_label)
  if (ncol(aux_matrix) > 0) {
    keep <- !colnames(aux_matrix) %in% drop_names
    aux_matrix <- aux_matrix[, keep, drop = FALSE]
  }

  el_validate_matrix_block(
    mat = aux_matrix,
    mask = mask,
    row_map = which(mask),
    label = "Auxiliary covariate",
    severity = "error"
  )

  aux_matrix
}

el_build_missingness_matrix <- function(fml,
                                        model_frame,
                                        n_rhs_parts,
                                        outcome_label,
                                        mask,
                                        respondent_indices) {
  n_resp <- sum(mask)
  if (n_resp == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  intercept_col <- matrix(1, nrow = n_resp, ncol = 1)
  colnames(intercept_col) <- "(Intercept)"

  response_vector <- Formula::model.part(fml, data = model_frame, lhs = 1, drop = TRUE)
  outcome_col <- matrix(response_vector[mask], ncol = 1)
  colnames(outcome_col) <- outcome_label

  rhs_predictors <- if (n_rhs_parts >= 2L) {
    rhs_formula <- stats::formula(fml, lhs = 0, rhs = 2)
    rhs_terms <- tryCatch(
      stats::terms(rhs_formula, data = model_frame),
      error = function(e) stop(sprintf("Formula evaluation failed for data: %s", conditionMessage(e)), call. = FALSE)
    )
    el_assert_no_offset(rhs_terms, "data", "response predictors")
    if (isTRUE(attr(rhs_terms, "intercept") == 0)) {
      warning("Missingness-model intercept is required; '-1' or '+0' has no effect.", call. = FALSE)
    }
    rhs_matrix <- tryCatch(
      stats::model.matrix(fml, data = model_frame, rhs = 2, na.action = stats::na.pass),
      error = function(e) stop(sprintf("Formula evaluation failed for data: %s", conditionMessage(e)), call. = FALSE)
    )
    rhs_matrix <- el_drop_intercept_columns(rhs_matrix)
    el_validate_matrix_block(
      mat = rhs_matrix,
      mask = mask,
      row_map = respondent_indices,
      label = "Missingness-model predictor",
      severity = "warn"
    )
  } else {
    matrix(nrow = n_resp, ncol = 0)
  }

  cbind(intercept_col, outcome_col, rhs_predictors)
}

#' Detect whether the RHS explicitly requests an intercept
#'
#' The Formula/terms machinery already controls implicit intercepts, but we
#' warn users who add `+ 1` manually. This helper scans the parsed language tree
#' for an explicit scalar `1` so we can emit a user-facing warning without
#' treating `.` expansions or other implicit intercepts as explicit requests.
#'
#' @keywords internal
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

#' Remove `(Intercept)` columns from a model matrix
#'
#' Auxiliary constraints should never include an intercept, and the
#' missingness design injects its own intercept. We therefore strip the
#' automatically generated `(Intercept)` column before running downstream
#' validation.
#'
#' @keywords internal
el_drop_intercept_columns <- function(mat) {
  if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0) return(mat)
  keep <- colnames(mat) != "(Intercept)"
  if (all(keep)) return(mat)
  mat[, keep, drop = FALSE]
}

#' Validate a design-matrix block and optionally return respondent rows
#'
#' Applies the respondent mask, enforces the "no NA, no zero-variance" policy,
#' and returns the subset that downstream routines should use.
#'
#' @keywords internal
el_validate_matrix_block <- function(mat, mask, row_map, label, severity = c("error", "warn")) {
  if (is.null(mat) || !is.matrix(mat)) return(mat)
  mat_sub <- if (is.null(mask)) mat else mat[mask, , drop = FALSE]
  if (nrow(mat_sub) == 0 || ncol(mat_sub) == 0) return(mat_sub)
  el_assert_no_na(mat_sub, row_map, label)
  el_check_constant_columns(mat_sub, label = label, severity = severity)
  mat_sub
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

el_outcome_label <- function(expr) {
  paste(deparse(expr), collapse = " ")
}

el_validate_design_spec <- function(design, data_nrow, context_label = "data") {
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
