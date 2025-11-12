el_parse_formula <- function(formula, data, require_na = TRUE) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  }

  outcome_expr <- formula[[2L]]
  outcome_vars <- all.vars(outcome_expr)
  if (length(outcome_vars) != 1L) {
    stop("The left-hand side must contain exactly one outcome variable.", call. = FALSE)
  }

  outcome_var <- outcome_vars[1]
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()

  el_validate_outcome(data, outcome_var, require_na)

  formula_obj <- Formula::as.Formula(as.list(formula))
  rhs_attr <- attr(formula_obj, "rhs")
  rhs_parts <- if (is.null(rhs_attr)) 1L else length(rhs_attr)
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
  respondent_indices <- which(respondent_mask)

  list(
    formula = formula,
    formula_obj = formula_obj,
    env = env,
    outcome_var = outcome_var,
    model_frame = model_frame,
    response_vector = response_vector,
    respondent_mask = respondent_mask,
    respondent_indices = respondent_indices,
    rhs_parts = rhs_parts
  )
}

#' Build EL specification via Formula workflows
#'
#' @description
#'   Parses the user formula with `Formula::as.Formula()`, builds the canonical
#'   model frame, constructs response and auxiliary design matrices, guards
#'   against unsupported terms (e.g., outcome on the auxiliary side), and
#'   augments the original data with the delta indicator used throughout the
#'   engine.
#'
#' @keywords internal
el_build_design <- function(parsed_spec, auxiliary_means = NULL) {
  response_matrix <- el_construct_response_matrix(
    response_vector = parsed_spec$response_vector,
    respondent_mask = parsed_spec$respondent_mask,
    formula_obj = parsed_spec$formula_obj,
    rhs_parts = parsed_spec$rhs_parts,
    model_frame = parsed_spec$model_frame,
    outcome_var = parsed_spec$outcome_var
  )
  el_validate_response_matrix(response_matrix, parsed_spec$respondent_indices)

  aux_design <- el_extract_auxiliary_design(
    formula_obj = parsed_spec$formula_obj,
    rhs_parts = parsed_spec$rhs_parts,
    model_frame = parsed_spec$model_frame,
    respondent_mask = parsed_spec$respondent_mask,
    outcome_var = parsed_spec$outcome_var,
    env = parsed_spec$env
  )

  el_validate_auxiliary_matrices(
    aux_resp = aux_design$respondent,
    aux_full = aux_design$full,
    respondent_indices = parsed_spec$respondent_indices,
    need_full_means = is.null(auxiliary_means)
  )

  structure(
    list(
      response = response_matrix,
      aux_resp = aux_design$respondent,
      aux_full = aux_design$full,
      mask = parsed_spec$respondent_mask,
      outcome = parsed_spec$outcome_var
    ),
    class = "el_design"
  )
}

el_construct_response_matrix <- function(response_vector,
                                         respondent_mask,
                                         formula_obj,
                                         rhs_parts,
                                         model_frame,
                                         outcome_var) {
  n_resp <- sum(respondent_mask)
  if (n_resp == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }

  intercept_col <- matrix(1, nrow = n_resp, ncol = 1)
  colnames(intercept_col) <- "(Intercept)"

  outcome_col <- matrix(response_vector[respondent_mask], ncol = 1)
  colnames(outcome_col) <- outcome_var

  resp_terms <- if (rhs_parts >= 2L) {
    mm <- tryCatch(
      stats::model.matrix(formula_obj, data = model_frame, rhs = 2L),
      error = el_rethrow_data_error
    )
    mm_resp <- mm[respondent_mask, , drop = FALSE]
    if ("(Intercept)" %in% colnames(mm_resp)) {
      mm_resp <- mm_resp[, colnames(mm_resp) != "(Intercept)", drop = FALSE]
    }
    mm_resp
  } else {
    matrix(nrow = n_resp, ncol = 0)
  }

  cbind(intercept_col, outcome_col, resp_terms)
}

el_extract_auxiliary_design <- function(formula_obj,
                                        rhs_parts,
                                        model_frame,
                                        respondent_mask,
                                        outcome_var,
                                        env) {
  n_total <- nrow(model_frame)
  n_resp <- sum(respondent_mask)

  if (n_total == 0) {
    return(list(
      full = matrix(nrow = 0, ncol = 0),
      respondent = matrix(nrow = 0, ncol = 0),
      has_aux = FALSE
    ))
  }

  if (rhs_parts >= 1L) {
    aux_formula <- stats::formula(formula_obj, lhs = 0, rhs = 1L)
  } else {
    aux_formula <- stats::as.formula("~ 0")
  }

  aux_terms <- tryCatch(
    stats::terms(aux_formula, data = model_frame),
    error = el_rethrow_data_error
  )
  term_labels <- attr(aux_terms, "term.labels")
  outcome_explicit <- el_expr_contains_var(aux_formula[[2L]], outcome_var)
  term_has_outcome <- if (length(term_labels) > 0) {
    vapply(
      term_labels,
      function(label) {
        vars <- all.vars(stats::as.formula(paste("~", label), env = env))
        outcome_var %in% vars
      },
      logical(1)
    )
  } else {
    logical(0)
  }

  if (outcome_explicit && any(term_has_outcome)) {
    bad_terms <- term_labels[term_has_outcome]
    stop(
      sprintf(
        "The outcome cannot appear in auxiliary constraints (invalid terms: %s).",
        paste(bad_terms, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  aux_full <- tryCatch(
    stats::model.matrix(formula_obj, data = model_frame, rhs = 1L),
    error = el_rethrow_data_error
  )

  if (ncol(aux_full) > 0 && outcome_var %in% colnames(aux_full)) {
    aux_full <- aux_full[, colnames(aux_full) != outcome_var, drop = FALSE]
  }

  intercept_present <- "(Intercept)" %in% colnames(aux_full)
  if (intercept_present) {
    aux_full <- aux_full[, colnames(aux_full) != "(Intercept)", drop = FALSE]
  }

  explicit_intercept <- el_formula_has_explicit_intercept(aux_formula[[2L]])
  if (intercept_present && explicit_intercept && ncol(aux_full) > 0) {
    warning(
      "Auxiliary constraints do not include an intercept; dropping the requested intercept (use '~ 0 + ...').",
      call. = FALSE
    )
  }

  has_aux <- ncol(aux_full) > 0

  if (ncol(aux_full) == 0) {
    return(list(
      full = matrix(nrow = n_total, ncol = 0),
      respondent = matrix(nrow = n_resp, ncol = 0)
    ))
  }

  aux_resp <- aux_full[respondent_mask, , drop = FALSE]
  list(
    full = aux_full,
    respondent = aux_resp
  )
}

el_formula_has_explicit_intercept <- function(node) {
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

el_expr_contains_var <- function(node, target_var) {
  if (is.null(node)) return(FALSE)
  if (is.symbol(node)) return(identical(as.character(node), target_var))
  if (is.call(node)) {
    op <- as.character(node[[1L]])
    if (op == ".") return(FALSE)
    for (i in seq_along(node)[-1]) {
      if (el_expr_contains_var(node[[i]], target_var)) return(TRUE)
    }
  }
  FALSE
}

el_validate_response_matrix <- function(response_matrix, respondent_indices) {
  if (!is.matrix(response_matrix) || ncol(response_matrix) == 0 || nrow(response_matrix) == 0) {
    return(invisible(NULL))
  }
  if (!anyNA(response_matrix)) return(invisible(NULL))

  na_loc <- which(is.na(response_matrix), arr.ind = TRUE)[1, , drop = TRUE]
  bad_col <- colnames(response_matrix)[na_loc[2]]
  bad_row <- if (length(respondent_indices) >= na_loc[1]) respondent_indices[na_loc[1]] else NA_integer_

  msg <- sprintf(
    "Response-model predictor '%s' contains NA values among respondents.%s",
    bad_col,
    if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else ""
  )
  stop(msg, call. = FALSE)
}

el_validate_auxiliary_matrices <- function(aux_resp,
                                           aux_full,
                                           respondent_indices,
                                           need_full_means) {
  if (is.null(aux_resp) || ncol(aux_resp) == 0) return(invisible(NULL))

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

  if (isTRUE(need_full_means) && anyNA(aux_full)) {
    na_loc <- which(is.na(aux_full), arr.ind = TRUE)[1, , drop = TRUE]
    bad_col <- colnames(aux_full)[na_loc[2]]
    bad_row <- na_loc[1]
    msg <- sprintf(
      "Auxiliary variables contain NA values in full data (needed to compute population means).%s\nEither provide auxiliary_means=... or remove NA values from auxiliary variables.",
      if (is.finite(bad_row)) sprintf("\nFirst NA in '%s' at row %d", bad_col, bad_row) else ""
    )
    stop(msg, call. = FALSE)
  }
  invisible(NULL)
}

el_make_delta_column <- function(data, outcome_var, respondent_mask = NULL) {
  if (is.null(respondent_mask)) {
    respondent_mask <- !is.na(data[[outcome_var]])
  }
  if (length(respondent_mask) != nrow(data)) {
    stop("Internal error: respondent mask must align with data.", call. = FALSE)
  }

  delta_name <- "..nmar_delta.."
  if (delta_name %in% names(data)) {
    i <- 1L
    while (paste0(delta_name, i) %in% names(data)) i <- i + 1L
    delta_name <- paste0(delta_name, i)
  }

  data[[delta_name]] <- as.integer(respondent_mask)
  list(data = data, delta_name = delta_name)
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
