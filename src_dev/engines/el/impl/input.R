#' Prepare inputs for EL estimation
#' @details Uses the Formula package to parse the user formula, enforce EL
#'   constraints (single outcome, no outcome on auxiliary side), and return
#'   normalized formulas plus the augmented data set containing the delta
#'   indicator.
#' @keywords internal
el_prepare_inputs <- function(formula, data, require_na = TRUE, auxiliary_means = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  }
  outcome_sym <- formula[[2L]]
  outcome_vars <- all.vars(outcome_sym)
  if (length(outcome_vars) != 1L) {
    stop("The left-hand side must contain exactly one outcome variable.", call. = FALSE)
  }
  outcome_var <- outcome_vars[1]
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()

  el_validate_outcome(data, outcome_var, require_na)
  respondent_mask <- !is.na(data[[outcome_var]])

  formula_obj <- Formula::Formula(stats::as.formula(formula))
  rhs_parts <- length(attr(formula_obj, "rhs"))
  if (rhs_parts > 2L) {
    stop("EL formulas support at most two RHS partitions (auxiliaries | response predictors).", call. = FALSE)
  }
  aux_rhs_formula <- if (rhs_parts >= 1L) {
    stats::formula(formula_obj, lhs = 0, rhs = 1L)
  } else {
    stats::as.formula("~ 0")
  }
  resp_rhs_formula <- if (rhs_parts >= 2L) stats::formula(formula_obj, lhs = 0, rhs = 2L) else NULL

  norm <- el_aux_formula_normalize(aux_rhs_formula, data = data, env = env, outcome_var = outcome_var)
  auxiliary_formula <- norm$formula_no_int
  if (isTRUE(norm$removed_intercept)) {
    warning(
      "Auxiliary constraints do not include an intercept; dropping the requested intercept (use '~ 0 + ...').",
      call. = FALSE
    )
  }
  el_validate_aux_no_na(
    data = data,
    aux_formula = auxiliary_formula,
    respondent_mask = respondent_mask,
    auxiliary_means = auxiliary_means
  )

# Response-model predictors must not contain NA among respondents.
  resp_predictor_expr <- NULL
  if (!is.null(resp_rhs_formula)) {
    resp_rhs_terms <- tryCatch(
      stats::terms(resp_rhs_formula, data = data),
      error = el_rethrow_data_error
    )
    resp_rhs_expanded <- stats::formula(resp_rhs_terms)
    resp_predictor_expr <- resp_rhs_expanded[[2L]]
  }
  rhs_resp_expr <- if (is.null(resp_predictor_expr)) outcome_sym else call("+", outcome_sym, resp_predictor_expr)
  resp_check_fml <- as.formula(call("~", call("+", 0, rhs_resp_expr)))
  environment(resp_check_fml) <- env
  resp_mf <- tryCatch(
    stats::model.frame(resp_check_fml, data = data, na.action = stats::na.pass),
    error = el_rethrow_data_error
  )
  if (ncol(resp_mf) > 0) {
    resp_mf_resp <- resp_mf[respondent_mask, , drop = FALSE]
    if (anyNA(resp_mf_resp)) {
      na_by_col <- vapply(seq_len(ncol(resp_mf_resp)), function(j) anyNA(resp_mf_resp[[j]]), logical(1))
      bad_col <- colnames(resp_mf_resp)[which(na_by_col)[1]]
      bad_row_local <- which(is.na(resp_mf_resp[[bad_col]]))[1]
      bad_row <- which(respondent_mask)[bad_row_local]
      msg <- sprintf(
        "Response-model predictor '%s' contains NA values among respondents.%s",
        bad_col,
        if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else ""
      )
      stop(msg, call. = FALSE)
    }
  }

  d <- el_make_delta_column(data, outcome_var)
  outcome_fml <- as.formula(call("~", outcome_sym, 1L))
  environment(outcome_fml) <- env
  response_fml <- as.formula(call("~", as.name(d$delta_name), rhs_resp_expr))
  environment(response_fml) <- env

  list(
    data = d$data,
    formula_list = list(
      outcome = outcome_fml,
      response = response_fml,
      auxiliary = auxiliary_formula
    )
  )
}

#' Normalize auxiliary RHS to a one-sided, no-intercept formula
#' @keywords internal
el_aux_formula_normalize <- function(aux_formula, data, env, outcome_var) {
  stopifnot(!is.null(env))
  if (is.null(aux_formula)) {
    return(list(formula_no_int = NULL, removed_intercept = FALSE))
  }
  aux_formula <- stats::as.formula(aux_formula)
  environment(aux_formula) <- env
  expr_vars <- all.vars(aux_formula)
  uses_dot <- "." %in% expr_vars
  if (outcome_var %in% expr_vars) {
    stop("The outcome cannot appear in auxiliary constraints.", call. = FALSE)
  }
  terms_obj <- tryCatch(
    stats::terms(aux_formula, data = data),
    error = el_rethrow_data_error
  )
  term_labels <- attr(terms_obj, "term.labels")
  if (uses_dot) {
    term_labels <- setdiff(term_labels, outcome_var)
  }
  has_terms <- length(term_labels) > 0
  expr <- aux_formula[[2L]]
  explicit_intercept <- el_expr_has_explicit_intercept(expr)
  removed_intercept <- explicit_intercept && isTRUE(attr(terms_obj, "intercept") == 1) && has_terms
  if (!has_terms) {
    return(list(formula_no_int = NULL, removed_intercept = FALSE))
  }
  aux_no_int <- stats::reformulate(term_labels, intercept = FALSE)
  environment(aux_no_int) <- env
  list(formula_no_int = aux_no_int, removed_intercept = removed_intercept)
}

#' Detect explicit inclusion of an intercept (e.g., "+ 1") in RHS expressions
#' @keywords internal
el_expr_has_explicit_intercept <- function(node) {
  if (is.null(node)) return(FALSE)
  if (is.numeric(node) && length(node) == 1L && isTRUE(all.equal(node, 1))) return(TRUE)
  if (is.call(node)) {
    op <- as.character(node[[1L]])
    if (op %in% c("+", "-", "~")) {
      for (i in 2L:length(node)) if (Recall(node[[i]])) return(TRUE)
    }
  }
  FALSE
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

#' Re-throw parsing/model-frame errors with a standardized message
#' @keywords internal
el_rethrow_data_error <- function(err) {
  msg <- conditionMessage(err)
  missing_pattern <- "object '([^']+)' not found"
  if (grepl(missing_pattern, msg, perl = TRUE)) {
    missing_var <- sub(missing_pattern, "\\1", msg, perl = TRUE)
    stop(sprintf("Variables not found in data: %s", missing_var), call. = FALSE)
  }
  stop(msg, call. = FALSE)
}
