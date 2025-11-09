#' Blueprint and terms builders
#'
#' Captures formula structure at parse time, including dot expansion, factor
#' level mappings, and the provenance of each RHS symbol. The frozen snapshot
#' ensures reproducible `model.matrix()` construction even if the data object
#' changes later. Variable expansion happens once and is reused for all
#' subsequent design-matrix calls via cached terms, `xlevels`, contrasts, and
#' data vs environment metadata.
#'
#' @keywords internal

nmar_build_terms_info <- function(formula, data, drop_response = FALSE) {
  if (is.null(formula)) {
    return(list(
      formula = NULL,
      terms = NULL,
      column_names = character(),
      source_variables = character(),
      source_variables_data = character(),
      source_variables_env = character(),
      xlevels = NULL,
      contrasts = NULL
    ))
  }
  vars_needed <- unique(setdiff(all.vars(formula), "."))
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  vars_missing <- vars_needed[!vapply(vars_needed, function(var) {
    (var %in% names(data)) || exists(var, envir = env, inherits = TRUE)
  }, logical(1L))]
  if (length(vars_missing) > 0) {
    stop("Variables not found in data: ", paste(vars_missing, collapse = ", "), call. = FALSE)
  }
  mf <- tryCatch(
    stats::model.frame(formula, data = data, na.action = stats::na.pass, drop.unused.levels = FALSE),
    error = function(e) {
      stop(
        "Failed to evaluate formula ",
        paste(deparse(formula, width.cutoff = 500L), collapse = " "),
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  tr <- attr(mf, "terms")
  tr_use <- tr
  if (isTRUE(drop_response) && attr(tr, "response") > 0) {
    tr_use <- stats::delete.response(tr)
  }
  mm <- stats::model.matrix(tr_use, mf)
# Derive RHS source variables from the response-less terms object.
# Using stats::formula(tr_use) ensures we only capture RHS symbols and
# not the LHS outcome. IMPORTANT: Dots (.) are expanded by model.frame at
# this point based on the data supplied, and the expansion is frozen in the
# terms object. We also record which symbols existed in the data so that
# validation can ignore helpers that live solely in the formula environment.
  rhs_vars <- unique(setdiff(all.vars(stats::formula(tr_use)), "."))
  rhs_vars_data <- rhs_vars[rhs_vars %in% names(data)]
  rhs_vars_env <- setdiff(rhs_vars, rhs_vars_data)
  list(
    formula = formula,
    terms = tr_use,
    column_names = colnames(mm),
    source_variables = rhs_vars,
    source_variables_data = rhs_vars_data,
    source_variables_env = rhs_vars_env,
    xlevels = attr(mf, "xlevels"),
    contrasts = attr(mm, "contrasts")
  )
}

#' @keywords internal
nmar_model_matrix_from_terms <- function(terms_info, data) {
  if (is.null(terms_info$terms)) {
    return(NULL)
  }
  mf <- stats::model.frame(
    terms_info$terms,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE,
    xlev = terms_info$xlevels
  )
  stats::model.matrix(terms_info$terms, mf, contrasts.arg = terms_info$contrasts)
}

#' @keywords internal
nmar_build_formula_blueprint <- function(outcome_vars,
                                         aux_expr,
                                         response_expr,
                                         data,
                                         env) {
  aux_formula <- nmar_make_aux_formula(outcome_vars, aux_expr, env)
  aux_terms <- nmar_build_terms_info(aux_formula, data, drop_response = TRUE)

  response_formula <- nmar_make_response_formula(outcome_vars, response_expr, env)
  response_info <- nmar_build_terms_info(response_formula, data, drop_response = TRUE)

  list(
    aux = aux_terms,
    response = response_info,
    outcome = outcome_vars
  )
}
