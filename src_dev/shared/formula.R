#' Formula partitioning and builders
#' @keywords internal

nmar_partition_rhs <- function(rhs_expr) {
  aux_expr <- rhs_expr
  response_expr <- NULL
  if (is.call(rhs_expr) && identical(rhs_expr[[1L]], as.name("|"))) {
    aux_expr <- rhs_expr[[2L]]
    response_expr <- rhs_expr[[3L]]
  }
  list(aux_expr = aux_expr, response_expr = response_expr)
}

#' Normalize auxiliary RHS expression according to engine traits
#'
#' Removes the intercept (when requested) using the original data so that dot
#' expansion matches `model.matrix()` behavior. Outcome variables are temporarily
#' placed on the LHS to prevent them from reappearing on the auxiliary side when
#' users rely on `.`.
#'
#' @keywords internal
nmar_normalize_auxiliary_expr <- function(aux_expr,
                                          drop_intercept = TRUE,
                                          data = NULL,
                                          env = parent.frame(),
                                          exclude = character()) {
  if (is.null(aux_expr)) {
    return(NULL)
  }
  if (!isTRUE(drop_intercept)) {
    return(aux_expr)
  }
  build_lhs <- function(vars) {
    if (!length(vars)) return(NULL)
    lhs <- as.name(vars[[1L]])
    if (length(vars) == 1L) return(lhs)
    for (nm in vars[-1L]) {
      lhs <- call("+", lhs, as.name(nm))
    }
    lhs
  }
  lhs_expr <- build_lhs(exclude)
  formula_call <- if (is.null(lhs_expr)) call("~", aux_expr) else call("~", lhs_expr, aux_expr)
  aux_formula <- stats::as.formula(formula_call)
  attr(aux_formula, ".Environment") <- env
# Prefer evaluating terms() with the original data so dot expansion follows
# the same rules as model.frame(); fall back to a data-free parse when needed
# (e.g., during early validation) to keep the workflow robust.
  terms_obj <- tryCatch(
    stats::terms(aux_formula, data = data, simplify = TRUE),
    error = function(e) stats::terms(aux_formula, simplify = TRUE)
  )
  attr(terms_obj, "intercept") <- 0L
  aux_formula_clean <- stats::formula(terms_obj)
  attr(aux_formula_clean, ".Environment") <- env
  rhs_index <- if (length(aux_formula_clean) >= 3L) 3L else 2L
  rhs_lang <- if (length(aux_formula_clean) >= rhs_index) aux_formula_clean[[rhs_index]] else NULL
  if (is.null(rhs_lang) || identical(rhs_lang, 0) || identical(rhs_lang, 0L)) {
    return(NULL)
  }
  rhs_lang
}

#' Rebuild a partitioned formula y ~ aux | response
#'
#' Constructs a partitioned formula from a base formula whose RHS is auxiliaries
#' only. Supplying `response_rhs_lang` inserts the `| response` partition. The
#' formula environment is preserved. No string manipulation is used; everything
#' is built with language objects. This keeps transforms/factors intact.
#' @keywords internal
nmar_rebuild_partitioned_formula <- function(base_formula,
                                            response_rhs_lang = NULL,
                                            aux_rhs_lang = NULL,
                                            env = NULL) {
  if (!inherits(base_formula, "formula") || length(base_formula) != 3L) {
    stop("`base_formula` must be a two-sided formula.", call. = FALSE)
  }
  lhs <- base_formula[[2L]]
  rhs_aux <- aux_rhs_lang %||% base_formula[[3L]]
  if (is.null(response_rhs_lang)) {
    f <- call("~", lhs, rhs_aux)
  } else {
    rhs_bar <- call("|", rhs_aux, response_rhs_lang)
    f <- call("~", lhs, rhs_bar)
  }
  if (is.null(env)) env <- environment(base_formula)
  if (is.null(env)) env <- parent.frame()
  class(f) <- "formula"
  attr(f, ".Environment") <- env
  f
}

#' Split a partitioned formula into components
#' @keywords internal
nmar_split_partitioned_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula.", call. = FALSE)
  }
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  outcome_var <- all.vars(formula[[2L]])
  if (length(outcome_var) != 1L) stop("LHS must be a single outcome variable.", call. = FALSE)

  rhs_parts <- nmar_partition_rhs(formula[[3L]])
  list(
    outcome_var = outcome_var[[1L]],
    aux_rhs_lang = rhs_parts$aux_expr,
    response_rhs_lang = rhs_parts$response_expr,
    env = env
  )
}

#' Build internal EL formulas (outcome, response, auxiliary)
#' @keywords internal
nmar_build_internal_formulas <- function(delta_name, outcome_var, aux_rhs_lang, response_rhs_lang, env) {
# Outcome model (used for outcome-side design): y ~ 1
  outcome_fml <- stats::as.formula(call("~", as.name(outcome_var), 1L))
  attr(outcome_fml, ".Environment") <- env
# Response model (missingness): delta ~ outcome (+ explicit RHS if present).
# The response intercept is handled at the design-matrix stage by
# model.matrix defaults or engine-specific logic; we do not inject a
# constant explicitly here.
  rhs_resp <- if (is.null(response_rhs_lang)) {
    as.name(outcome_var)
  } else {
    call("+", as.name(outcome_var), response_rhs_lang)
  }
  response_fml <- stats::as.formula(call("~", as.name(delta_name), rhs_resp))
  attr(response_fml, ".Environment") <- env
  auxiliary_fml <- NULL
  if (!is.null(aux_rhs_lang)) {
# Build auxiliaries without an intercept to match population-moment logic.
# Using 0 + RHS is the base-R idiom to suppress the intercept in a formula.
    rhs_aux <- call("+", 0, aux_rhs_lang)
    auxiliary_fml <- stats::as.formula(call("~", rhs_aux))
    attr(auxiliary_fml, ".Environment") <- env
  }
  list(outcome = outcome_fml, response = response_fml, auxiliary = auxiliary_fml)
}

#' Generate a unique column name that does not collide with names(data)
#' @keywords internal
nmar_make_unique_colname <- function(base, data_names) {
  nm <- base
  if (!(nm %in% data_names)) return(nm)
  i <- 1L
  repeat {
    cand <- paste0(base, i)
    if (!(cand %in% data_names)) return(cand)
    i <- i + 1L
  }
}

#' Build auxiliary and response formulas for blueprints
#' @keywords internal
nmar_make_aux_formula <- function(outcome_vars, aux_expr, env) {
  if (is.null(aux_expr)) return(NULL)
  if (length(all.vars(aux_expr)) == 0L) return(NULL)
  lhs <- nmar_make_lhs_expr(outcome_vars)
  rhs <- call("+", 0, aux_expr)
  f <- if (is.null(lhs)) call("~", rhs) else call("~", lhs, rhs)
  stats::as.formula(f, env = env)
}

#' @keywords internal
nmar_make_response_formula <- function(outcome_vars, response_expr, env) {
  if (is.null(response_expr)) return(NULL)
  if (length(all.vars(response_expr)) == 0L) return(NULL)
  lhs <- nmar_make_lhs_expr(outcome_vars)
  rhs <- call("+", 0, response_expr)
  f <- if (is.null(lhs)) call("~", rhs) else call("~", lhs, rhs)
  stats::as.formula(f, env = env)
}

#' @keywords internal
nmar_make_lhs_expr <- function(outcome_vars) {
  if (length(outcome_vars) == 0L) return(NULL)
  lhs <- as.name(outcome_vars[[1L]])
  if (length(outcome_vars) == 1L) return(lhs)
  for (var in outcome_vars[-1L]) {
    lhs <- call("+", lhs, as.name(var))
  }
  lhs
}
