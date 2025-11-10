#' Minimal input helpers (internal)
#'
#' Tiny, focused utilities to support a single formula-first API while engines
#' remain in charge of semantics. These helpers do just enough to build model
#' matrices via base R (model.frame/terms/model.matrix) without blueprints,
#' environment freezing, or traits. Keep them small and predictable.
#'
#' @keywords internal

# Split a RHS possibly containing a partition bar `|` into
# auxiliaries (left of `|`) and response-only predictors (right of `|`).
nmar_split_rhs <- function(rhs_expr) {
  aux_expr <- rhs_expr
  resp_expr <- NULL
  if (is.call(rhs_expr) && identical(rhs_expr[[1L]], as.name("|"))) {
    aux_expr <- rhs_expr[[2L]]
    resp_expr <- rhs_expr[[3L]]
  }
  list(aux = aux_expr, resp = resp_expr)
}

# Evaluate the LHS outcome expression once using model.frame(~ expr, na.pass).
# Returns a numeric vector or a 1-column matrix (treated as vector), the
# canonical outcome column name, and the (possibly updated) data object in case
# of survey designs where variables live under $variables.
nmar_eval_lhs <- function(formula, data) {
  stopifnot(inherits(formula, "formula"), length(formula) == 3L)
  lhs <- formula[[2L]]
  env <- environment(formula) %||% parent.frame()

# Fast-path: symbol pointing to a data column
  if (!is.call(lhs) && is.symbol(lhs)) {
    nm <- as.character(lhs)
    if (inherits(data, "survey.design")) {
      if (!nm %in% names(data$variables)) stop("Outcome variable not found in design$variables.")
      y <- data$variables[[nm]]
      return(list(y = as.numeric(y), outcome_name = nm, data = data))
    } else {
      if (!nm %in% names(data)) stop("Outcome variable not found in data.")
      y <- data[[nm]]
      return(list(y = as.numeric(y), outcome_name = nm, data = data))
    }
  }

# General case: evaluate transform via model.frame(~ expr)
  f <- stats::as.formula(call("~", lhs), env = env)
  mf <- if (inherits(data, "survey.design")) {
    stats::model.frame(f, data = data$variables, na.action = stats::na.pass, drop.unused.levels = FALSE)
  } else {
    stats::model.frame(f, data = data, na.action = stats::na.pass, drop.unused.levels = FALSE)
  }
  if (ncol(mf) == 0) stop("Outcome expression did not produce any values.")
  val <- mf[[1L]]
  if (!is.numeric(val)) stop("Outcome transformation must evaluate to numeric.")
# Assign a stable internal column name for transformed outcomes
  outcome_name <- "..nmar_outcome.."
  if (inherits(data, "survey.design")) {
    data$variables[[outcome_name]] <- as.numeric(val)
  } else {
    data[[outcome_name]] <- as.numeric(val)
  }
  list(y = as.numeric(val), outcome_name = outcome_name, data = data)
}

# Build a model matrix from a RHS language (may be NULL). Intercept control via
# adding + 0 when intercept = FALSE. Returns NULL when expr is NULL or when the
# resulting matrix has zero columns.
nmar_mm <- function(expr, data, intercept = FALSE) {
  if (is.null(expr)) return(NULL)
  rhs <- if (isTRUE(intercept)) expr else call("+", 0, expr)
  f <- stats::as.formula(call("~", rhs))
  mf <- if (inherits(data, "survey.design")) {
    stats::model.frame(f, data = data$variables, na.action = stats::na.pass, drop.unused.levels = FALSE)
  } else {
    stats::model.frame(f, data = data, na.action = stats::na.pass, drop.unused.levels = FALSE)
  }
  mm <- stats::model.matrix(attr(mf, "terms"), mf)
  if (is.null(mm) || ncol(mm) == 0) return(NULL)
  mm
}

# Build an auxiliary formula from outcome and aux RHS language; kept for
# testing utilities and internal callers that need a standalone formula.
# Returns NULL if `aux_rhs_lang` is NULL or has no variables.
nmar_make_aux_formula <- function(outcome, aux_rhs_lang, env = parent.frame()) {
  if (is.null(aux_rhs_lang) || length(all.vars(aux_rhs_lang)) == 0) return(NULL)
  rhs <- call("+", 0, aux_rhs_lang)
  f <- stats::as.formula(call("~", rhs))
  environment(f) <- env
  f
}

# Retrieve design weights: for survey designs, stats::weights(design);
# for data.frames, a vector of ones of length nrow(data).
nmar_get_weights <- function(data) {
  if (inherits(data, "survey.design")) {
    w <- stats::weights(data)
    if (is.null(w)) stop("Unable to retrieve weights from survey design.")
    return(as.numeric(w))
  }
  rep(1, nrow(data))
}
