#' Minimal input helpers (internal)
#'
#' Tiny, focused utilities to support a formula-first API while engines
#' remain in charge of modeling semantics. Keep helpers minimal and only where
#' they add clarity or testability.
#'
#' @keywords internal

#' Build an auxiliary formula from RHS language; kept for testing utilities.
#' Returns NULL if `aux_rhs_lang` is NULL or has no variables.
nmar_make_aux_formula <- function(outcome, aux_rhs_lang, env = parent.frame()) {
  if (is.null(aux_rhs_lang) || length(all.vars(aux_rhs_lang)) == 0) return(NULL)
  rhs <- call("+", 0, aux_rhs_lang)
  f <- stats::as.formula(call("~", rhs))
  environment(f) <- env
  f
}
