#' Internal: Adapt NMAR formula-list to EL inputs
#' @details Converts the NMAR `formula` list used by `nmar()` into a two-sided
#'   EL formula and a character vector of response-model predictors. Terms in
#'   `covariates_missingness` become response-only predictors; they do not need
#'   to appear on the RHS of the outcome formula (auxiliaries). Only variables
#'   on the outcome RHS are treated as auxiliaries and require known population
#'   means when `auxiliary_means` is supplied. See Qin, Leung and Shao (2002).
#' @keywords internal
nmar_formula_to_el <- function(formula_list) {
  stopifnot(is.list(formula_list))
  lhs <- all.vars(formula_list$outcome)
  if (length(lhs) != 1) stop("Outcome formula must specify exactly one LHS variable.")
  rhs_terms <- if (!is.null(formula_list$covariates_outcome)) {
    all.vars(formula_list$covariates_outcome)
  } else {
    character(0)
  }
  rhs <- if (length(rhs_terms)) paste(rhs_terms, collapse = " + ") else "1"
  el_formula <- stats::as.formula(paste(lhs, "~", rhs))

  miss_terms <- if (!is.null(formula_list$covariates_missingness)) {
    all.vars(formula_list$covariates_missingness)
  } else {
    character(0)
  }
  resp_preds <- if (length(miss_terms)) miss_terms else NULL

  list(formula = el_formula, response_predictors = resp_preds)
}
