#' @title Not Missing at Random (NMAR) Estimation
#'
#' @description Unified interface for NMAR estimation using configurable engines.
#'
#' @param formula Named list of formulas specifying model components
#' @param data Data frame containing the variables
#' @param engine Engine configuration object created by engine constructor
#' @return Estimation results specific to the engine
#' @export
nmar <- function(formula, data, engine) {
  outcome_variable <- all.vars(formula$outcome)
  covariates_for_outcome <- all.vars(formula$covariates_outcome)
  covariates_for_missingness <- all.vars(formula$covariates_missingness)

  validate_nmar(
    data = data,
    outcome_variable = outcome_variable,
    covariates_for_outcome = covariates_for_outcome,
    covariates_for_missingness = covariates_for_missingness
  )

  if (!inherits(engine, "nmar_engine")) {
    stop("Engine must be created by an engine constructor function")
  }
  run_engine(engine, formula, data)
}

run_engine <- function(engine, ...) {
  UseMethod("run_engine")
}

