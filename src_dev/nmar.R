#' @title Not Missing at Random (NMAR) Estimation
#'
#' @description Provides a unified interface for Not Missing at Random (NMAR) estimation.
#'   This function orchestrates the estimation process by validating inputs and
#'   dispatching to the appropriate engine based on the provided `engine` object.
#'   It ensures that all necessary data and model specifications are correctly
#'   formatted before computation begins.
#'
#' @param formula A named list of formulas specifying the model components.
#'   It must contain:
#'   \itemize{
#'     \item `outcome`: A formula for the outcome variable (e.g., `~ Y`).
#'     \item `covariates_outcome`: A formula for covariates affecting the outcome (e.g., `~ x1 + x2`).
#'     \item `covariates_missingness`: A formula for covariates affecting the missingness mechanism.
#'       Use `NULL` if no specific covariates are used for missingness (e.g., in some NMAR models,
#'       missingness depends implicitly on the outcome itself).
#'   }
#' @param data A data frame containing all variables specified in the `formula`.
#' @param engine An engine configuration object, typically created by an
#'   engine constructor function like `exptilt()`. This object defines the
#'   specific NMAR estimation method and its parameters. It must inherit from
#'   class `nmar_engine`.
#'
#' @return An object containing the estimation results, whose structure will be
#'   specific to the `engine` used. This might include estimated parameters,
#'   convergence information, and other relevant output from the chosen NMAR method.
#'
#' @export
nmar <- function(formula, data, engine,response_predictors=NULL) {
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(engine, "nmar_engine"))

  outcome_variable <- as.vector(all.vars(formula[[2]]))
  covariates_for_outcome <- as.vector(all.vars(formula[[3]]))

  # validate_nmar(
  #   data = data,
  #   outcome_variable = outcome_variable,
  #   covariates_for_outcome = covariates_for_outcome,
  #   covariates_for_missingness = covariates_for_missingness
  # )


  res <- run_engine(engine, formula, data,response_predictors)
  # validate_nmar_result(res)

  return(res)

}

run_engine <- function(engine, formula, data,response_predictors) {
  UseMethod("run_engine")
}
