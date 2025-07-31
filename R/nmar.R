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
#'
#' @examples
#' # First, create an engine configuration for Exponential Tilting
#' # (assuming 'exptilt' function is defined as below)
#' exptilt_config <- exptilt(
#'   prob_model_type = 'probit',
#'   y_dens = 'normal',
#'   tol_value = 0.00001,
#'   min_iter = 30,
#'   max_iter = 100
#' )
#'
#' # Create some example data
#' set.seed(123)
#' n_obs <- 100
#' example_data <- data.frame(
#'   Y = rnorm(n_obs),
#'   x1 = rnorm(n_obs),
#'   x2 = runif(n_obs)
#' )
#'
#' # Define the formulas
#' model_formula <- list(
#'   outcome = ~ Y,
#'   covariates_outcome = ~ x1 + x2,
#'   covariates_missingness = NULL # Or specify if you have covariates for missingness
#' )
#'
#' # Run the NMAR estimation
#' # Note: 'validate_nmar' and 'run_engine' are internal functions here.
#' # The example focuses on the user-facing call.
#' # res <- nmar(formula = model_formula, data = example_data, engine = exptilt_config)
#'
#' # Due to the internal nature of 'validate_nmar' and 'run_engine' not being available,
#' # the example cannot be fully executable without the full package.
#' # However, this demonstrates the intended usage pattern.
#' # For a real scenario, you would expect 'res' to contain estimation results.
nmar <- function(formula, data, engine) {
  outcome_variable <- all.vars(formula$outcome)
  covariates_for_outcome <- all.vars(formula$covariates_outcome)
  covariates_for_missingness <- all.vars(formula$covariates_missingness)

  # validate_nmar(
  #   data = data,
  #   outcome_variable = outcome_variable,
  #   covariates_for_outcome = covariates_for_outcome,
  #   covariates_for_missingness = covariates_for_missingness
  # )

  if (!inherits(engine, "nmar_engine")) {
    stop("Engine must be created by an engine constructor function")
  }
  res <- run_engine(engine, formula, data)
  # validate_nmar_result(res)

  return(res)

}

run_engine <- function(engine, ...) {
  UseMethod("run_engine")
}

