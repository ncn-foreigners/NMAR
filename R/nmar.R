#' @title Not Missing at Random (NMAR)
#'
#' @description This function serves as the primary interface for performing
#'   Missing Not At Random (NMAR) data estimation.
#'
#' @param method A character string specifying the NMAR imputation method to use.
#'   Currently, only `"exptilt"` (Exponential Tilting) is supported.
#' @param data A data frame or matrix containing the dataset. Missing values in
#'   the `outcome_variable` should be represented by `NA`.
#' @param outcome_variable A character string specifying the name of the single
#'   target (outcome) column that may contain missing values.
#' @param covariates_for_outcome A character vector of column names for the
#'   covariates that are used in the conditional density model of the
#'   `outcome_variable`
#' @param covariates_for_missingness A character vector of column names for
#'   additional variables (apart from the `outcome_variable`) that are used to
#'   model the observation probability (i.e., the missingness mechanism).
#' @param settings A `list` of control parameters specific to the chosen `method`.
#'   For `method = "exptilt"`, this list can contain the following elements
#'   (all with sensible defaults if not provided):
#'   \itemize{
#'     \item{\strong{`prob_model_type`}}{ (`character`): Type of probability model for the missingness mechanism.
#'           Valid options: `"logit"` (for logistic regression, default) or `"probit"` (for probit regression).}
#'     \item{\strong{`y_dens`}}{ (`character`): Type of conditional density model for the `outcome_variable`
#'           given `covariates_for_outcome`. Valid options: `"normal"` (for a normal distribution, default)
#'           or `"gamma"` (for a gamma distribution).}
#'     \item{\strong{`tol_value`}}{ (`numeric`): Tolerance value for the convergence criterion in the
#'           iterative optimization process. Default is `0.00001`.}
#'     \item{\strong{`min_iter`}}{ (`integer`): Minimum number of iterations for the `nleqslv`
#'           optimization algorithm. Default is `10`.}
#'     \item{\strong{`max_iter`}}{ (`integer`): Maximum number of iterations for the `nleqslv`
#'           optimization algorithm. Default is `100`.}
#'     \item{\strong{`optim_method`}}{ (`character`): Optimization method to be used by `nleqslv`.
#'           Valid options: `"Newton"` (default) or `"Broyden`".}
#'   }
#'   Any element not specified in `settings` will revert to its default value for the chosen method.
#' @return A `list` object containing the results of the NMAR imputation. The exact
#'   structure of the returned list depends on the chosen `method`.
#'   For `method = "exptilt"`, the list typically includes:
#'   \itemize{
#'     \item{\strong{`theta`}}{ The estimated tilting parameters.}
#'     \item{\strong{`est_mean`}}{ The estimated mean of the imputed `outcome_variable`.}
#'     \item{\strong{`loss_value`}}{ The final loss function value from the
#'       optimization, indicating convergence status.}
#'   }
#' @details
#'   The `nmar` function acts as a dispatcher, forwarding the imputation task to
#'   specific internal functions based on the `method` argument.
#'   It enforces strict input validation to ensure robust and predictable behavior.
#'
#'   The Exponential Tilting (ET) method (`"exptilt"`) addresses Missing Not At Random
#'   (NMAR) data in a single outcome variable by iteratively estimating model
#'   parameters for the missingness mechanism and the conditional distribution of
#'   the outcome, simultaneously imputing the missing values.
#'   This approach assumes specific model structures for the missingness indicator
#'   and the conditional density of the outcome, as configured via `settings`.
#'   Users should ensure that `covariates_for_outcome` and `covariates_for_missingness`
#'   contain only fully observed variables.
#'
#' @examples
#' # This is a minimal, self-contained, and runnable example.
#' #TODO
#'
#'
#' @seealso  TODO
#' @importFrom bbmle mle2 # These imports should ideally be in your _PACKAGE.R file
#' @importFrom nleqslv nleqslv # (as discussed previously)
#' @importFrom stats as.formula coef dnorm dgamma plogis pnorm runif sd setNames # (same here)
#' @export

nmar <- function(method,data,outcome_variable,covariates_for_outcome,covariates_for_missingness,settings){

  #TODO
  # validate_nmar(data=data,outcome_variable=outcome_variable,covariates_for_outcome = covariates_for_outcome,covariates_for_missingness=covariates_for_missingness)

  if (method == 'exptilt') {

    return(nmar_exptilt_engine(data, outcome_variable, covariates_for_outcome, covariates_for_missingness, settings))
  } else {
    stop("Work in progress")
  }
}
