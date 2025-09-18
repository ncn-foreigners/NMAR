#' @title Exponential Tilting Engine for NMAR Estimation
#'
#' @description Creates a configuration object for the Exponential Tilting (ET) method
#'   used in Not Missing at Random (NMAR) estimation. This function allows users to
#'   specify key parameters for the ET algorithm, such as the probability model type
#'   for missingness, the density model for the outcome variable, and convergence criteria.
#'   The returned object can then be passed to the `nmar()` function to perform
#'   the NMAR estimation using the ET approach.
#'
#' @param family A character string specifying the probability model for
#'   data missingness. Valid options are `"logit"` for a logistic regression model
#'   or `"probit"` for a probit regression model. Defaults to a value loaded from
#'   package settings (likely `"probit"`).
#' @param y_dens A character string specifying the density model for the outcome
#'   variable. Valid options are `"normal"` for a normal distribution or `"gamma"`
#'   for a gamma distribution. Defaults to a value loaded from package settings
#'   (likely `"normal"`).
#' @param tol_value A numeric value representing the tolerance for convergence
#'   of the optimization algorithm. The algorithm stops when the change in
#'   parameters or likelihood falls below this value. Defaults to `0.00001`.
#' @param min_iter An integer specifying the minimum number of iterations the
#'   optimization algorithm will perform, regardless of convergence. Defaults to `10`.
#' @param max_iter An integer specifying the maximum number of iterations the
#'   optimization algorithm will perform. The algorithm will stop if this limit
#'   is reached, even if convergence has not been achieved. Defaults to `100`.
#' @param optim_method A character string specifying the optimization method to
#'   be used for parameter estimation. Valid options include `"Newton"` for
#'   Newton-Raphson or `"Broyden"` for Broyden's method. Defaults to a value
#'   loaded from package settings.
#'
#' @return An engine configuration object of S3 class `c("nmar_engine_exptilt", "nmar_engine")`.
#'   This object encapsulates all the specified parameters for the Exponential
#'   Tilting method and is ready to be used with the `nmar()` function.
#'
#' @importFrom jsonlite read_json
#' @importFrom utils modifyList
#' @export
exptilt_engine <- function(
    standardize=TRUE,
    on_failure = c("return","error"),
    bootstrap_reps=10,
    supress_warnings=FALSE,
    auxiliary_means = NULL,
    control=list(),
    family=c("logit", "probit"),
    y_dens=c("auto","normal", "gamma"),
    variance_method=c("delta","bootstrap"),
    min_iter=10,#todo move to control
    max_iter=100,#todo move to control
    optim_method=c("Newton","Broyden"), #todo move to control
    tol_value=1e-5 #todo move to control

    ) {

  on_failure <- match.arg(on_failure)
  family <- match.arg(family)
  y_dens <- match.arg(y_dens)
  variance_method <- match.arg(variance_method)
  y_dens <- match.arg(y_dens)
  optim_method <- match.arg(optim_method)

  validator$assert_logical(standardize, name = "standardize")
  validator$assert_choice(on_failure, choices = c("return", "error"), name = "on_failure")
  validator$assert_positive_integer(bootstrap_reps, name = "bootstrap_reps",is.finite=TRUE)
  validator$assert_logical(supress_warnings, name = "supress_warnings")
  validator$assert_choice(family, choices = c("logit", "probit"), name = "family")
  validator$assert_choice(y_dens, choices = c("auto", "normal", "gamma"), name = "y_dens")
  validator$assert_choice(variance_method, choices = c("delta", "bootstrap"), name = "variance_method")
  validator$assert_positive_integer(min_iter, name = "min_iter")
  validator$assert_positive_integer(max_iter, name = "max_iter")
  validator$assert_choice(optim_method, choices = c("Newton", "Broyden"), name = "optim_method")
  validator$assert_number(tol_value, name = "tol_value", min = 0, max = Inf)

  engine <- list(
    prob_model_type = family,
    y_dens = y_dens,
    tol_value = tol_value,
    min_iter = min_iter,
    max_iter = max_iter,
    auxiliary_means = auxiliary_means,
    standardize=standardize,
    optim_method = optim_method,
    bootstrap_reps=bootstrap_reps,
    variance_method=variance_method

  )

  class(engine) <- c("nmar_engine_exptilt", "nmar_engine")
  return(engine)
}
