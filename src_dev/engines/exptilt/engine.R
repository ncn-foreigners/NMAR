#' @title Exponential Tilting Engine for NMAR Estimation
#'
#' @description Creates a configuration object for the Exponential Tilting (ET) method
#'   used in Not Missing at Random (NMAR) estimation. This function allows users to
#'   specify key parameters for the ET algorithm, such as the probability model type
#'   for missingness, the density model for the outcome variable, and convergence criteria.
#'   The returned object can then be passed to the `nmar()` function to perform
#'   the NMAR estimation using the ET approach.
#'
#' @param prob_model_type A character string specifying the probability model for
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
#'
#' @examples
#' # Create a default Exponential Tilting engine configuration
#' default_exptilt_config <- exptilt()
#' print(default_exptilt_config)
#'
#' # Create an Exponential Tilting engine with custom parameters
#' custom_exptilt_config <- exptilt(
#'   prob_model_type = 'logit',
#'   y_dens = 'gamma',
#'   tol_value = 1e-6,
#'   min_iter = 50,
#'   max_iter = 200,
#'   optim_method = 'Broyden'
#' )
#' print(custom_exptilt_config)
#'
#' # This configuration object can then be passed to the 'nmar' function:
#' # nmar_results <- nmar(formula = my_formula, data = my_data, engine = custom_exptilt_config)
exptilt <- function(
    prob_model_type = get_json_param_info(all_schemas, "exptilt", "prob_model_type")$default,
    y_dens = get_json_param_info(all_schemas, "exptilt", "y_dens")$default,
    tol_value = get_json_param_info(all_schemas, "exptilt", "tol_value")$default,
    min_iter = get_json_param_info(all_schemas, "exptilt", "min_iter")$default,
    max_iter = get_json_param_info(all_schemas, "exptilt", "max_iter")$default,
    optim_method = get_json_param_info(all_schemas, "exptilt", "optim_method")$default) {

  # Load the JSON file inside the function
  json_path <- system.file("extdata", "method_params.json", package = "nmar")

  if (!file.exists(json_path)) {
    stop("Method parameters JSON file (method_params.json) not found in 'inst/extdata/'. Package might be corrupted or not installed correctly.")
  }

  all_schemas <- jsonlite::read_json(json_path, simplifyVector = TRUE)

  config <- list(
    prob_model_type = prob_model_type,
    y_dens = y_dens,
    tol_value = tol_value,
    min_iter = min_iter,
    max_iter = max_iter,
    optim_method = optim_method
  )

  #TODO
  # config <- validate_method_settings(
  #   method_name = "exptilt",
  #   arguments = config,
  #   schemas = all_schemas
  # )

  class(config) <- c("nmar_engine_exptilt", "nmar_engine")
  return(config)
}
