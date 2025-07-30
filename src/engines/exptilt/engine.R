#' Exponential Tilting Engine for NMAR Estimation
#'
#' @param prob_model_type Probability model for data missingness ("logit" or "probit")
#' @param y_dens Density model for outcome variable ("normal" or "gamma")
#' @param tol_value Tolerance in convergence (default = 0.00001)
#' @param min_iter Minimum number of iterations (default = 10)
#' @param max_iter Maximum number of iterations (default = 100)
#' @param optim_method Optimization method ("Newton" or "Broyden")
#' @return Engine configuration object of class `nmar_engine_exptilt`
#' @export

#' @importFrom jsonlite read_json
#' @importFrom utils modifyList

#load inst/extdata/method_params.json
json_path <- system.file("extdata", "method_params.json", package = "nmar")

if (!file.exists(json_path)) {
  stop("Method parameters JSON file (method_params.json) not found in 'inst/extdata/'. Package might be corrupted or not installed correctly.")
}

all_schemas <- jsonlite::read_json(json_path, simplifyVector = TRUE)

exptilt <- function(
    prob_model_type = get_json_param_info(all_schemas, "exptilt", "prob_model_type")$default,
    y_dens =  get_json_param_info(all_schemas, "exptilt", "y_dens")$default,
    tol_value = get_json_param_info(all_schemas, "exptilt", "tol_value")$default,
    min_iter =  get_json_param_info(all_schemas, "exptilt", "min_iter")$default,
    max_iter =  get_json_param_info(all_schemas, "exptilt", "max_iter")$default,
    optim_method = get_json_param_info(all_schemas, "exptilt", "optim_method")$default) {

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
  # browser()
  class(config) <- c("nmar_engine_exptilt", "nmar_engine")
  return(config)
}

