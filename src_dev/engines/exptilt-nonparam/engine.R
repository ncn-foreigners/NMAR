#' @title Nonparametric Exponential Tilting Engine for NMAR Estimation
#'
#' @description Creates a configuration object for the nonparametric Exponential Tilting (ET) method
#'   used in Not Missing at Random (NMAR) estimation. This version uses a nonparametric
#'   approach to estimate the missing data mechanism without distributional assumptions.
#'   The returned object can be passed to the `nmar()` function for NMAR estimation.
#'
#' @param outcome_cols A character vector specifying the names of columns containing
#'   the outcome variables in the data.
#' @param refusal_col A character string specifying the name of the column containing
#'   refusal/missingness counts.
#' @param max_iter An integer specifying the maximum number of iterations for the
#'   EM algorithm. Defaults to 100.
#' @param tol_value A numeric value representing the convergence tolerance for the
#'   EM algorithm. The algorithm stops when parameter changes are below this value.
#'   Defaults to 1e-6.
#'
#' @return An engine configuration object of S3 class `c("nmar_engine_exptilt_nonparam", "nmar_engine")`.
#'   This object contains all specified parameters for the nonparametric ET method.
#'
#' @importFrom jsonlite read_json

#'
#' @examples
#' # Create a nonparametric ET configuration with default parameters
#' config <- exptilt_nonparam(
#'   outcome_cols = c("vote_A", "vote_B", "abstain"),
#'   refusal_col = "refused"
#' )
#'
#' # Create with custom iteration parameters
#' custom_config <- exptilt_nonparam(
#'   outcome_cols = c("yes", "no"),
#'   refusal_col = "missing",
#'   max_iter = 200,
#'   tol_value = 1e-8
#' )
#'
#' # Use with nmar() function:
#' # results <- nmar(data = survey_data, engine = config)
#' @export
exptilt_nonparam_engine <- function(
    # outcome_cols,
    refusal_col,
    # max_iter = get_json_param_info(all_schemas, "nonparametric_em", "max_iter")$default,
    # tol_value = get_json_param_info(all_schemas, "nonparametric_em", "tol_value")$default
    #TODO
    max_iter=100,
    tol_value=1e-6
) {
  # json_path <- system.file("extdata", "method_params.json", package = "nmar")

  # if (!file.exists(json_path)) {
  #   stop("Method parameters JSON file (method_params.json) not found in 'inst/extdata/'. Package might be corrupted or not installed correctly.")
  # }

  # all_schemas <- jsonlite::read_json(json_path, simplifyVector = TRUE)

  config <- list(
    # outcome_cols = outcome_cols,
    refusal_col = refusal_col,
    max_iter = 100, #hotfix
    tol_value = 1e-6 #hotfix
  )

  # TODO: add validation
  # config <- validate_method_settings(...)

  class(config) <- c("nmar_engine_exptilt_nonparam", "nmar_engine")
  return(config)
}
