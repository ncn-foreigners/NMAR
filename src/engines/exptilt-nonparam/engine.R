exptilt_nonparam <- function(
    outcome_cols,
    refusal_col,
    # max_iter = get_json_param_info(all_schemas, "nonparametric_em", "max_iter")$default,
    # tol_value = get_json_param_info(all_schemas, "nonparametric_em", "tol_value")$default
    #TODO
    max_iter=100,
    tol_value=1e-6
) {
  json_path <- system.file("extdata", "method_params.json", package = "nmar")

  if (!file.exists(json_path)) {
    stop("Method parameters JSON file (method_params.json) not found in 'inst/extdata/'. Package might be corrupted or not installed correctly.")
  }

  all_schemas <- jsonlite::read_json(json_path, simplifyVector = TRUE)

  config <- list(
    outcome_cols = outcome_cols,
    refusal_col = refusal_col,
    max_iter = get_json_param_info(all_schemas, "exptilt_nonparam", "max_iter")$default,
    tol_value = get_json_param_info(all_schemas, "exptilt_nonparam", "tol_value")$default
  )

  # TODO: add validation
  # config <- validate_method_settings(...)

  class(config) <- c("nmar_engine_exptilt_nonparam", "nmar_engine")
  return(config)
}
