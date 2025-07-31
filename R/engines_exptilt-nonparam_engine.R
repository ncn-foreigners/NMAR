exptilt_nonparam <- function(
    outcome_cols,
    refusal_col,
    # max_iter = get_json_param_info(all_schemas, "nonparametric_em", "max_iter")$default,
    # tol_value = get_json_param_info(all_schemas, "nonparametric_em", "tol_value")$default
    #TODO
    max_iter=100,
    tol_value=1e-6
) {
  config <- list(
    outcome_cols = outcome_cols,
    refusal_col = refusal_col,
    max_iter = max_iter,
    tol_value = tol_value
  )

  # TODO: add validation
  # config <- validate_method_settings(...)

  class(config) <- c("nmar_engine_exptilt_nonparam", "nmar_engine")
  return(config)
}
