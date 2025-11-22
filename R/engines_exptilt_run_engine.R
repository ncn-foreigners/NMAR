#' @exportS3Method run_engine nmar_engine_exptilt
run_engine.nmar_engine_exptilt <- function(engine, formula, data, trace_level) {

  args <- list(
    data = data,
    formula = formula,
    trace_level = trace_level,
    standardize = engine$standardize,
    prob_model_type = engine$prob_model_type,
    y_dens = engine$y_dens,
    variance_method = engine$variance_method,
    bootstrap_reps = engine$bootstrap_reps,
    control = engine$control,
    stopping_threshold = engine$stopping_threshold,
    on_failure = engine$on_failure,
    supress_warnings = engine$supress_warnings,
    sample_size = engine$sample_size
  )

  res <- do.call(exptilt, args)
  if (!inherits(res, "nmar_result_exptilt")) {
    stop("Exptilt engine did not return an 'nmar_result_exptilt' object.")
  }
  res
}

estim_mean <- function(model) {
  UseMethod("estim_mean", model)
}
estim_var <- function(model) {
  UseMethod("estim_var", model)
}
generate_Odds <- function(model, theta) {
  UseMethod("generate_Odds", model)
}
s_function <- function(model, delta, x, theta) {
  UseMethod("s_function", model)
}
# validate_df <- function(model, covariate_outcome, covariates_aux, covariates_missingness,X,Y,Z) {
#   # Dispatch on `model` class only. Additional args are for method signatures.
#   UseMethod("validate_df", model)
# }
