#' @exportS3Method run_engine nmar_engine_exptilt
run_engine.nmar_engine_exptilt <- function(engine, task) {
# All ET fits flow through the shared design prep so data-frame and survey
# paths inherit the same scaling, auxiliary moments, and weight handling as EL
  design_info <- prepare_nmar_design(
    task,
    standardize = engine$standardize,
    auxiliary_means = engine$auxiliary_means,
    include_response = TRUE,
    include_auxiliary = TRUE
  )

# Reconstruct a formula carrying response-only predictors to the right of `|`
  f_use <- nmar_rebuild_partitioned_formula(
    base_formula = task$formula,
    response_predictors = design_info$response_predictors,
    env = task$environment
  )

  args <- list(
    data = design_info$survey_design %||% design_info$data,
    formula = f_use,
    auxiliary_means = design_info$auxiliary_means,
    standardize = design_info$standardize,
    prob_model_type = engine$prob_model_type,
    y_dens = engine$y_dens,
    variance_method = engine$variance_method,
    bootstrap_reps = engine$bootstrap_reps,
    control = engine$control,
    stopping_threshold = engine$stopping_threshold,
    on_failure = engine$on_failure,
    supress_warnings = engine$supress_warnings,
    verbose = engine$verbose,
    trace_level = engine$trace_level
  )

  if (!isTRUE(design_info$is_survey)) {
# For data.frames we pass explicit design weights (all ones unless supplied).
# For survey designs the method extracts weights via stats::weights(), so we
# avoid passing a duplicate argument to prevent ambiguity
    args$design_weights <- design_info$weights
  }

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
generate_Odds <- function(model, ...) {
  UseMethod("generate_Odds", model)
}
s_function <- function(model, ...) {
  UseMethod("s_function", model)
}
