#' @exportS3Method run_engine nmar_engine_exptilt
run_engine.nmar_engine_exptilt <- function(engine, spec) {
  outcome_variable <- spec$outcome[[1]]
  covariates_for_outcome <- spec$auxiliary_vars
  covariates_for_missingness <- spec$response_predictors

  data_required <- unique(c(outcome_variable, covariates_for_outcome, covariates_for_missingness))

  if (spec$is_survey) {
    # Keep the original survey.design object intact, but trim the variables
    # so downstream code sees exactly the columns required by this engine
    data_object <- spec$original_data
    data_object$variables <- data_object$variables[, data_required, drop = FALSE]
  } else {
    data_object <- spec$data[, data_required, drop = FALSE]
  }

  model <- structure(
    list(
      data = if (spec$is_survey) spec$original_data else data_object,
      col_y = outcome_variable,
      cols_y_observed = covariates_for_outcome,
      cols_delta = covariates_for_missingness,
      prob_model_type =engine$prob_model_type,
      y_dens =engine$y_dens,
      tol_value =engine$tol_value,
      min_iter =engine$min_iter,
      auxiliary_means =engine$auxiliary_means,
      standardize =engine$standardize,
      max_iter =engine$max_iter,
      optim_method =engine$optim_method,
      variance_method = engine$variance_method,
      bootstrap_reps = engine$bootstrap_reps
    ),
    class = "nmar_exptilt"
  )

  if (spec$is_survey) {
    # Flag the model as survey-backed and retain the design so variance/diagnostics
    # can recover replicate structures when requested (e.g., bootstrap path)
    model$is_survey <- TRUE
    model$design <- data_object
  }

  model$family <- if (model$prob_model_type == "logit") {
    logit_family()
  } else if (model$prob_model_type == "probit") {
    probit_family()
  }

  # Keep a pristine copy of the pre-fit model so bootstrap replicates can reuse
  # the exact same starting values without inheriting stateful mutations
  model$original_params <- unserialize(serialize(model, NULL))
  model <- exptilt(data_object, model)
  if (!inherits(model, "nmar_result_exptilt")) {
    stop("Exptilt engine did not return an 'nmar_result_exptilt' object.")
  }

  return(model)
}

estim_mean <- function(model) {
  UseMethod("estim_mean", model)
}
estim_var <- function(model) {
  UseMethod("estim_var", model)
}
generate_Odds <- function(model,...) {
  UseMethod("generate_Odds", model)
}
s_function <- function(model, ...) {
  UseMethod("s_function", model)
}
