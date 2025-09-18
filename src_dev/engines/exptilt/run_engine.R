#' @exportS3Method NULL
run_engine.nmar_engine_exptilt <- function(engine, formula, data,response_predictors) {
  outcome_variable <- as.vector(all.vars(formula[[2]]))
  covariates_for_outcome <- as.vector(all.vars(formula[[3]]))
  covariates_for_missingness <- response_predictors

   validate_data(data, outcome_variable, covariates_for_outcome, covariates_for_missingness)

  data_ <- data[,c(outcome_variable,covariates_for_outcome,covariates_for_missingness)]
  model <- structure(
    list(
      data = data_,
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

  model$family <- if (model$prob_model_type == "logit") {
    logit_family()
  } else if (model$prob_model_type == "probit") {
    probit_family()
  }

  model$original_params <- model
  model <- exptilt(data_,model)
  if (!inherits(model, "nmar_result_exptilt")) {
    stop("Exptilt engine did not return an 'nmar_result_el' object.")
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


