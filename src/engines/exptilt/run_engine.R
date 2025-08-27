run_engine.nmar_engine_exptilt <- function(engine, formula, data) {
  outcome_variable <- all.vars(formula$outcome)
  covariates_for_outcome <- all.vars(formula$covariates_outcome)
  covariates_for_missingness <- all.vars(formula$covariates_missingness)

  model <- structure(
    list(
      x = data,
      col_y = outcome_variable,
      cols_y_observed = covariates_for_outcome,
      cols_delta = covariates_for_missingness,
      prob_model_type =engine$prob_model_type,
      y_dens =engine$y_dens,
      tol_value =engine$tol_value,
      min_iter =engine$min_iter,
      max_iter =engine$max_iter,
      optim_method =engine$optim_method
    ),
    class = "nmar_exptilt"
  )

  model <- run_nmar_exptilt(model)

  return(structure(list(theta = model$theta
              ,est_mean=estim_mean(model)
              ,loss_value=model$loss_value
  ),class = "nmar_result"))
}
# run <-function(model) {
#   UseMethod("run")
# }

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

# step_func <- function(model,..){
#   UseMethod("step_func", model)
# }
