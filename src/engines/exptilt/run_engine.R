run_engine.nmar_engine_exptilt <- function(engine, formula, data) {
  outcome_variable <- all.vars(formula$outcome)
  covariates_for_outcome <- all.vars(formula$covariates_outcome)
  covariates_for_missingness <- all.vars(formula$covariates_missingness)


  settings_final <- engine
  #todo - unnecessary
  model <- .nmar_exptilt_create(
    x = data,
    col_y = outcome_variable,
    cols_y_observed = covariates_for_outcome,
    cols_delta = covariates_for_missingness,
    prob_model_type = settings_final$prob_model_type,
    y_dens = settings_final$y_dens,
    tol_value = settings_final$tol_value,
    min_iter = settings_final$min_iter,
    max_iter = settings_final$max_iter,
    optim_method = settings_final$optim_method
  )

  model <- .nmar_exptilt_run(model)

  return(list(theta = model$theta
              ,est_mean=estim_mean(model)
              ,loss_value=model$loss_value
  ))
}
