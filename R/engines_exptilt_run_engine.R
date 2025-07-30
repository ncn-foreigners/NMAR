run_engine.nmar_engine_exptilt <- function(engine, formula, data) {
#TODO
  # settings_final <- validate_method_settings(method_name='exptilt',settings=settings)
  # x,col_y,cols_y_observed=c(),cols_delta=c(),prob_model_type='logit',y_dens='normal',tol_value=0.00001,min_iter=10,max_iter=100,optim_method='Newton'

  outcome_variable <- all.vars(formula$outcome)
  covariates_for_outcome <- all.vars(formula$covariates_outcome)
  covariates_for_missingness <- all.vars(formula$covariates_missingness)


  settings_final <- engine
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
  # model$theta=c(1,1,1)
  model <- .nmar_exptilt_run(model)



  return(list(theta = model$theta
              ,est_mean=estim_mean(model)
              ,loss_value=model$loss_value
  ))
}
