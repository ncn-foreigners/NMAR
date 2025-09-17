#' @importFrom nleqslv nleqslv
#' @importFrom stats as.formula coef dnorm dgamma sd setNames
#' @export
#' #todo on_failure logic
exptilt.data.frame <- function(data,model,on_failure=c('return')){
  model$x=data
  # model$x_1 <- model$x[!is.na(model$x[,model$col_y]),,drop=FALSE] #observed
  # model$x_0 <- model$x[is.na(model$x[,model$col_y]),,drop=FALSE] #unobserved
  # model$y_1 <- model$x_1[,model$col_y,drop=TRUE] #observed y
  #
  # model$x_for_y_obs <- model$x_1[,model$cols_y_observed,drop=FALSE]
  # model$x_for_y_unobs <- model$x_0[,model$cols_y_observed,drop=FALSE]

  has_aux = length(model$coly) > 0#TODO make sure not col_y
  mu_x_unscaled <- if (has_aux) auxiliary_means else NULL

  scaling_result <- validate_and_apply_nmar_scaling(
    standardize = model$standardize,
    has_aux = has_aux,#TODO make sure not col_y
    response_model_matrix_unscaled = model$x[,c(model$col_y,model$cols_delta),drop=FALSE],
    auxiliary_matrix_unscaled = model$x[,model$cols_y_observed,drop=FALSE],
    mu_x_unscaled = model$auxiliary_means
  )
  model$nmar_scaling_recipe <- scaling_result$nmar_scaling_recipe
  response_model_matrix_scaled <- scaling_result$response_model_matrix_scaled
  auxiliary_matrix_scaled <- scaling_result$auxiliary_matrix_scaled
  mu_x_scaled <- scaling_result$mu_x_scaled

  model$x_1 <- response_model_matrix_scaled[!is.na(response_model_matrix_scaled[,model$col_y]),,drop=FALSE] #observed
  model$x_0 <- response_model_matrix_scaled[is.na(response_model_matrix_scaled[,model$col_y]),,drop=FALSE] #unobserved
  model$y_1 <- model$x_1[,model$col_y,drop=TRUE] #observed y
  model$x_for_y_obs <- auxiliary_matrix_scaled[!is.na(response_model_matrix_scaled[,model$col_y]),,drop=FALSE]
  model$x_for_y_unobs <- auxiliary_matrix_scaled[is.na(response_model_matrix_scaled[,model$col_y]),,drop=FALSE]

  # model$respondent_weights <- weights(model$x_1)

  model$respondent_weights <- rep(1, nrow(model$x_1))
  #
  stopifnot(
    nrow(model$x_0)>0,
    nrow(model$x_1)>0
  )

  model$theta=stats::runif(length(model$cols_delta)+2,0,0.1)

  dens_response <- generate_conditional_density(model)
  # #peek gradients func
  model$density_fun <- dens_response$density_function
  model$density_fun_gradient <- dens_response$density_function_grad
  model$density_fun_hess <- dens_response$density_function_hess
  model$density_num_of_coefs <- dens_response$num_of_coefs
  model$O_matrix_nieobs <- generate_Odds(model)
  #const
  model$f_matrix_nieobs <- generate_conditional_density_matrix(model)
  model$C_matrix_nieobs <- generate_C_matrix(model)

  # cat(model)



  target_function <- function(theta) {
    model$theta <<- theta  # Natychmiastowa aktualizacja
    step_func(model,theta, model$O_matrix_nieobs)
    # step_func_2(theta, O_matrix_nieobs, f_matrix_nieobs, C_matrix_nieobs)
  }


  solution <- nleqslv(
    x = model$theta,
    fn = target_function,
    method = "Newton",
    control = list(
      maxit = 1

    )
  )

  theta_prev=model$theta
  model$theta <- solution$x
  model$loss_value <- solution$fvec
  iter=0

  while(sum((model$theta-theta_prev)^2) > model$tol_value || (iter < model$min_iter && iter < model$max_iter)) {

    solution <- nleqslv(
      x = model$theta,
      fn = target_function,
      method = "Newton",
      control = list(
        maxit = 1

      )

    )
    # cat("Iter:", iter, "Theta:", model$theta, "Change:", sum((model$theta-theta_prev)^2), "\n")

    theta_prev <- model$theta
    model$theta <- solution$x
    model$loss_value <- solution$fvec
    # O_matrix_nieobs <- generate_Odds(theta_prev,x_0[,cols_delta],y_1)
    #
    model$O_matrix_nieobs <- generate_Odds(model)
    iter<-iter+1



  }
  model$theta <- if (model$standardize) unscale_coefficients(model$theta, estim_var(model)$vcov, model$nmar_scaling_recipe)$coefficients else model$theta
  # names(beta_hat_unscaled) <- colnames(response_model_matrix_unscaled)
  model$x_1 <- model$x[!is.na(model$x[,model$col_y]),,drop=FALSE] #observed
  model$x_0 <- model$x[is.na(model$x[,model$col_y]),,drop=FALSE] #unobserved
  model$y_1 <- model$x_1[,model$col_y,drop=TRUE] #observed y

  model$x_for_y_obs <- model$x_1[,model$cols_y_observed,drop=FALSE]
  model$x_for_y_unobs <- model$x_0[,model$cols_y_observed,drop=FALSE]


    se_final=NaN
  if(model$variance_method=="bootstrap"){

    model$original_params$variance_method='delta'
    cat(model$bootstrap_reps)
    test_var <- as.list(model$original_params)
    bootstrap_results <- do.call(
      bootstrap_variance,
      list(data=model$x,estimator=exptilt,point_estimate=estim_mean(model), bootstrap_reps=model$bootstrap_reps,as.list(model$original_params)
           ),
      # as.list(model$original_params)
    )
    se_final <- bootstrap_results$se

  }
  else{
    se_final <- sqrt(estim_var(model)$var_est)
  }


  return(validate_nmar_result(new_nmar_result_exptilt(
    y_hat=estim_mean(model),
    se_final = se_final,
    weights=NULL, #TODO
    coefficients = model$theta,
    vcov = estim_var(model)$vcov,
    converged = TRUE, #TODO
    class="nmar_result_exptilt"

  ), "nmar_result_exptilt"))
  # return(model)
}


