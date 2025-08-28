#' @importFrom nleqslv nleqslv
#' @importFrom stats as.formula coef dnorm dgamma sd setNames # Ensure these are imported via roxygen2
run_nmar_exptilt <- function(model){
  model$x_1 <- model$x[!is.na(model$x[,model$col_y]),,drop=FALSE] #observed
  model$x_0 <- model$x[is.na(model$x[,model$col_y]),,drop=FALSE] #unobserved
  model$y_1 <- model$x_1[,model$col_y,drop=TRUE] #observed y

  model$x_for_y_obs <- model$x_1[,model$cols_y_observed,drop=FALSE]
  model$x_for_y_unobs <- model$x_0[,model$cols_y_observed,drop=FALSE]
  # browser()
  stopifnot(
    nrow(model$x_0)>0,
    nrow(model$x_1)>0
  )

  model$theta=stats::runif(length(model$cols_delta)+2,0,0.1)

  dens_response <- generate_conditional_density(model)
  # browser()#peek gradients func
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
    # browser()
    model$O_matrix_nieobs <- generate_Odds(model)
    iter<-iter+1



  }
  return(model)
}


