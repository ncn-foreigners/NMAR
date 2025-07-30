estim_mean.nmar_exptilt <- function(model){
  numerator<-sum(model$y_1/pi_func(model, model$x_1[,model$cols_delta], model$x_1[,model$col_y], func = "reg"))
  denominator<-sum(1/pi_func(model, model$x_1[,model$cols_delta], model$x_1[,model$col_y], func = "reg"))

  return(sum(numerator/denominator))
}

