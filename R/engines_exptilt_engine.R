


nmar_exptilt <- function(){
    cat('Hello exptilt')
  }
# nmar_exptilt() <- function(x,col_y,cols_y_observed=c(),cols_delta=c(),prob_model_type='logit',y_dens='normal',tol_value=0.00001,min_iter=10,max_iter=100,optim_method='Newton') {
#   cat('=========')
#   cat(cols_y_observed)
#   # stopifnot(
#   #   "At least one observed variable must be specified in 'cols_y_observed'" =
#   #     length(cols_y_observed) > 0,
#   #
#   #   "'prob_model_type' must be either 'logit' or 'probit'" =
#   #     prob_model_type %in% c('logit','probit'),
#   #
#   #   "'y_dens' must be either 'normal' or 'gamma'" =
#   #     y_dens %in% c('normal','gamma'),
#   #
#   #   "'min_iter' must be less than 'max_iter'" =
#   #     min_iter < max_iter,
#   #
#   #   "min_iter must be greater than 0" =
#   #     min_iter > 0,
#   #
#   #   "'optim_method' must be either 'Newton' or 'Broyden'" =
#   #     optim_method %in% c('Newton','Broyden'),
#   #
#   #   "Exactly one outcome variable must be specified in 'col_y'" =
#   #     length(col_y) == 1,
#   #
#   #   "Column Y must be present in the data frame" =
#   #     col_y %in% colnames(x),
#   #
#   #   "Columns in 'cols_y_observed' must be present in the data frame" =
#   #     all(cols_y_observed %in% colnames(x)),
#   #
#   #   "Column Y must not be in 'cols_y_observed'" =
#   #     !(col_y %in% cols_y_observed),
#   # )
#   model <- .nmar_exptilt_create(
#     x = x,
#     col_y = col_y,
#     cols_y_observed = cols_y_observed,
#     cols_delta = cols_delta,
#     prob_model_type = prob_model_type,
#     y_dens = y_dens,
#     tol_value = tol_value,
#     min_iter = min_iter,
#     max_iter = max_iter,
#     optim_method = optim_method
#   )
#   # model$theta=c(1,1,1)
#   model <- .nmar_exptilt_run(model)
#
#
#
#   return(list(theta = model$theta
#               ,est_mean=estim_mean(model)
#               ,loss_value=model$loss_value
#   ))
# }
