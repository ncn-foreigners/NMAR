#' Exponential Tilting Engine for NMAR Estimation
#'
#' @param prob_model_type Probability model for data missingness ("logit" or "probit")
#' @param y_dens Density model for outcome variable ("normal" or "gamma")
#' @param tol_value Tolerance in convergence (default = 0.00001)
#' @param min_iter Minimum number of iterations (default = 10)
#' @param max_iter Maximum number of iterations (default = 100)
#' @param optim_method Optimization method ("Newton" or "Broyden")
#' @return Engine configuration object of class `nmar_engine_exptilt`
#' @export
exptilt <- function(
    prob_model_type = "logit",
    y_dens = "normal",
    tol_value = 0.00001,
    min_iter = 10,
    max_iter = 100,
    optim_method = "Newton") {

  # Walidacja parametrów
  if (!prob_model_type %in% c("logit", "probit")) {
    stop("prob_model_type must be 'logit' or 'probit'")
  }
  if (!y_dens %in% c("normal", "gamma")) {
    stop("y_dens must be 'normal' or 'gamma'")
  }
  if (!optim_method %in% c("Newton", "Broyden")) {
    stop("optim_method must be 'Newton' or 'Broyden'")
  }

  # Tworzenie obiektu konfiguracyjnego
  config <- list(
    prob_model_type = prob_model_type,
    y_dens = y_dens,
    tol_value = tol_value,
    min_iter = min_iter,
    max_iter = max_iter,
    optim_method = optim_method
  )

  class(config) <- c("nmar_engine_exptilt", "nmar_engine")
  return(config)
}




#
#
#
# # nmar_exptilt <- function(){
# #     cat('Hello exptilt')
# #   }
# nmar_exptilt_engine <- function(data, outcome_variable, covariates_for_outcome, covariates_for_missingness, settings) {
# #TODO
#   # settings_final <- validate_method_settings(method_name='exptilt',settings=settings)
#   # x,col_y,cols_y_observed=c(),cols_delta=c(),prob_model_type='logit',y_dens='normal',tol_value=0.00001,min_iter=10,max_iter=100,optim_method='Newton'
#
#   settings_final <- settings
#   model <- .nmar_exptilt_create(
#     x = data,
#     col_y = outcome_variable,
#     cols_y_observed = covariates_for_outcome,
#     cols_delta = covariates_for_missingness,
#     prob_model_type = settings_final$prob_model_type,
#     y_dens = settings_final$y_dens,
#     tol_value = settings_final$tol_value,
#     min_iter = settings_final$min_iter,
#     max_iter = settings_final$max_iter,
#     optim_method = settings_final$optim_method
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
