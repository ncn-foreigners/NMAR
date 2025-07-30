#' @importFrom stats coef runif pnorm dnorm dgamma plogis as.formula coef sd setNames

.nmar_exptilt_create <- function(x, col_y, cols_y_observed, cols_delta, prob_model_type, y_dens, tol_value, min_iter, max_iter, optim_method) {
  # stopifnot(length(cols_y_observed) > 0,
  #           length(cols_delta) >= 0,
  #           is.matrix(x) || is.data.frame(x),
  #           length(col_y) == 1,
  #           prob_model_type %in% c("logit", "probit"),
  #           y_dens %in% c("normal", "gamma"),
  #           tol_value > 0,
  #           min_iter >= 0,
  #           max_iter >= min_iter,
  #           optim_method %in% c("Newton", "BFGS"),
  #           is.numeric(x[, col_y]),
  #           # is.numeric(x[, cols_y_observed]),
  #           all(cols_delta %in% colnames(x)),
  #           all(col_y %in% colnames(x)),
  #           all(cols_y_observed %in% colnames(x)),
  #           any(is.na(x[,col_y]))
#
  # )
  structure(
    list(
      x = x,
      col_y = col_y,
      cols_y_observed = cols_y_observed,
      cols_delta = cols_delta,
      prob_model_type = prob_model_type,
      y_dens = y_dens,
      tol_value = tol_value,
      min_iter = min_iter,
      max_iter = max_iter,
      optim_method = optim_method
    ),
    class = "nmar_exptilt"
  )

}
