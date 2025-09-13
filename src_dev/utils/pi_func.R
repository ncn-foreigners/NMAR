#' Calculate probability or its derivative for logistic/probit models
#'
#' @description
#' Internal function used to compute either:
#' - predicted probabilities (for `func = "reg"`), or
#' - derivatives of the link function (for `func = "deriv"`).
#'
#' @param model A model object containing:
#'   - `theta` (numeric vector): Model coefficients.
#'   - `y_1` (numeric): Response variable value.
#'   - `prob_model_type` (character): Either `"logit"` or `"probit"`.
#' @param x Numeric matrix/vector of predictor variables.
#' @param func Character specifying computation type:
#'   - `"reg"`: Returns probabilities (default).
#'   - `"deriv"`: Returns derivatives of the link function.
#' @param theta Optional override of model coefficients (default: `model$theta`).
#'
#' @return
#' - If `func = "reg"`: Numeric vector of probabilities.
#' - If `func = "deriv"`: Matrix of derivatives (same nrow as `x`).
#'
#' @details
#' For internal use only. The function augments `x` with an intercept and `y_1` before computation.
#'
#' @examples
#' \dontrun{
#' # Example for a logistic model
#' model <- list(
#'   theta = c(0.5, -1.2, 0.3),
#'   y_1 = c(0,1,1),
#'   prob_model_type = "logit"
#' )
#' x <- matrix(c(0.1, 0.5, 0.8), ncol = 1)
#' pi_func(model, x, func = "reg")  # probabilities
#' pi_func(model, x, func = "deriv")  # derivatives
#' }
#' @keywords internal


pi_func <- function(model,x,func='reg',theta=model$theta) {
  # pi_func <- function(theta, x_for_delta, y, type = "logit", func = "reg") {
  # stopifnot(
  #   !any(is.na(y)),
  #   !any(is.na(theta)),
  #   !any(is.na(x_for_delta)),
  #   length(theta) == ncol(x_for_delta) + 2, # intercept + parametry dla x + parametr dla y
  #   type %in% c("logit", "probit"),
  #   func %in% c("reg", "deriv")
  # )

  x_mat <- as.matrix(x)
  x_aug <- cbind(1, x_mat, model$y_1)

  eta <- as.vector(x_aug %*% model$theta)

  if (model$prob_model_type == "logit") {
    if (func == "reg") {
      return(plogis(eta))
    } else {

      p <- plogis(eta)
      deriv_common <- p * (1 - p)
      return(sweep(x_aug, MARGIN = 1, STATS = deriv_common, FUN = "*"))
    }
  } else if (model$prob_model_type == "probit") {
    if (func == "reg") {
      return(pnorm(eta))
    } else {
      deriv_common <- dnorm(eta)
      return(sweep(x_aug, MARGIN = 1, STATS = deriv_common, FUN = "*"))
    }
  } else {
    stop("type must be 'logit' or 'probit'")
  }
}
