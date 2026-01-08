#' Nonparametric Exponential Tilting (Internal Generic)
#'
#' @param data A data.frame or survey.design object
#' @param ... Other arguments passed to methods
#' @return An engine-specific NMAR result object for the nonparametric
#'   exponential tilting estimator.
#' @keywords internal
exptilt_nonparam <- function(data, ...) {
  UseMethod("exptilt_nonparam")
}
