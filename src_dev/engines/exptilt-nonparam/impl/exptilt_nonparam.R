#' Nonparametric Exponential Tilting (Internal Generic)
#'
#' @param data A data.frame or survey.design object
#' @param ... Other arguments passed to methods
#' @keywords internal
exptilt_nonparam <- function(data, ...) {
  UseMethod("exptilt_nonparam")
}
