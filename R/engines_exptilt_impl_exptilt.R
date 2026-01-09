#' Exponential tilting estimator
#' @description Generic for the exponential tilting (ET) estimator under NMAR.
#'   Methods are provided for `data.frame` and `survey.design`.
#' @param data A `data.frame` or a `survey.design`.
#' @param ... Passed to class-specific methods.
#' @return An engine-specific NMAR result object (for example
#'   \code{nmar_result_exptilt}).
#' @seealso `exptilt.data.frame()`, `exptilt.survey.design()`, `exptilt_engine()`
#' @keywords internal
exptilt <- function(data, ...) {
  UseMethod("exptilt")
}
