#' Exponential tilting estimator
#' @description Generic for the exponential tilting (ET) estimator under NMAR.
#'   Methods are provided for `data.frame` and `survey.design`.
#' @param data A `data.frame` or a `survey.design`.
#' @param ... Passed to class-specific methods.
#' @seealso `exptilt.data.frame()`, `exptilt.survey.design()`, `exptilt_engine()`
#' @export
exptilt <- function(data, ...) {
  UseMethod("exptilt")
}
