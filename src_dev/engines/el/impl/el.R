#' Empirical likelihood estimator
#' @description Generic for the empirical likelihood (EL) estimator under NMAR.
#'   Methods are provided for `data.frame` and `survey.design`.
#' @param data A `data.frame` or a `survey.design`.
#' @param ... Passed to class-specific methods.
#' @seealso `el.data.frame()`, `el.survey.design()`, `el_engine()`
#' @keywords internal
el <- function(data, ...) {
  UseMethod("el")
}
