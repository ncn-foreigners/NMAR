#' Empirical likelihood estimator
#' @description Generic for the empirical likelihood (EL) estimator under NMAR.
#'   Methods are provided for \code{data.frame} and \code{survey.design}.
#' @param data A \code{data.frame} or a \code{survey.design}.
#' @param ... Passed to class-specific methods.
#' @seealso \code{\link{el_engine}}
#' @keywords internal
el <- function(data, ...) {
  UseMethod("el")
}
