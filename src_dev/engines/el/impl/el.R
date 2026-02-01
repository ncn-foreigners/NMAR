#' Empirical likelihood estimator
#'
#' @param data A \code{data.frame} or a \code{survey.design}.
#' @param ... Passed to class-specific methods.
#'
#' @keywords internal
el <- function(data, ...) {
  UseMethod("el")
}
