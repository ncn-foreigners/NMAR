#' Null-coalescing operator
#'
#' Returns the left-hand side if it is not NULL, otherwise returns the
#' right-hand side. This is an internal utility used throughout the NMAR
#' package for providing default values.
#'
#' @param a First value
#' @param b Default value if `a` is NULL
#' @return `a` if not NULL, otherwise `b`
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' NULL %||% "default" # Returns "default"
#' "value" %||% "default" # Returns "value"
#' NA %||% "default" # Returns NA (not NULL, so returns a)
#' }
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}
