#' @keywords internal
validate_nmar_result<- function(x,class_name) {
  browser()
  stopifnot(is.list(x), inherits(x, class_name))
  if (isTRUE(x$converged)) {
    stopifnot(is.finite(x$y_hat), is.numeric(x$se))
  }
  x
}

# validate_nmar_result <- function(x, ...) {
#   UseMethod("validate_nmar_result")
# }
