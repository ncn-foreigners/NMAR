# validate_nmar_result <- function(res) {
#   required_components <- c("theta", "est_mean", "loss_value")
#
#   if (!all(required_components %in% names(res))) {
#     missing_comps <- setdiff(required_components, names(res))
#     stop(paste("Result missing required components:",
#                paste(missing_comps, collapse = ", ")))
#   }
#
#   if (!is.numeric(res$theta)) stop("theta must be numeric")
#   if (!is.numeric(res$est_mean) || length(res$est_mean) != 1) {
#     stop("est_mean must be a single numeric value")
#   }
#   # Modified to handle vector loss_value
#   if (!is.numeric(res$loss_value)) {
#     stop("loss_value must be numeric")
#   }
#
#   invisible(TRUE)
# }
