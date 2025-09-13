#' Error/condition helpers
#'
#' Developer note:
#' - Add lightweight condition constructors here for common error types to keep
#'   signaling consistent across engines and shared utilities.
#' - Keep these unexported; they are meant for internal flow control and clear
#'   diagnostics (e.g., `convergenceError`).
#' - Prefer simple S3 condition classes like `c("<type>", "error", "condition")`.
#'
#' @keywords internal
convergenceError <- function(message, call = NULL) {
  structure(list(message = message, call = call), class = c("convergenceError", "error", "condition"))
}
