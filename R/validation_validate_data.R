#' @title Validate Data for NMAR Analysis
#' @description Little sanity-check for data
#'
#' @param data A data frame or a survey object.

#' @return Returns `invisible(NULL)` on success, stopping with a descriptive error on failure.
#' @keywords internal
validate_data <- function(data) {
# Validate data object type
  if (!inherits(data, c("data.frame", "survey.design"))) {
    stop("'data' must be a data.frame or survey.design object. Received: ", class(data)[1], call. = FALSE)
  }

# Extract variables from survey object
  if (inherits(data, "survey.design")) {
    data <- data$variables
  }

# Check for empty data
  if (nrow(data) == 0) {
    stop("Input dataset is empty (0 rows).", call. = FALSE)
  }


  invisible(NULL)
}
