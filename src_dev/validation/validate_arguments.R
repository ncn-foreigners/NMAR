#' Custom Validation Module
#'
#' This module provides a set of internal helper functions for argument validation.
#' The functions are designed to throw a clear error message if a validation check fails.
#'
#' @keywords internal
validator <- new.env()

# --- Internal helper function for positive numbers ---
#' Checks if a value is a single positive number.
#'
#' @param x A numeric value.
#' @return A logical value.
#' @keywords internal
#' @name validator_is_positive
is_positive <- function(x) {
  is.numeric(x) && x > 0
}

# --- Custom assert_choice ---
#' @title Assert Choice
#' @description Checks if a value is one of the allowed choices.
#' @param x The value to check.
#' @param choices A vector of allowed choices.
#' @param name A string representing the name of the argument being checked.
#' @return Returns nothing on success, stops with an error on failure.
#' @name validator_assert_choice
validator$assert_choice <- function(x, choices, name) {
  if (!x %in% choices) {
    stop(
      paste0(
        "Argument '", name, "' should be one of (",
        paste(choices, collapse = ", "),
        "). Value is '", x, "'."
      )
    )
  }
}

# --- Custom assert_number (for min/max) ---
#' @title Assert Number
#' @description Checks if a value is a number within a specified range.
#' @param x The value to check.
#' @param name A string representing the name of the argument being checked.
#' @param min The minimum allowed value (inclusive). Defaults to -Inf.
#' @param max The maximum allowed value (inclusive). Defaults to Inf.
#' @return Returns nothing on success, stops with an error on failure.
#' @name validator_assert_number
validator$assert_number <- function(x, name, min = -Inf, max = Inf) {
  if (!is.numeric(x) || x < min || x > max) {
    stop(
      paste0(
        "Argument '", name, "' should be a number in range (",
        min, ", ", max, "). Value is '", x, "'."
      )
    )
  }
}

#' @title Assert Positive Integer
#' @description Checks if a value is a single, positive, and optionally finite integer.
#' @param x The value to check.
#' @param name A string representing the name of the argument being checked.
#' @param is.finite A logical value. If TRUE (default), checks if the number is finite.
#' @return Returns nothing on success. Stops with a clear error on failure.
#' @name validator_assert_positive_integer
validator$assert_positive_integer <- function(x, name, is.finite = TRUE) {
  if (!is.numeric(x) || x != as.integer(x)) {
    stop(paste0("Argument '", name, "' must be an integer. Value is ", x, "."))
  }

  if (x <= 0) {
    stop(paste0("Argument '", name, "' must be a positive integer. Value is ", x, "."))
  }

  if (is.finite && !base::is.finite(x)) {
    stop(paste0("Argument '", name, "' must be a finite integer. Value is ", x, "."))
  }
}

# --- Custom assert_logical ---
#' @title Assert Logical
#' @description Checks if a value is a single logical value.
#' @param x The value to check.
#' @param name A string representing the name of the argument being checked.
#' @return Returns nothing on success, stops with an error on failure.
#' @name validator_assert_logical
validator$assert_logical <- function(x, name) {
  if (!is.logical(x)) {
    stop(
      paste0(
        "Argument '", name, "' should be a logical value (TRUE/FALSE). Value is '", x, "'."
      )
    )
  }
}
