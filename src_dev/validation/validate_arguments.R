#' Custom Validation Module
#'
#' This module provides a set of internal helper functions for argument validation.
#' The functions are designed to throw a clear error message if a validation check fails.
#'
#' @keywords internal
validator <- new.env()

#' Validation helpers for NMAR (internal)
#' @name nmar_validator_helpers
#' @keywords internal
NULL

#' Checks if a value is a single positive number.
#'
#' @param x A numeric value.
#' @return A logical value.
#' @keywords internal
#' @rdname nmar_validator_helpers
#' @noRd
is_positive <- function(x) {
  is.numeric(x) && x > 0
}

#' @title Assert Choice
#' @description Checks if a value is one of the allowed choices.
#' @param x The value to check.
#' @param choices A vector of allowed choices.
#' @param name A string representing the name of the argument being checked.
#' @return Returns nothing on success, stops with an error on failure.
#' @rdname nmar_validator_helpers
#' @noRd
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

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_list <- function(x, name) {
  if (!is.list(x)) {
    stop(paste0("Argument '", name, "' must be a list."), call. = FALSE)
  }
}

#' @title Assert Number
#' @description Checks if a value is a number within a specified range.
#' @param x The value to check.
#' @param name A string representing the name of the argument being checked.
#' @param min The minimum allowed value (inclusive). Defaults to -Inf.
#' @param max The maximum allowed value (inclusive). Defaults to Inf.
#' @return Returns nothing on success, stops with an error on failure.
#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_number <- function(x, name, min = -Inf, max = Inf) {
  if (!is.numeric(x) || any(x < min) || any(x > max)) {
    stop(
      paste0(
        "Argument '", name, "' should be a number in range (",
        min, ", ", max, "). Value is '", x, "'."
      )
    )
  }
}

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_positive_number <- function(x, name, allow_infinite = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(paste0("Argument '", name, "' must be a single numeric value."), call. = FALSE)
  }
  if (!allow_infinite && !is.finite(x)) {
    stop(paste0("Argument '", name, "' must be finite."), call. = FALSE)
  }
  if (!(is.finite(x) && x > 0) && !(allow_infinite && is.infinite(x) && x > 0)) {
    stop(paste0("Argument '", name, "' must be positive", if (allow_infinite) " or Inf" else "", "."), call. = FALSE)
  }
}

#' @title Assert Positive Integer
#' @description Checks if a value is a single, positive, and optionally finite integer.
#' @param x The value to check.
#' @param name A string representing the name of the argument being checked.
#' @param is.finite A logical value. If TRUE (default), checks if the number is finite.
#' @return Returns nothing on success. Stops with a clear error on failure.
#' @rdname nmar_validator_helpers
#' @noRd
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

#' @title Assert Logical
#' @description Checks if a value is a single logical value.
#' @param x The value to check.
#' @param name A string representing the name of the argument being checked.
#' @return Returns nothing on success, stops with an error on failure.
#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_logical <- function(x, name) {
  if (!is.logical(x)) {
    stop(
      paste0(
        "Argument '", name, "' should be a logical value (TRUE/FALSE). Value is '", x, "'."
      )
    )
  }
}

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_scalar_numeric <- function(x, name, allow_na = FALSE, finite = TRUE) {
  ok_type <- is.numeric(x) && length(x) == 1L
  ok_na <- allow_na || (!is.na(x))
  ok_fin <- (!finite) || is.finite(x) || (allow_na && is.na(x))
  if (!(ok_type && ok_na && ok_fin)) {
    stop(paste0("Argument '", name, "' must be a ", if (finite) "finite " else "", "numeric scalar",
                if (allow_na) " (NA allowed)" else "", "."), call. = FALSE)
  }
}

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_scalar_character <- function(x, name, allow_na = FALSE, non_empty = TRUE) {
  ok_type <- is.character(x) && length(x) == 1L
  ok_na <- allow_na || (!is.na(x))
  ok_ne <- (!non_empty) || (nzchar(x))
  if (!(ok_type && ok_na && ok_ne)) {
    stop(paste0("Argument '", name, "' must be a character scalar",
                if (allow_na) " (NA allowed)" else "",
                if (non_empty) " and non-empty" else "", "."), call. = FALSE)
  }
}

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_scalar_logical <- function(x, name) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    stop(paste0("Argument '", name, "' must be a logical scalar (TRUE/FALSE)."), call. = FALSE)
  }
}

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_matrix_ncol <- function(X, expected, name) {
  if (!is.matrix(X)) {
    stop(paste0("Argument '", name, "' must be a matrix."), call. = FALSE)
  }
  if (ncol(X) != expected) {
    stop(paste0("Matrix '", name, "' has ", ncol(X), " columns but expected ", expected, "."), call. = FALSE)
  }
}

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_formula_two_sided <- function(fml, name) {
  if (!inherits(fml, "formula") || length(fml) != 3L) {
    stop(paste0("'", name, "' must be a two-sided formula."), call. = FALSE)
  }
}

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_data_frame_or_survey <- function(x, name) {
  if (!inherits(x, c("data.frame", "survey.design"))) {
    stop(paste0("'", name, "' must be a data.frame or survey.design object."), call. = FALSE)
  }
}

#' @rdname nmar_validator_helpers
#' @noRd
validator$assert_named_numeric <- function(x, name, allow_null = TRUE) {
  if (is.null(x)) {
    if (allow_null) return(invisible(TRUE))
    stop(paste0("Argument '", name, "' must be a named numeric vector (not NULL)."), call. = FALSE)
  }
  if (!is.numeric(x) || is.null(names(x)) || anyNA(names(x))) {
    stop(paste0("Argument '", name, "' must be a named numeric vector with non-NA names."), call. = FALSE)
  }
  invisible(TRUE)
}
