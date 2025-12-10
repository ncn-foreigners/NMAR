#' S3 helpers for NMAR engine objects
#'
#' Lightweight, user-facing methods for engine configuration objects
#' (class `nmar_engine`). These improve discoverability and provide a
#' consistent print surface across engines while keeping the objects as
#' simple lists internally.
#'
#' @section Design:
#' - `engine_name()` returns a canonical identifier used across the
#'   package (e.g., in `nmar_result$meta$engine_name`).
#' - `print.nmar_engine()` provides a concise, readable summary of the
#'   engine configuration; engine-specific classes reuse the parent
#'   method unless they need to override it.
#' - `engine_config()` returns the underlying configuration as a named
#'   list for programmatic inspection.
#'
#' @name nmar_engine_helpers
#' @keywords engine_view
NULL

#' Canonical engine name
#'
#' Returns a stable, machine-friendly identifier for an engine object. This
#' identifier is also used in `nmar_result$meta$engine_name` to keep a
#' consistent naming scheme between configurations and results.
#'
#' @param x An object inheriting from class `nmar_engine`.
#' @return A single character string, e.g. "empirical_likelihood".
#' @keywords engine_view
#' @export
engine_name <- function(x) {
  UseMethod("engine_name")
}

#' @keywords engine_view
#' @export
engine_name.nmar_engine <- function(x) "nmar_engine"

#' @keywords engine_view
#' @export
engine_name.nmar_engine_el <- function(x) "empirical_likelihood"

#' @keywords engine_view
#' @export
engine_name.nmar_engine_exptilt <- function(x) "exponential_tilting"

#' @keywords engine_view
#' @export
engine_name.nmar_engine_exptilt_nonparam <- function(x) "exponential_tilting_nonparam"


#' @keywords internal
s3_engine_label <- function(name) {
  switch(name,
    empirical_likelihood = "Empirical Likelihood (EL)",
    exponential_tilting = "Exponential Tilting (ET)",
    exponential_tilting_nonparam = "Exponential Tilting (nonparametric)",
    name
  )
}

#' @keywords internal
s3_is_scalar <- function(x) is.atomic(x) && length(x) == 1L

#' @keywords internal
s3_fmt_val <- function(x) {
  if (is.null(x)) return("<NULL>")
  if (is.numeric(x) && s3_is_scalar(x)) {
    if (is.nan(x)) return("NaN")
    if (is.infinite(x)) return(if (x > 0) "Inf" else "-Inf")
    return(format(x, digits = 6))
  }
  if (is.logical(x) && s3_is_scalar(x)) return(if (x) "TRUE" else "FALSE")
  if (is.character(x) && s3_is_scalar(x)) return(x)
  if (is.list(x) && !is.null(x$name) && is.character(x$name)) return(x$name)
  if (is.function(x)) return("<fn>")
  if (is.atomic(x) && length(x) > 1L) return(sprintf("[%d]", length(x)))
  if (is.list(x)) return(sprintf("<list:%d>", length(x)))
  "<value>"
}

#' @keywords internal
s3_engine_display_keys <- function(x) {
  cls <- class(x)
  out <- c("standardize", "variance_method")
  if ("nmar_engine_el" %in% cls) {
    out <- c(out, "family", "trim_cap")
  } else if ("nmar_engine_exptilt" %in% cls) {
    out <- c(out, "prob_model_type", "y_dens", "optim_method", "min_iter", "max_iter", "tol_value")
  } else if ("nmar_engine_exptilt_nonparam" %in% cls) {
    out <- c(out, "refusal_col", "max_iter", "tol_value")
  }
  out
}

#' Print method for NMAR engines
#'
#' Provides a compact, human-friendly summary for `nmar_engine` objects.
#' Child classes inherit this method; they can override it if they need a
#' different presentation.
#'
#' @param x An engine object inheriting from `nmar_engine`.
#' @param ... Unused.
#' @keywords engine_view
#' @export
print.nmar_engine <- function(x, ...) {
  nm <- engine_name(x)
  cat("NMAR engine: ", s3_engine_label(nm), "\n", sep = "")

  keys <- intersect(s3_engine_display_keys(x), names(x))
  if (length(keys)) {
    for (k in keys) {
      val <- x[[k]]
      if (identical(k, "variance_method") && is.character(val) && val == "bootstrap" && !is.null(x$bootstrap_reps)) {
        cat(sprintf("  %-18s %s (reps = %s)\n", paste0(k, ":"), s3_fmt_val(val), s3_fmt_val(x$bootstrap_reps)))
      } else if (identical(k, "trim_cap") && isTRUE(is.infinite(val))) {
        cat(sprintf("  %-18s %s\n", paste0(k, ":"), "Inf"))
      } else {
        cat(sprintf("  %-18s %s\n", paste0(k, ":"), s3_fmt_val(val)))
      }
    }
  }

# Brief indicator if a custom start was provided
  st <- x$start
  if (!is.null(st) && is.list(st)) {
    kinds <- character()
    if (!is.null(st$beta)) kinds <- c(kinds, "beta")
    if (!is.null(st$z) && !is.null(st$W)) {
      kinds <- c(kinds, "z/W")
    } else if (!is.null(st$z)) {
      kinds <- c(kinds, "z")
    } else if (!is.null(st$W)) {
      kinds <- c(kinds, "W")
    }
    if (!is.null(st$lambda)) kinds <- c(kinds, "lambda")
    if (length(kinds)) {
      cat(sprintf("\nstart: %s\n", paste(kinds, collapse = ", ")))
    }
  }

  cat("\nUse engine_config(x) for full parameters.\n")
  invisible(x)
}

#' One-line formatter for NMAR engines
#'
#' Returns a single concise line summarizing an engine configuration.
#'
#' @param x An engine object inheriting from `nmar_engine`.
#' @param ... Unused.
#' @return A length-1 character vector.
#' @keywords engine_view
#' @export
format.nmar_engine <- function(x, ...) {
  nm <- s3_engine_label(engine_name(x))
  keys <- intersect(s3_engine_display_keys(x), names(x))
  parts <- character()
  for (k in keys) {
    val <- x[[k]]
    v <- if (identical(k, "variance_method") && is.character(val) && val == "bootstrap" && !is.null(x$bootstrap_reps)) {
      sprintf("bootstrap(%s)", s3_fmt_val(x$bootstrap_reps))
    } else if (identical(k, "trim_cap") && isTRUE(is.infinite(val))) {
      "Inf"
    } else {
      s3_fmt_val(val)
    }
    parts <- c(parts, sprintf("%s=%s", k, v))
  }
  paste0(nm, ": ", paste(parts, collapse = ", "))
}

#' Extract engine configuration
#'
#' Returns the underlying configuration of an engine as a named list. This is
#' intended for programmatic inspection (e.g., parameter tuning, logging). The
#' returned object should be treated as read-only.
#'
#' @param x An object inheriting from class `nmar_engine`.
#' @return A named list of configuration fields.
#' @keywords engine_view
#' @export
engine_config <- function(x) {
  UseMethod("engine_config")
}

#' @export
engine_config.nmar_engine <- function(x) {
# Return a shallow copy to discourage in-place mutation of the source object
  out <- x
  attr(out, "class") <- NULL
  out
}
