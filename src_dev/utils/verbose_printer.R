#' Create Verbose Printer Factory
#'
#' Creates a verbose printing function based on trace level settings.
#' Messages are printed only if their level is <= trace_level.
#'
#' @param verbose Logical; if FALSE, returns a no-op function
#' @param trace_level Integer 1-3; controls verbosity detail:
#'   - 1: Major steps only (initialization, convergence)
#'   - 2: Moderate detail (iteration summaries, key diagnostics)
#'   - 3: Full detail (all diagnostics, intermediate values)
#'
#' @return A function with signature:
#'   `verboser(msg, level = 1, type = c("info", "step", "detail", "result"))`
#'
#' @keywords internal
create_verboser <- function(verbose = FALSE, trace_level = 1) {
# Validate inputs
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }
  if (!is.numeric(trace_level) || length(trace_level) != 1 || !trace_level %in% 1:3) {
    stop("trace_level must be 1, 2, or 3")
  }

  if (!verbose) {
# Return a no-op function when verbose is FALSE
    return(function(...) invisible(NULL))
  }

# Return the actual verbose printer function when verbose is TRUE
  function(msg = NULL, level = 1, type = c("info", "step", "detail", "result"),
           obj = NULL, max_print = 10) {

    type <- match.arg(type)

# Only print if the message level is <= current trace_level
    if (level > trace_level) return(invisible(NULL))

# Prefix based on type and level
    prefix <- switch(type,
      "info"   = "[INFO]",
      "step"   = paste0("[STEP-", level, "]"),
      "detail" = paste0("[DETAIL-", level, "]"),
      "result" = "[RESULT]"
    )

# Print message if provided
    if (!is.null(msg)) {
      if (is.character(msg) && length(msg) == 1) {
        cat(prefix, msg, "\n")
      } else {
        cat(prefix, "\n")
        print(msg)
      }
    }

# Print object if provided
    if (!is.null(obj)) {
      if (is.numeric(obj) && length(obj) == 1) {
# Single numeric value
        cat("  Value:", format(obj, digits = 6), "\n")
      } else if (is.vector(obj) && length(obj) <= max_print) {
# Short vector - print all
        cat("  Values:", paste(format(obj, digits = 4), collapse = ", "), "\n")
      } else if (is.vector(obj)) {
# Long vector - print summary
        cat("  Length:", length(obj), "\n")
        cat("  Range: [", format(min(obj, na.rm = TRUE), digits = 4), ",",
            format(max(obj, na.rm = TRUE), digits = 4), "]\n")
        cat("  Mean:", format(mean(obj, na.rm = TRUE), digits = 4), "\n")
      } else if (is.matrix(obj)) {
# Matrix - print dimensions and summary
        cat("  Dimensions:", nrow(obj), "x", ncol(obj), "\n")
        if (nrow(obj) * ncol(obj) <= max_print) {
          print(obj)
        } else {
          cat("  Range: [", format(min(obj, na.rm = TRUE), digits = 4), ",",
              format(max(obj, na.rm = TRUE), digits = 4), "]\n")
        }
      } else {
# Other objects - use default print
        print(obj)
      }
    }

    invisible(NULL)
  }
}
