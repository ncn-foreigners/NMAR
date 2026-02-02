#' Create Verbose Printer Factory
#'
#' Creates a verbose printing function based on trace level settings.
#' Messages are printed only if their level is <= trace_level.
#'
#' @param trace_level Integer 0-3; controls verbosity detail:
#'   - 0: No output (silent mode)
#'   - 1: Major steps only (initialization, convergence)
#'   - 2: Moderate detail (iteration summaries, key diagnostics)
#'   - 3: Full detail (all diagnostics, intermediate values)
#'
#' @return A function with signature:
#'   `verboser(msg, level = 1, type = c("info", "step", "detail", "result"))`
#'
#' @keywords internal
create_verboser <- function(trace_level = 0) {
  validator_assert_choice(trace_level, choices = 0:3, name = "trace_level")

  if (trace_level == 0) {
# no display when trace_level is 0
    return(function(...) invisible(NULL))
  }

# printer function when trace_level > 0
  function(msg = NULL, level = 1, type = c("info", "step", "detail", "result"),
           obj = NULL, max_print = 10) {

    type <- match.arg(type)

# only print if the message level is <= current trace_level
    if (level > trace_level) return(invisible(NULL))

    prefix <- switch(type,
      "info"   = "[INFO]",
      "step"   = paste0("[STEP-", level, "]"),
      "detail" = paste0("[DETAIL-", level, "]"),
      "result" = "[RESULT]"
    )

    if (!is.null(msg)) {
      if (is.character(msg) && length(msg) == 1) {
        cat(prefix, msg, "\n")
      } else {
        cat(prefix, "\n")
        print(msg)
      }
    }


    if (!is.null(obj)) {
      if (is.numeric(obj) && length(obj) == 1) {

        cat("  Value:", format(obj, digits = 6), "\n")
      } else if (is.vector(obj) && length(obj) <= max_print) {

        cat("  Values:", paste(format(obj, digits = 4), collapse = ", "), "\n")
      } else if (is.vector(obj)) {

        cat("  Length:", length(obj), "\n")
        cat("  Range: [", format(min(obj, na.rm = TRUE), digits = 4), ",",
            format(max(obj, na.rm = TRUE), digits = 4), "]\n")
        cat("  Mean:", format(mean(obj, na.rm = TRUE), digits = 4), "\n")
      } else if (is.matrix(obj)) {

        cat("  Dimensions:", nrow(obj), "x", ncol(obj), "\n")
        if (nrow(obj) * ncol(obj) <= max_print) {
          print(obj)
        } else {
          cat("  Range: [", format(min(obj, na.rm = TRUE), digits = 4), ",",
              format(max(obj, na.rm = TRUE), digits = 4), "]\n")
        }
      } else {

        print(obj)
      }
    }

    invisible(NULL)
  }
}
