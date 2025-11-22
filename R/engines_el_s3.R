#' @import generics
#' @importFrom generics tidy glance
#' @importFrom stats weights fitted
NULL

#' Print method for EL results
#' @description Compact print for objects of class `nmar_result_el`.
#' @param x An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @keywords result_view
#' @export
print.nmar_result_el <- function(x, ...) {
  meta <- x$meta %||% list()
# Print an abridged call line instead of the full captured call to avoid
# dumping large data objects
  call_line <- nmar_format_call_line(x)
  if (is.character(call_line) && nzchar(call_line)) {
    cat(call_line, "\n\n", sep = "")
  }

  NextMethod()

  diagnostics <- nmar_result_get_diagnostics(x)
  method_label <- meta$engine_name %||% "Empirical Likelihood (EL)"
  cat("\nMethod: ", method_label, "\n", sep = "")
  if (!isTRUE(x$converged)) {
    msg <- diagnostics$message %||% NA_character_
    if (!is.na(msg)) cat("Convergence message: ", msg, "\n", sep = "")
    return(invisible(x))
  }

  if (!is.null(diagnostics$max_equation_residual)) {
    cat(sprintf("Max equation residual: %.3e\n", diagnostics$max_equation_residual))
  }
  if (!is.null(diagnostics$constraint_sum_W)) {
    cat(sprintf("Constraint sum (W): %.3e\n", diagnostics$constraint_sum_W))
  }
  if (!is.null(diagnostics$constraint_sum_aux) && length(diagnostics$constraint_sum_aux) > 0) {
    cat("Constraint sums (aux):\n")
    print(diagnostics$constraint_sum_aux)
  }
  invisible(x)
}

#' Summary method for EL results
#' @description Summarize estimation, standard error and missingness-model coefficients.
#' @param object An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @keywords result_view
#' @export
summary.nmar_result_el <- function(object, ...) {
  base <- NextMethod()
  model <- nmar_result_get_model(object)
  base$response_model <- model$coefficients
  base$response_vcov <- model$vcov
# Keep an abridged call string for printing; avoid storing the full call
# (which may embed the entire data frame)
  base$call_line <- nmar_format_call_line(object)
  base$df <- nmar_result_get_inference(object)$df
  class(base) <- c("summary_nmar_result_el", class(base))
  base
}

#' @keywords result_view
#' @export
print.summary_nmar_result_el <- function(x, ...) {
  NextMethod()
  if (!is.null(x$call_line) && isTRUE(getOption("nmar.show_call", TRUE))) {
    cat(x$call_line, "\n", sep = "")
  }
  if (!is.null(x$response_model)) {
    cat("\nMissingness-model coefficients:\n")
    beta <- x$response_model
    if (!is.null(x$response_vcov) && is.matrix(x$response_vcov)) {
      se <- sqrt(diag(x$response_vcov))
      stat <- beta / se
      df <- x$df %||% NA_real_
      use_t <- is.finite(df)
      pval <- if (use_t) 2 * stats::pt(-abs(stat), df = df) else 2 * stats::pnorm(-abs(stat))
      stat_label <- if (use_t) "t value" else "z value"
      p_label <- if (use_t) "Pr(>|t|)" else "Pr(>|z|)"
# Apply nmar.digits formatting for visual consistency
      d <- nmar_get_digits()
      tab <- data.frame(
        Estimate = nmar_fmt_num(beta, d),
        `Std. Error` = nmar_fmt_num(se, d),
        check.names = FALSE
      )
      tab[[stat_label]] <- nmar_fmt_num(stat, d)
      tab[[p_label]] <- nmar_fmt_num(pval, d)
      rownames(tab) <- names(beta)
      print(tab, row.names = TRUE)
    } else {
      d <- nmar_get_digits()
      tb <- data.frame(Estimate = nmar_fmt_num(beta, d))
      rownames(tb) <- names(beta)
      print(tb, row.names = TRUE)
    }
  }
  invisible(x)
}

## Engine-specific methods beyond parent defaults are not required here.
