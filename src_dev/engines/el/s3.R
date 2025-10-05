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
  call_obj <- meta$call %||% x$call
  if (!is.null(call_obj)) {
    cat("Call:\n")
    print(call_obj)
    cat("\n")
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
#' @description Summarize estimation, standard error and response-model coefficients.
#' @param object An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @keywords result_view
#' @export
summary.nmar_result_el <- function(object, ...) {
  base <- NextMethod()
  model <- nmar_result_get_model(object)
  base$response_model <- model$coefficients
  base$response_vcov <- model$vcov
  base$call <- object$meta$call %||% object$call
  base$df <- nmar_result_get_inference(object)$df
  class(base) <- c("summary_nmar_result_el", class(base))
  base
}

#' @keywords result_view
#' @export
print.summary_nmar_result_el <- function(x, ...) {
  NextMethod()
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  if (!is.null(x$response_model)) {
    cat("\nResponse-model coefficients:\n")
    beta <- x$response_model
    if (!is.null(x$response_vcov) && is.matrix(x$response_vcov)) {
      se <- sqrt(diag(x$response_vcov))
      stat <- beta / se
      df <- x$df %||% NA_real_
      if (is.finite(df)) {
        pval <- 2 * stats::pt(-abs(stat), df = df)
      } else {
        pval <- 2 * stats::pnorm(-abs(stat))
      }
      tab <- data.frame(Estimate = beta, `Std. Error` = se, `z value` = stat, `Pr(>|z|)` = pval, check.names = FALSE)
      print(tab)
    } else {
      print(data.frame(Estimate = beta))
    }
  }
  invisible(x)
}

## Engine-specific methods beyond parent defaults are not required here.
