#' Print method for nmar_result
#'
#' @param x nmar_result object
#' @param ... Additional parameters
#' @export

print.nmar_result <- function(x, ...) {
  est <- nmar_result_get_estimate(x)
  se <- nmar_result_get_std_error(x)
  nm <- nmar_result_get_estimate_name(x)
  inference <- nmar_result_get_inference(x)
  sample <- nmar_result_get_sample(x)
  meta <- x$meta %||% list()

  cat("NMAR Result\n")
  cat("------------\n")
  if (is.finite(est)) {
    cat(sprintf("%s: %0.6f\n", nm, est))
  } else {
    cat(sprintf("%s: NA\n", nm))
  }
  if (is.finite(se)) cat(sprintf("Std. Error: %0.6f\n", se))
  cat("Converged:", isTRUE(x$converged), "\n")
  if (!is.null(inference$variance_method)) {
    cat("Variance method:", inference$variance_method, "\n")
  }
  if (!is.null(meta$engine_name)) {
    cat("Estimator:", meta$engine_name, "\n")
  }
  if (is.finite(sample$n_total)) {
    cat("Sample size:", sample$n_total)
    if (is.finite(sample$n_respondents)) {
      cat(sprintf(" (respondents: %d)", sample$n_respondents))
    }
    cat("\n")
  }
  invisible(x)
}

#' Summary method for nmar_result
#'
#' @param object nmar_result object
#' @param conf.level Confidence level for intervals.
#' @param ... Additional parameters
#' @export

summary.nmar_result <- function(object, conf.level = 0.95, ...) {
  est <- nmar_result_get_estimate(object)
  se <- nmar_result_get_std_error(object)
  nm <- nmar_result_get_estimate_name(object)
  inference <- nmar_result_get_inference(object)
  sample <- nmar_result_get_sample(object)
  diagnostics <- object$diagnostics %||% list()
  ci <- confint(object, level = conf.level)

  structure(
    list(
      estimate = as.numeric(est),
      estimate_name = nm,
      std_error = se,
      conf_int = ci,
      converged = isTRUE(object$converged),
      variance_method = inference$variance_method,
      variance_message = inference$message,
      sample = sample,
      diagnostics = diagnostics,
      meta = object$meta %||% list(),
      conf.level = conf.level
    ),
    class = "summary_nmar_result"
  )
}

#' Print method for summary.nmar_result
#'
#' @param x summary_nmar_result object
#' @param ... Additional parameters
#' @export

print.summary_nmar_result <- function(x, ...) {
  cat("NMAR Model Summary\n")
  cat("=================\n")
  cat(sprintf("%s estimate: %0.6f\n", x$estimate_name, x$estimate))
  if (is.finite(x$std_error)) {
    cat(sprintf("Std. Error: %0.6f\n", x$std_error))
  }
  if (!anyNA(x$conf_int)) {
    cat(sprintf("%g%% CI: (%0.6f, %0.6f)\n", 100 * x$conf.level, x$conf_int[1, 1], x$conf_int[1, 2]))
  }
  cat("Converged:", x$converged, "\n")
  if (!is.null(x$variance_method) && !is.na(x$variance_method)) {
    cat("Variance method:", x$variance_method, "\n")
  }
  if (!is.null(x$variance_message) && !is.na(x$variance_message)) {
    cat("Variance notes:", x$variance_message, "\n")
  }
  if (is.list(x$sample)) {
    if (is.finite(x$sample$n_total)) {
      cat("Total units:", x$sample$n_total, "\n")
    }
    if (is.finite(x$sample$n_respondents)) {
      cat("Respondents:", x$sample$n_respondents, "\n")
    }
  }
  invisible(x)
}
