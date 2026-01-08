#' @import generics
#' @importFrom generics tidy glance
#' @importFrom stats weights fitted
NULL
## Exptilt-specific print/summary implementations (src_dev)

#' Print method for Exponential Tilting results (engine-specific)
#'
#' This print method is tailored for `nmar_result_exptilt` objects and shows a
#' concise, human-friendly summary of the estimation result together with
#' exptilt-specific diagnostics (loss, iterations) and a compact view of the
#' response coefficients stored in the fitted model.
#'
#' @param x An object of class `nmar_result_exptilt`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @keywords result_view
#' @export
print.nmar_result_exptilt <- function(x, ...) {
  est <- nmar_result_get_estimate(x)
  se <- nmar_result_get_se(x)
  nm <- nmar_result_get_estimate_name(x)
  inference <- nmar_result_get_inference(x)
  sample <- nmar_result_get_sample(x)
  meta <- x$meta %||% list()
  diagnostics <- nmar_result_get_diagnostics(x)

  cat("NMAR Result (Exponential tilting)\n")
  cat("-------------------------------\n")
  d <- nmar_get_digits()
  if (is.finite(est) && is.finite(se)) {
    cat(sprintf("%s mean: %s (%s)\n", nm, nmar_fmt_num(est, d), nmar_fmt_num(se, d)))
  } else if (is.finite(est)) {
    cat(sprintf("%s mean: %s\n", nm, nmar_fmt_num(est, d)))
  } else {
    cat(sprintf("%s mean: NA\n", nm))
  }
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

# Sampling info
  if (isTRUE(diagnostics$sampling_performed)) {
    cat("Stratified sampling:\n")
    cat(sprintf("  Original size: %d (resp: %d, non-resp: %d)\n",
                diagnostics$original_n_total,
                diagnostics$original_n_resp,
                diagnostics$original_n_nonresp))
    cat(sprintf("  Sampled size: %d (preserving ratio)\n", sample$n_total))
  }

# Exptilt-specific diagnostics
  cat("\nExptilt diagnostics:\n")
  if (!is.null(diagnostics$loss_value)) cat("  Loss value:", nmar_fmt_num(diagnostics$loss_value, d), "\n")
  if (!is.null(diagnostics$iterations)) cat("  Iterations:", as.integer(diagnostics$iterations), "\n")
  if (!is.null(diagnostics$variance_method)) cat("  Variance method:", diagnostics$variance_method, "\n")
  if (!is.null(diagnostics$bootstrap_reps)) cat("  Bootstrap reps:", diagnostics$bootstrap_reps, "\n")
  if (!is.null(diagnostics$stopping_threshold)) cat("  Stopping threshold:", nmar_fmt_num(diagnostics$stopping_threshold, d), "\n")

  invisible(x)
}


#' Summary method for Exponential Tilting results (engine-specific)
#' @description Summarize estimation, standard error and model coefficients.
#' @param object An object of class `nmar_result_exptilt`.
#' @param conf.level Confidence level for confidence interval (default 0.95).
#' @param ... Ignored.
#' @return An object of class `summary_nmar_result_exptilt`.
#' @keywords result_view
#' @export
summary.nmar_result_exptilt <- function(object, conf.level = 0.95, ...) {
  est <- nmar_result_get_estimate(object)
  se <- nmar_result_get_se(object)
  nm <- nmar_result_get_estimate_name(object)
  inference <- nmar_result_get_inference(object)
  sample <- nmar_result_get_sample(object)
  diagnostics <- nmar_result_get_diagnostics(object)

  ci <- confint(object, level = conf.level)

  model <- nmar_result_get_model(object)

  out <- list(
    y_hat = as.numeric(est),
    estimate_name = nm,
    se = se,
    conf_int = ci,
    converged = isTRUE(object$converged),
    variance_method = inference$variance_method,
    variance_message = inference$message,
    sample = sample,
    diagnostics = diagnostics,
    meta = object$meta %||% list(),
    conf.level = conf.level
  )

# exptilt-specific attachments
  theta_val <- model$coefficients %||% NULL
  if (is.null(theta_val) && !is.null(object$extra) && is.list(object$extra$raw)) {
    theta_val <- object$extra$raw$model$theta %||% NULL
  }
  out$theta <- theta_val
  out$theta_vcov <- model$vcov
  out$call_line <- nmar_format_call_line(object)
  out$df <- inference$df

  class(out) <- c("summary_nmar_result_exptilt", "summary_nmar_result")
  out
}


#' @keywords result_view
#' @export
print.summary_nmar_result_exptilt <- function(x, ...) {
# Print base summary content
  cat("NMAR Model Summary (Exponential tilting)\n")
  cat("=================================\n")
  d <- nmar_get_digits()
  if (is.finite(x$y_hat) && is.finite(x$se)) {
    cat(sprintf("%s mean: %s (%s)\n", x$estimate_name, nmar_fmt_num(x$y_hat, d), nmar_fmt_num(x$se, d)))
  } else if (is.finite(x$y_hat)) {
    cat(sprintf("%s mean: %s\n", x$estimate_name, nmar_fmt_num(x$y_hat, d)))
  } else {
    cat(sprintf("%s mean: NA\n", x$estimate_name))
  }
  if (!anyNA(x$conf_int)) {
    cat(sprintf("%g%% CI: (%s, %s)\n", 100 * x$conf.level, nmar_fmt_num(x$conf_int[1, 1], d), nmar_fmt_num(x$conf_int[1, 2], d)))
  }
  cat("Converged:", x$converged, "\n")
  if (!is.null(x$variance_method) && !is.na(x$variance_method)) cat("Variance method:", x$variance_method, "\n")
  if (!is.null(x$variance_message) && !is.na(x$variance_message)) cat("Variance notes:", x$variance_message, "\n")
  if (is.list(x$sample)) {
    if (is.finite(x$sample$n_total)) cat("Total units:", x$sample$n_total, "\n")
    if (is.finite(x$sample$n_respondents)) cat("Respondents:", x$sample$n_respondents, "\n")
  }

  if (!is.null(x$call_line) && isTRUE(getOption("nmar.show_call", TRUE))) {
    cat(x$call_line, "\n", sep = "")
  }

# Theta coefficients (response model) - print inline for clarity
  theta <- x$theta
  if (!is.null(theta) && length(theta) > 0) {
    cat("\nResponse-model (theta) coefficients:\n")
    d <- nmar_get_digits()
    if (!is.null(x$theta_vcov) && is.matrix(x$theta_vcov)) {
      se <- sqrt(diag(x$theta_vcov))
      stat <- as.numeric(theta) / se
      df <- x$df %||% NA_real_
      use_t <- is.finite(df)
      pval <- if (use_t) 2 * stats::pt(-abs(stat), df = df) else 2 * stats::pnorm(-abs(stat))
      for (i in seq_along(theta)) {
        nm <- names(theta)[i] %||% paste0("theta", i)
        cat(sprintf("  %-20s : %s (SE=%s, %s=%s, p=%s)\n",
                    nm,
                    nmar_fmt_num(as.numeric(theta[i]), d),
                    nmar_fmt_num(se[i], d),
                    if (use_t) "t" else "z",
                    nmar_fmt_num(stat[i], d),
                    nmar_fmt_num(pval[i], d)))
      }
    } else {
      for (i in seq_along(theta)) {
        nm <- names(theta)[i] %||% paste0("theta", i)
        cat(sprintf("  %-20s : %s\n", nm, nmar_fmt_num(as.numeric(theta[i]), d)))
      }
    }
  } else {
    cat("\nNo response-model coefficients available.\n")
  }

  invisible(x)
}
