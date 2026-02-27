#' @import generics
#' @importFrom generics tidy glance
#' @importFrom stats fitted
NULL

#' Print method for induced-logit results
#'
#' @description Print for objects of class `nmar_result_induced_logit`.
#' @param x An object of class `nmar_result_induced_logit`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @keywords result_view
#' @export
print.nmar_result_induced_logit <- function(x, ...) {
  meta <- x$meta %||% list()
  call_line <- nmar_format_call_line(x)
  if (is.character(call_line) && nzchar(call_line)) {
    cat(call_line, "\n\n", sep = "")
  }

  NextMethod()

  diagnostics <- nmar_result_get_diagnostics(x)
  method_label <- "Induced Logistic Regression"
  method_id <- meta$engine_name %||% "induced_logistic"
  cat("\nMethod: ", method_label, " (", method_id, ")\n", sep = "")

  if (!isTRUE(x$converged)) {
    msg <- diagnostics$message %||% nmar_result_get_inference(x)$message %||% NA_character_
    if (!is.na(msg)) cat("Convergence message: ", msg, "\n", sep = "")
    return(invisible(x))
  }

  d <- nmar_get_digits()
  eta_hat <- diagnostics$eta_hat %||% NULL
  if (!is.null(eta_hat) && is.finite(eta_hat)) {
    cat(sprintf("Response rate (eta_hat): %s\n", nmar_fmt_num(eta_hat, d)))
  }
  gamma_hat <- diagnostics$gamma_hat_paper %||% NULL
  if (!is.null(gamma_hat) && is.finite(gamma_hat)) {
    cat(sprintf("Gamma (paper; = -coef(mu_hat)): %s\n", nmar_fmt_num(gamma_hat, d)))
  }
  alpha0_hat <- diagnostics$alpha0_hat_paper %||% NULL
  if (!is.null(alpha0_hat) && is.finite(alpha0_hat)) {
    cat(sprintf("Alpha0 (paper; derived): %s\n", nmar_fmt_num(alpha0_hat, d)))
  }
  m2_over_m1 <- diagnostics$m2_over_m1 %||% NULL
  if (!is.null(m2_over_m1) && is.finite(m2_over_m1)) {
    cat(sprintf("M2/M1: %s\n", nmar_fmt_num(m2_over_m1, d)))
  }

  warn <- diagnostics$warnings %||% list()
  n_mu <- length(warn$mu %||% character())
  n_glm <- length(warn$induced_glm %||% character())
  if (n_mu > 0L || n_glm > 0L) {
    cat(sprintf("Warnings captured: mu=%d, induced_glm=%d\n", n_mu, n_glm))
  }

  invisible(x)
}

#' Summary method for induced-logit results
#'
#' @description Summarize estimation, standard error and missingness-model coefficients.
#' @param object An object of class `nmar_result_induced_logit`.
#' @param ... Ignored.
#' @return An object of class `summary_nmar_result_induced_logit`.
#' @keywords result_view
#' @export
summary.nmar_result_induced_logit <- function(object, ...) {
  base <- NextMethod()
  model <- nmar_result_get_model(object)
  diagnostics <- nmar_result_get_diagnostics(object)
  base$response_model <- model$coefficients
  base$response_vcov <- model$vcov
  base$call_line <- nmar_format_call_line(object)
  base$df <- nmar_result_get_inference(object)$df
  base$gamma_hat_paper <- diagnostics$gamma_hat_paper %||% NA_real_
  class(base) <- c("summary_nmar_result_induced_logit", class(base))
  base
}

#' @keywords result_view
#' @export
print.summary_nmar_result_induced_logit <- function(x, ...) {
  NextMethod()
  if (!is.null(x$call_line) && isTRUE(getOption("nmar.show_call", TRUE))) {
    cat(x$call_line, "\n", sep = "")
  }

  diag <- x$diagnostics %||% list()
  d <- nmar_get_digits()
  cat("\nInduced-logit diagnostics:\n")
  eta_hat <- diag$eta_hat %||% NULL
  if (!is.null(eta_hat) && is.finite(eta_hat)) {
    cat("  eta_hat:", nmar_fmt_num(eta_hat, d), "\n")
  }
  mu_bar <- diag$mu_bar %||% NULL
  if (!is.null(mu_bar) && is.finite(mu_bar)) {
    cat("  mu_bar:", nmar_fmt_num(mu_bar, d), "\n")
  }
  m2_over_m1 <- diag$m2_over_m1 %||% NULL
  if (!is.null(m2_over_m1) && is.finite(m2_over_m1)) {
    cat("  M2/M1:", nmar_fmt_num(m2_over_m1, d), "\n")
  }
  gamma_hat <- diag$gamma_hat_paper %||% NULL
  if (!is.null(gamma_hat) && is.finite(gamma_hat)) {
    cat("  gamma (paper):", nmar_fmt_num(gamma_hat, d), "\n")
  }
  alpha0_hat <- diag$alpha0_hat_paper %||% NULL
  if (!is.null(alpha0_hat) && is.finite(alpha0_hat)) {
    cat("  alpha0 (paper; derived):", nmar_fmt_num(alpha0_hat, d), "\n")
  }

  if (!is.null(x$response_model)) {
    cat("\nMissingness-model coefficients (GLM, logit Pr(R=1|x)):\n")
    cat("  Note: paper parameters have the opposite sign; gamma_paper = -coef(mu_hat).\n")
    beta <- x$response_model
    nm <- names(beta) %||% rep("", length(beta))
    nm <- sub("^\\.\\.nmar_mu_hat\\.\\.$", "mu_hat", nm)
    names(beta) <- nm

    if (!is.null(x$response_vcov) && is.matrix(x$response_vcov)) {
      vc <- x$response_vcov
      rn <- rownames(vc) %||% names(beta)
      cn <- colnames(vc) %||% names(beta)
      rownames(vc) <- sub("^\\.\\.nmar_mu_hat\\.\\.$", "mu_hat", rn)
      colnames(vc) <- sub("^\\.\\.nmar_mu_hat\\.\\.$", "mu_hat", cn)

      se <- sqrt(diag(vc))
      stat <- beta / se
      df <- x$df %||% NA_real_
      use_t <- is.finite(df)
      pval <- if (use_t) 2 * stats::pt(-abs(stat), df = df) else 2 * stats::pnorm(-abs(stat))
      stat_label <- if (use_t) "t value" else "z value"
      p_label <- if (use_t) "Pr(>|t|)" else "Pr(>|z|)"

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
      tb <- data.frame(Estimate = nmar_fmt_num(beta, d))
      rownames(tb) <- names(beta)
      print(tb, row.names = TRUE)
    }

    if (is.finite(x$gamma_hat_paper)) {
      cat(sprintf("\nGamma (paper; = -coef(mu_hat)): %s\n", nmar_fmt_num(x$gamma_hat_paper, d)))
    }
  }

  invisible(x)
}
