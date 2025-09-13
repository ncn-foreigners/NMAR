#' @import generics
#' @importFrom generics tidy glance
#' @importFrom stats weights fitted
NULL

#' Print method for EL results
#' @description Compact print for objects of class `nmar_result_el`.
#' @param x An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @export
print.nmar_result_el <- function(x, ...) {
  cat("Call:\n")
  if (!is.null(x$call)) print(x$call)
  cat("\n--- NMAR Estimation Result ---\n")
  method <- if (!is.null(x$data_info$method)) x$data_info$method else "Empirical Likelihood (EL)"
  cat("Method:", method, "\n\n")
  if (isTRUE(x$converged)) {
    nm <- x$data_info$outcome_var %||% "y"
    cat("Population Mean Estimate\n")
    est <- x$y_hat
    names(est) <- nm
    print(est)
  } else {
    cat("--> Estimation failed to converge.\n")
    if (!is.null(x$diagnostics$message)) cat("--> Message:", x$diagnostics$message, "\n")
  }
  invisible(x)
}

#' Summary method for EL results
#' @description Summarize estimation, standard error and response-model coefficients.
#' @param object An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @export
summary.nmar_result_el <- function(object, ...) {
  if (!isTRUE(object$converged)) {
    cat("Model did not converge.\n")
    return(invisible(object))
  }
  is_svy <- isTRUE(object$data_info$is_survey)
  df <- if (is_svy && !is.null(object$data_info$design)) tryCatch(survey::degf(object$data_info$design), error = function(e) Inf) else Inf
  crit <- if (is.finite(df) && is_svy) stats::qt(0.975, df = df) else stats::qnorm(0.975)
  ci <- object$y_hat + c(-1, 1) * crit * object$se

  cat("Call:\n")
  if (!is.null(object$call)) print(object$call)
  cat("\n")
  cat("--- Population Mean Estimate ---\n")
  cat(sprintf("Estimate of mean(%s): %.4f\n", object$data_info$outcome_var, object$y_hat))
  cat(sprintf("Std. Error: %.4f\n", object$se))
  vm <- object$data_info$variance_method %||% "unknown"
  cat(sprintf("95%% CI (Wald, %s): (%.4f, %.4f)\n", vm, ci[1], ci[2]))

  # Coefficients (if available)
  if (is.list(object$coefficients) && !is.null(object$coefficients$response_model)) {
    cat("\n--- Response Model Coefficients ---\n")
    beta <- object$coefficients$response_model
    if (is.matrix(object$vcov)) {
      se <- sqrt(diag(object$vcov))
      stat <- beta / se
      p <- if (is_svy) 2 * stats::pt(-abs(stat), df = df) else 2 * stats::pnorm(-abs(stat))
      out <- data.frame(Estimate = beta, `Std. Error` = se, `z/t value` = stat, `Pr(>|z/t|)` = p, check.names = FALSE)
      print(out)
    } else {
      print(data.frame(Estimate = beta))
    }
  }
  invisible(object)
}

#' Estimate for EL results
#' @param x An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @export
estimate.nmar_result_el <- function(x, ...) {
  est <- x$y_hat
  nm <- x$data_info$outcome_var %||% "y"
  names(est) <- nm
  est
}

#' @export
coef.nmar_result_el <- function(object, type = c("response", "estimand"), ...) {
  type <- match.arg(type)
  if (type == "response") {
    return(object$coefficients$response_model)
  }
  est <- object$y_hat
  names(est) <- object$data_info$outcome_var
  est
}

#' Variance-covariance for the primary estimand (EL)
#' @description Variance-covariance for the primary estimand (population mean).
#' @details Returns a 1x1 matrix for Var(\eqn{\hat Y}). The covariance matrix
#'   for the response-model coefficients is stored in `object$vcov` and is used
#'   by `tidy()`; it is not returned by this `vcov()` method.
#' @param object An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @export
vcov.nmar_result_el <- function(object, ...) {
  if (!isTRUE(object$converged) || !is.finite(object$se)) {
    mat <- matrix(NA_real_, 1, 1)
  } else {
    mat <- matrix(object$se^2, 1, 1)
  }
  dimnames(mat) <- list(object$data_info$outcome_var, object$data_info$outcome_var)
  mat
}

#' Confidence interval for the primary estimand (EL)
#' @description Wald confidence interval for the primary estimand (population mean).
#' @details Uses normal quantiles for IID data and t-quantiles with survey
#'   design degrees-of-freedom for `survey.design` inputs.
#' @param object An object of class `nmar_result_el`.
#' @param parm Not used; included for S3 compatibility.
#' @param level Confidence level in (0,1). Default 0.95.
#' @param ... Ignored.
#' @export
confint.nmar_result_el <- function(object, parm, level = 0.95, ...) {
  if (!isTRUE(object$converged) || !is.finite(object$se)) {
    ci <- c(NA_real_, NA_real_)
  } else {
    is_svy <- isTRUE(object$data_info$is_survey)
    df <- if (is_svy && !is.null(object$data_info$design)) tryCatch(survey::degf(object$data_info$design), error = function(e) Inf) else Inf
    alpha <- 1 - level
    crit <- if (is_svy && is.finite(df)) stats::qt(1 - alpha / 2, df = df) else stats::qnorm(1 - alpha / 2)
    ci <- object$y_hat + c(-1, 1) * crit * object$se
  }
  m <- matrix(ci, nrow = 1)
  colnames(m) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), " %")
  rownames(m) <- object$data_info$outcome_var
  m
}

#' Fitted probabilities for EL results
#' @param object An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @export
fitted.nmar_result_el <- function(object, ...) {
  fv <- object$fitted_values
  if (is.null(fv) || length(fv) == 0) {
    return(numeric(0))
  }
  as.numeric(fv)
}

#' Weights for EL results
#' @param object An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @export
weights.nmar_result_el <- function(object, ...) {
  w <- object$weights
  if (is.null(w)) {
    return(numeric(0))
  }
  w <- as.numeric(w)
  tr <- object$diagnostics$trimmed_fraction
  attr(w, "trimmed_fraction") <- tr
  w
}

#' Estimation formula for EL results
#' @param x An object of class `nmar_result_el`.
#' @param ... Ignored.
#' @export
formula.nmar_result_el <- function(x, ...) {
  if (!is.null(x$data_info$formula)) {
    return(x$data_info$formula)
  }
  NULL
}

# small infix helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Additional S3 methods for user workflows

#' @export
autoplot.nmar_result <- function(object, type = c("weights", "fitted", "constraints"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for autoplot.nmar_result", call. = FALSE)
  if (type == "weights") {
    w <- weights(object)
    df <- data.frame(w = as.numeric(w))
    ggplot2::ggplot(df, ggplot2::aes(x = w)) +
      ggplot2::geom_histogram(color = "white", fill = "gray") +
      ggplot2::labs(title = "EL weights", x = "weight")
  } else if (type == "fitted") {
    fv <- fitted(object)
    df <- data.frame(p = as.numeric(fv))
    ggplot2::ggplot(df, ggplot2::aes(x = p)) +
      ggplot2::geom_histogram(color = "white", fill = "gray") +
      ggplot2::labs(title = "Fitted response probabilities", x = "p_hat")
  } else {
    vals <- c(eqW = object$diagnostics$constraint_sum_W)
    if (!is.null(object$diagnostics$constraint_sum_aux) && length(object$diagnostics$constraint_sum_aux) > 0) vals <- c(vals, object$diagnostics$constraint_sum_aux)
    df <- data.frame(term = names(vals), value = as.numeric(vals))
    ggplot2::ggplot(df, ggplot2::aes(x = term, y = value)) +
      ggplot2::geom_col(fill = "gray") +
      ggplot2::geom_hline(yintercept = 0, color = "darkgray") +
      ggplot2::labs(title = "Constraint sums", x = NULL, y = "sum") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
}
