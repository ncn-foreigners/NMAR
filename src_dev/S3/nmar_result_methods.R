#' S3 surface for base `nmar_result`
#'
#' @title Parent S3 surface for NMAR results
#' @description Methods that apply to the parent `nmar_result` class and are not
#'   specific to a particular engine (e.g., EL). Engines return a child class
#'   (e.g., `nmar_result_el`) that inherits from `nmar_result` and may override
#'   or extend behavior.
#' @details
#'   Result objects expose a universal schema:
#'   - `estimate`, `estimate_name`, `std_error`, `converged`.
#'   - `model`: list with `coefficients`, `vcov`, plus optional extras.
#'   - `weights_info`: list with respondent weights and trimming metadata.
#'   - `sample`: list with total units, respondent count, survey flag, and `design`.
#'   - `inference`: variance metadata (`variance_method`, `df`, diagnostic flags).
#'   - `diagnostics`, `meta`, and `extra` for estimator-specific details.
#'
#'   New engines should populate these components in their constructors and rely
#'   on the `nmar_result_get_*` utilities when implementing child-specific S3
#'   methods.
#'
#' @name nmar_result_s3_parent
#' @keywords internal
#' @import generics
#' @importFrom generics tidy glance
#' @importFrom stats fitted weights coef
NULL

#' Extract the primary estimate
#' @description Generic to extract the quantity of interest (typically the population mean).
#' @param x A fitted object.
#' @param ... Ignored.
#' @export
estimate <- function(x, ...) UseMethod("estimate")

#' Estimate for base NMAR results
#' @param x An object of class `nmar_result`.
#' @param ... Ignored.
#' @export
estimate.nmar_result <- function(x, ...) {
  est <- nmar_result_get_estimate(x)
  nm <- nmar_result_get_estimate_name(x)
  if (is.na(est)) return(NA_real_)
  setNames(as.numeric(est), nm)
}

#' Variance-covariance for base NMAR results
#' @param object An object of class `nmar_result`.
#' @param ... Ignored.
#' @export
vcov.nmar_result <- function(object, ...) {
  se <- nmar_result_get_std_error(object)
  if (length(se) == 1 && is.finite(se)) {
    mat <- matrix(as.numeric(se)^2, 1, 1)
  } else {
    mat <- matrix(NA_real_, 1, 1)
  }
  nm <- nmar_result_get_estimate_name(object)
  dimnames(mat) <- list(nm, nm)
  mat
}

#' Wald confidence interval for base NMAR results
#' @param object An object of class `nmar_result`.
#' @param parm Ignored.
#' @param level Confidence level.
#' @param ... Ignored.
#' @export
confint.nmar_result <- function(object, parm, level = 0.95, ...) {
  se <- nmar_result_get_std_error(object)
  nm <- nmar_result_get_estimate_name(object)
  est <- nmar_result_get_estimate(object)
  if (!is.finite(se)) {
    return(matrix(c(NA_real_, NA_real_), nrow = 1, dimnames = list(nm, c("2.5 %", "97.5 %"))))
  }
  alpha <- 1 - level
  inference <- nmar_result_get_inference(object)
  sample <- nmar_result_get_sample(object)
  if (isTRUE(sample$is_survey) && is.finite(inference$df)) {
    crit <- stats::qt(1 - alpha / 2, df = inference$df)
  } else {
    crit <- stats::qnorm(1 - alpha / 2)
  }
  ci <- as.numeric(est) + c(-1, 1) * crit * se
  m <- matrix(ci, nrow = 1)
  colnames(m) <- paste0(format(100 * c(alpha / 2, 1 - alpha / 2)), " %")
  rownames(m) <- nm
  m
}

#' Tidy summary for NMAR results
#' @description Return a data frame with the primary estimate and (if available) response-model coefficients.
#' @param x An object of class `nmar_result`.
#' @param conf.level Confidence level for the primary estimate.
#' @param ... Ignored.
#' @export
tidy.nmar_result <- function(x, conf.level = 0.95, ...) {
  est <- nmar_result_get_estimate(x)
  se <- nmar_result_get_std_error(x)
  nm <- nmar_result_get_estimate_name(x)
  inference <- nmar_result_get_inference(x)
  sample <- nmar_result_get_sample(x)
  if (isTRUE(sample$is_survey) && is.finite(inference$df)) {
    crit <- stats::qt(1 - (1 - conf.level) / 2, df = inference$df)
  } else {
    crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  }
  ci <- c(NA_real_, NA_real_)
  if (is.finite(se)) ci <- as.numeric(est + c(-1, 1) * crit * se)
  rows <- list(data.frame(
    term = nm,
    estimate = as.numeric(est),
    std.error = se,
    conf.low = ci[1],
    conf.high = ci[2],
    component = "estimand",
    statistic = NA_real_,
    p.value = NA_real_,
    check.names = FALSE
  ))

  model <- nmar_result_get_model(x)
  beta <- model$coefficients
  if (!is.null(beta) && length(beta) > 0) {
    beta_vec <- as.numeric(beta)
    se_beta <- rep(NA_real_, length(beta_vec))
    stat <- pval <- rep(NA_real_, length(beta_vec))
    if (!is.null(model$vcov) && is.matrix(model$vcov)) {
      se_beta <- sqrt(diag(model$vcov))
      stat <- beta_vec / se_beta
      if (isTRUE(sample$is_survey) && is.finite(inference$df)) {
        pval <- 2 * stats::pt(-abs(stat), df = inference$df)
      } else {
        pval <- 2 * stats::pnorm(-abs(stat))
      }
    }
    beta_names <- names(beta)
    if (is.null(beta_names) || length(beta_names) != length(beta_vec) || anyNA(beta_names)) {
      beta_names <- paste0("coef", seq_along(beta_vec))
    }
    rows[[2]] <- data.frame(
      term = beta_names,
      estimate = beta_vec,
      std.error = se_beta,
      statistic = stat,
      p.value = pval,
      conf.low = rep(NA_real_, length(beta_vec)),
      conf.high = rep(NA_real_, length(beta_vec)),
      component = rep("response", length(beta_vec)),
      check.names = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Glance summary for NMAR results
#' @description One-row diagnostics for NMAR fits.
#' @param x An object of class `nmar_result`.
#' @param ... Ignored.
#' @export
glance.nmar_result <- function(x, ...) {
  est <- nmar_result_get_estimate(x)
  se <- nmar_result_get_std_error(x)
  inference <- nmar_result_get_inference(x)
  sample <- nmar_result_get_sample(x)
  diagnostics <- nmar_result_get_diagnostics(x)
  if (isTRUE(sample$is_survey) && is.finite(inference$df)) {
    crit <- stats::qt(0.975, df = inference$df)
  } else {
    crit <- stats::qnorm(0.975)
  }
  ci <- c(NA_real_, NA_real_)
  if (is.finite(se)) ci <- as.numeric(est + c(-1, 1) * crit * se)
  data.frame(
    estimate = as.numeric(est),
    std.error = se,
    conf.low = ci[1],
    conf.high = ci[2],
    converged = isTRUE(x$converged),
    trimmed_fraction = diagnostics$trimmed_fraction %||% NA_real_,
    variance_method = inference$variance_method,
    solver_jacobian = diagnostics$solver_jacobian %||% NA_character_,
    jacobian_source = diagnostics$jacobian_source %||% NA_character_,
    jacobian_condition_number = diagnostics$jacobian_condition_number %||% NA_real_,
    jacobian_rel_diff = diagnostics$jacobian_rel_diff %||% NA_real_,
    max_equation_residual = diagnostics$max_equation_residual %||% NA_real_,
    min_denominator = diagnostics$min_denominator %||% NA_real_,
    fraction_small_denominators = diagnostics$fraction_small_denominators %||% NA_real_,
    used_pseudoinverse = inference$used_pseudoinverse,
    nobs = sample$n_total,
    nobs_resp = sample$n_respondents,
    is_survey = isTRUE(sample$is_survey),
    check.names = FALSE
  )
}

#' Base plotting for NMAR results
#' @description Quick base plots for weights, fitted probabilities, constraints or diagnostics.
#' @param x An object of class `nmar_result`.
#' @param which Which plot: one of `"weights"`, `"fitted"`, `"constraints"`, `"diagnostics"`.
#' @param ... Ignored.
#' @export
plot.nmar_result <- function(x, which = c("weights", "fitted", "constraints", "diagnostics"), ...) {
  which <- match.arg(which)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  weights_info <- nmar_result_get_weights_info(x)
  diagnostics <- nmar_result_get_diagnostics(x)
  if (which == "weights") {
    w <- as.numeric(weights_info$values)
    w <- w[is.finite(w)]
    if (length(w) < 2) {
      message("No weights available to plot.")
      return(invisible(x))
    }
    br <- max(2L, ceiling(sqrt(length(w))))
    graphics::hist(w, breaks = br, main = "NMAR weights (respondents)", xlab = "weight", col = "gray")
    graphics::abline(v = max(w), col = "firebrick", lty = 2)
    tr <- weights_info$trimmed_fraction
    if (length(tr) == 1 && is.finite(tr)) graphics::mtext(sprintf("trimmed_fraction = %.2f%%", 100 * tr), side = 3, cex = 0.8)
  } else if (which == "fitted") {
    fv <- as.numeric(fitted(x))
    fv <- fv[is.finite(fv)]
    if (length(fv) < 2) {
      message("No fitted values available to plot.")
      return(invisible(x))
    }
    br <- max(2L, ceiling(sqrt(length(fv))))
    graphics::hist(fv, breaks = br, main = "Fitted response probabilities (respondents)", xlab = "p_hat", col = "gray")
  } else if (which == "constraints") {
    vals <- c(eqW = diagnostics$constraint_sum_W %||% NA_real_)
    if (!is.null(diagnostics$constraint_sum_aux) && length(diagnostics$constraint_sum_aux) > 0) vals <- c(vals, diagnostics$constraint_sum_aux)
    graphics::barplot(as.numeric(vals), names.arg = names(vals), main = "Constraint sums at solution", ylab = "sum", las = 2)
    graphics::abline(h = 0, col = "darkgray")
  } else if (which == "diagnostics") {
    graphics::plot.new()
    txt <- c(
      sprintf("converged: %s", isTRUE(x$converged)),
      sprintf("variance_method: %s", (nmar_result_get_inference(x)$variance_method)),
      sprintf("Jacobian source: %s  cond: %s", diagnostics$jacobian_source %||% NA_character_, format(diagnostics$jacobian_condition_number %||% NA_real_, digits = 3)),
      sprintf("Jacobian rel diff: %s", format(diagnostics$jacobian_rel_diff %||% NA_real_, digits = 3)),
      sprintf("Max eq residual: %s", format(diagnostics$max_equation_residual %||% NA_real_, digits = 3)),
      sprintf("Min denominator: %s (frac<1e-6: %.2f%%)", format(diagnostics$min_denominator %||% NA_real_, digits = 3), 100 * (diagnostics$fraction_small_denominators %||% 0)),
      sprintf("Trimmed fraction: %.2f%%", 100 * (weights_info$trimmed_fraction %||% 0)),
      sprintf("used_pseudoinverse: %s", isTRUE(nmar_result_get_inference(x)$used_pseudoinverse))
    )
    graphics::text(0.02, 0.95, adj = c(0, 1), labels = paste(txt, collapse = "\n"), cex = 0.9)
  }
  invisible(x)
}

#' ggplot2 autoplot generic
#' @description Generic for autoplot; methods provide plotting for NMAR results.
#' @param object An object.
#' @param ... Passed to methods.
#' @export
autoplot <- function(object, ...) UseMethod("autoplot")

#' Default ggplot2 autoplot for NMAR results
#' @description Quick ggplot2 visualizations for result objects.
#' @param object An object of class `nmar_result` or a subclass.
#' @param type One of "weights", "fitted", or "constraints".
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @export
autoplot.nmar_result <- function(object, type = c("weights", "fitted", "constraints"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for autoplot.nmar_result", call. = FALSE)
  diagnostics <- nmar_result_get_diagnostics(object)
  if (type == "weights") {
    w <- stats::weights(object)
    df <- data.frame(w = as.numeric(w))
    ggplot2::ggplot(df, ggplot2::aes(x = w)) +
      ggplot2::geom_histogram(color = "white", fill = "gray") +
      ggplot2::labs(title = "NMAR weights", x = "weight")
  } else if (type == "fitted") {
    fv <- stats::fitted(object)
    df <- data.frame(p = as.numeric(fv))
    ggplot2::ggplot(df, ggplot2::aes(x = p)) +
      ggplot2::geom_histogram(color = "white", fill = "gray") +
      ggplot2::labs(title = "Fitted response probabilities", x = "p_hat")
  } else {
    vals <- c(eqW = diagnostics$constraint_sum_W %||% NA_real_)
    if (!is.null(diagnostics$constraint_sum_aux) && length(diagnostics$constraint_sum_aux) > 0) vals <- c(vals, diagnostics$constraint_sum_aux)
    df <- data.frame(term = names(vals), value = as.numeric(vals))
    ggplot2::ggplot(df, ggplot2::aes(x = term, y = value)) +
      ggplot2::geom_col(fill = "gray") +
      ggplot2::geom_hline(yintercept = 0, color = "darkgray") +
      ggplot2::labs(title = "Constraint sums", x = NULL, y = "sum") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
}

#' Default coefficients for NMAR results
#' @description Returns response-model coefficients if available.
#' @param object An `nmar_result` object.
#' @param ... Ignored.
#' @return A named numeric vector or `NULL`.
#' @export
coef.nmar_result <- function(object, ...) {
  nmar_result_get_model(object)$coefficients
}

#' Default fitted values for NMAR results
#' @description Returns fitted response probabilities if available.
#' @param object An `nmar_result` object.
#' @param ... Ignored.
#' @return A numeric vector (possibly length 0).
#' @export
fitted.nmar_result <- function(object, ...) {
  fv <- object$extra$fitted_values %||% object$fitted_values
  if (is.null(fv) || length(fv) == 0) return(numeric(0))
  as.numeric(fv)
}

#' Default weights for NMAR results
#' @description Returns respondent weights with trimmed_fraction attribute when present.
#' @param object An `nmar_result` object.
#' @param ... Ignored.
#' @return A numeric vector (possibly length 0); attribute `trimmed_fraction` may be set.
#' @export
weights.nmar_result <- function(object, ...) {
  info <- nmar_result_get_weights_info(object)
  w <- info$values
  if (is.null(w)) return(numeric(0))
  w <- as.numeric(w)
  attr(w, "trimmed_fraction") <- info$trimmed_fraction
  w
}

#' Default formula for NMAR results
#' @description Returns the estimation formula if available.
#' @param x An `nmar_result` object.
#' @param ... Ignored.
#' @return A formula or `NULL`.
#' @export
formula.nmar_result <- function(x, ...) {
  x$meta$formula %||% NULL
}
