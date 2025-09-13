#' S3 surface for base `nmar_result`
#'
#' @title Parent S3 surface for NMAR results
#' @description Methods that apply to the parent `nmar_result` class and are not
#'   specific to a particular engine (e.g., EL). Engines return a child class
#'   (e.g., `nmar_result_el`) that inherits from `nmar_result` and may override
#'   or extend behavior.
#' @details
#'   Result objects are expected to include at least:
#'   - `y_hat`: the primary estimand (typically population mean).
#'   - `se`: standard error for `y_hat`.
#'   - `data_info`: a list created by `new_nmar_data_info()` with fields such as
#'     `outcome_var`, `is_survey`, and (optionally) `design`.
#'   - `diagnostics`: a list created by `new_nmar_diagnostics()` (solver/Jacobian
#'     diagnostics, constraint sums, etc.).
#'   Optionally, engines may include:
#'   - `coefficients$response_model` and `vcov` for the response model; these are
#'     consumed by `tidy()` for reporting.
#'
#'   This separation keeps shared S3 surface concise and allows new engines to
#'   reuse these methods with minimal glue. To add a new NMAR estimator, return a
#'   child class `nmar_result_*` with these fields and engine‑specific S3 where
#'   appropriate (e.g., `weights()`, `fitted()`).
#'
#' @name nmar_result_s3_parent
#' @keywords internal
#' @import generics
#' @importFrom generics tidy glance
#' @importFrom stats fitted weights
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
  if (!is.null(x$est_mean)) {
    return(x$est_mean)
  }
  if (!is.null(x$y_hat)) {
    est <- x$y_hat
    nm <- if (!is.null(x$data_info$outcome_var)) x$data_info$outcome_var else "y"
    names(est) <- nm
    return(est)
  }
  NA_real_
}

#' Variance-covariance for base NMAR results
#' @param object An object of class `nmar_result`.
#' @param ... Ignored.
#' @export
vcov.nmar_result <- function(object, ...) {
  se <- object$se %||% object$est_var %||% NA_real_
  if (length(se) == 1 && is.finite(se)) {
    mat <- matrix(as.numeric(se)^2, 1, 1)
  } else {
    mat <- matrix(NA_real_, 1, 1)
  }
  nm <- if (!is.null(object$data_info$outcome_var)) object$data_info$outcome_var else "y"
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
  se <- object$se %||% NA_real_
  if (!is.finite(se)) {
    nm <- if (!is.null(object$data_info$outcome_var)) object$data_info$outcome_var else "y"
    return(matrix(c(NA, NA), nrow = 1, dimnames = list(nm, c("2.5 %", "97.5 %"))))
  }
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)
  ci0 <- object$est_mean %||% object$y_hat
  ci <- as.numeric(ci0) + c(-1, 1) * crit * se
  m <- matrix(ci, nrow = 1)
  colnames(m) <- paste0(format(100 * c(alpha / 2, 1 - alpha / 2)), " %")
  rownames(m) <- if (!is.null(object$data_info$outcome_var)) object$data_info$outcome_var else "y"
  m
}

#' Tidy summary for NMAR results
#' @description Return a data frame with the primary estimate and (if available) response-model coefficients.
#' @param x An object of class `nmar_result`.
#' @param conf.level Confidence level for the primary estimate.
#' @param ... Ignored.
#' @export
tidy.nmar_result <- function(x, conf.level = 0.95, ...) {
  se <- x$se
  est <- x$y_hat
  nm <- if (!is.null(x$data_info$outcome_var)) x$data_info$outcome_var else "y"
  is_svy <- isTRUE(x$data_info$is_survey)
  if (is_svy) {
    df <- tryCatch(survey::degf(x$data_info$design), error = function(e) Inf)
    crit <- stats::qt(1 - (1 - conf.level) / 2, df = df)
  } else {
    crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  }
  ci <- c(NA_real_, NA_real_)
  if (is.finite(se)) ci <- as.numeric(est + c(-1, 1) * crit * se)
  rows <- list(data.frame(term = nm, estimate = est, std.error = se, conf.low = ci[1], conf.high = ci[2], component = "estimand", statistic = NA_real_, p.value = NA_real_, check.names = FALSE))
  if (is.list(x$coefficients) && !is.null(x$coefficients$response_model)) {
    beta <- x$coefficients$response_model
    se_beta <- rep(NA_real_, length(beta))
    stat <- pval <- rep(NA_real_, length(beta))
    if (is.matrix(x$vcov)) {
      se_beta <- sqrt(diag(x$vcov))
      stat <- beta / se_beta
      pval <- if (is_svy) 2 * stats::pt(-abs(stat), df = tryCatch(survey::degf(x$data_info$design), error = function(e) Inf)) else 2 * stats::pnorm(-abs(stat))
    }
    rows[[2]] <- data.frame(term = names(beta), estimate = as.numeric(beta), std.error = se_beta, statistic = stat, p.value = pval, conf.low = NA_real_, conf.high = NA_real_, component = "response", check.names = FALSE)
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
  se <- x$se
  est <- x$y_hat
  is_svy <- isTRUE(x$data_info$is_survey)
  if (is_svy) {
    df <- tryCatch(survey::degf(x$data_info$design), error = function(e) Inf)
    crit <- stats::qt(0.975, df = df)
  } else {
    crit <- stats::qnorm(0.975)
  }
  ci <- c(NA_real_, NA_real_)
  if (is.finite(se)) ci <- as.numeric(est + c(-1, 1) * crit * se)
  data.frame(
    estimate = est, std.error = se, conf.low = ci[1], conf.high = ci[2],
    converged = x$converged, trimmed_fraction = x$diagnostics$trimmed_fraction,
    variance_method = x$data_info$variance_method, solver_jacobian = x$diagnostics$solver_jacobian,
    jacobian_source = x$diagnostics$jacobian_source, jacobian_condition_number = x$diagnostics$jacobian_condition_number,
    jacobian_rel_diff = x$diagnostics$jacobian_rel_diff, max_equation_residual = x$diagnostics$max_equation_residual,
    min_denominator = x$diagnostics$min_denominator, fraction_small_denominators = x$diagnostics$fraction_small_denominators,
    used_pseudoinverse = isTRUE(x$diagnostics$used_pseudoinverse), nobs = x$data_info$nobs, nobs_resp = x$data_info$nobs_resp,
    is_survey = is_svy, check.names = FALSE
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
  if (which == "weights") {
    w <- as.numeric(weights(x))
    w <- w[is.finite(w)]
    if (length(w) < 2) {
      message("No weights available to plot.")
      return(invisible(x))
    }
    br <- max(2L, ceiling(sqrt(length(w))))
    graphics::hist(w, breaks = br, main = "EL weights (respondents)", xlab = "weight", col = "gray")
    graphics::abline(v = max(w), col = "firebrick", lty = 2)
    tr <- attr(w, "trimmed_fraction")
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
    vals <- c(eqW = x$diagnostics$constraint_sum_W)
    if (!is.null(x$diagnostics$constraint_sum_aux) && length(x$diagnostics$constraint_sum_aux) > 0) vals <- c(vals, x$diagnostics$constraint_sum_aux)
    graphics::barplot(as.numeric(vals), names.arg = names(vals), main = "Constraint sums at solution", ylab = "sum", las = 2)
    graphics::abline(h = 0, col = "darkgray")
  } else if (which == "diagnostics") {
    graphics::plot.new()
    js <- x$diagnostics$jacobian_source
    jc <- x$diagnostics$jacobian_condition_number
    jrd <- x$diagnostics$jacobian_rel_diff
    mer <- x$diagnostics$max_equation_residual
    md <- x$diagnostics$min_denominator
    fsmall <- x$diagnostics$fraction_small_denominators
    tr <- x$diagnostics$trimmed_fraction
    txt <- c(sprintf("converged: %s", x$converged), sprintf("variance_method: %s", x$data_info$variance_method), sprintf("Jacobian source: %s  cond: %s", js, format(jc, digits = 3)), sprintf("Jacobian rel diff: %s", format(jrd, digits = 3)), sprintf("Max eq residual: %s", format(mer, digits = 3)), sprintf("Min denominator: %s (frac<1e-6: %.2f%%)", format(md, digits = 3), 100 * fsmall), sprintf("Trimmed fraction: %.2f%%", 100 * tr), sprintf("used_pseudoinverse: %s", isTRUE(x$diagnostics$used_pseudoinverse)))
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
