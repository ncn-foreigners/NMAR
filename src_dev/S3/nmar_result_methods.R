#' S3 surface for base `nmar_result`
#'
#' @title Parent S3 surface for NMAR results
#' @description Methods that apply to the parent `nmar_result` class and are not
#'   specific to a particular engine (e.g., EL). Engines return a child class
#'   (e.g., `nmar_result_el`) that inherits from `nmar_result` and may override
#'   or extend behavior.
#' @details
#'   Result objects expose a universal schema:
#'   - `y_hat`, `estimate_name`, `se`, `converged`.
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



#' Variance-covariance for base NMAR results
#' @param object An object of class `nmar_result`.
#' @param ... Ignored.
#' @keywords result_param
#' @export
vcov.nmar_result <- function(object, ...) {
  se <- nmar_result_get_se(object)
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
#' @keywords result_param
#' @export
confint.nmar_result <- function(object, parm, level = 0.95, ...) {
  se <- nmar_result_get_se(object)
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
#' @keywords result_view
#' @method tidy nmar_result
#' @export
tidy.nmar_result <- function(x, conf.level = 0.95, ...) {
  est <- nmar_result_get_estimate(x)
  se <- nmar_result_get_se(x)
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
    y_hat = as.numeric(est),
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
      y_hat = beta_vec,
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
#' @keywords result_view
#' @method glance nmar_result
#' @export
glance.nmar_result <- function(x, ...) {
  est <- nmar_result_get_estimate(x)
  se <- nmar_result_get_se(x)
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
    y_hat = as.numeric(est),
    std.error = se,
    conf.low = ci[1],
    conf.high = ci[2],
    converged = isTRUE(x$converged),
    trimmed_fraction = diagnostics$trimmed_fraction %||% NA_real_,
    variance_method = inference$variance_method,
    jacobian_condition_number = diagnostics$jacobian_condition_number %||% NA_real_,
    max_equation_residual = diagnostics$max_equation_residual %||% NA_real_,
    min_denominator = diagnostics$min_denominator %||% NA_real_,
    fraction_small_denominators = diagnostics$fraction_small_denominators %||% NA_real_,

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
#' @keywords internal
#' @method plot nmar_result
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
      sprintf("Jacobian cond: %s", format(diagnostics$jacobian_condition_number %||% NA_real_, digits = 3)),
      sprintf("Max eq residual: %s", format(diagnostics$max_equation_residual %||% NA_real_, digits = 3)),
      sprintf("Min denominator: %s (frac<1e-6: %.2f%%)", format(diagnostics$min_denominator %||% NA_real_, digits = 3), 100 * (diagnostics$fraction_small_denominators %||% 0)),
      sprintf("Trimmed fraction: %.2f%%", 100 * (weights_info$trimmed_fraction %||% 0)),
      NULL
    )
    txt <- txt[!vapply(txt, is.null, logical(1))]
    graphics::text(0.02, 0.95, adj = c(0, 1), labels = paste(unlist(txt), collapse = "\n"), cex = 0.9)
  }
  invisible(x)
}

#' Autoplot generic
#'
#' Generic function for creating ggplot2 plots. If ggplot2 is loaded, its
#' generic will be used; otherwise, this minimal generic enables S3 dispatch.
#' @param object An object
#' @param ... Additional arguments passed to methods
#' @export
autoplot <- function(object, ...) {
  UseMethod("autoplot")
}

#' ggplot2 autoplot method for NMAR results
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
#' @keywords result_param
#' @export
coef.nmar_result <- function(object, ...) {
  nmar_result_get_model(object)$coefficients
}

#' Default fitted values for NMAR results
#' @description Returns fitted response probabilities if available.
#' @param object An `nmar_result` object.
#' @param ... Ignored.
#' @return A numeric vector (possibly length 0).
#' @keywords result_param
#' @export
fitted.nmar_result <- function(object, ...) {
  fv <- object$extra$fitted_values %||% object$fitted_values
  if (is.null(fv) || length(fv) == 0) {
    return(numeric(0))
  }
  as.numeric(fv)
}

#' Extract Weights from NMAR Result
#'
#' @param object An object of class \code{nmar_result}
#' @param scale Character: \code{"probability"} (default) or \code{"population"}.
#'   \describe{
#'     \item{\code{"probability"}}{
#'       Returns p_i = w_tilde_i / sum_j w_tilde_j where sum_i p_i = 1 (exact).
#'       This is the paper's canonical form (QLS 2002, Eq. 11).
#'       Use for computing means: y_bar = sum_i p_i * y_i
#'     }
#'     \item{\code{"population"}}{
#'       Returns w_i = N_pop * p_i where sum_i w_i = N_pop (exact).
#'       This follows survey package conventions.
#'       Use for computing totals: T_hat = sum_i w_i * y_i = N_pop * y_bar
#'     }
#'   }
#' @param ... Additional arguments (ignored)
#'
#' @details
#' The empirical likelihood estimator computes unnormalized masses
#' w_tilde_i = d_i / D_i that satisfy the constraint sum_i w_tilde_i = sum_i d_i (without trimming).
#' This method provides two standardized representations:
#'
#' \strong{Mathematical guarantees} (hold even with trimming):
#' \itemize{
#'   \item \code{sum(weights(object, scale = "probability")) = 1} (within machine precision)
#'   \item \code{sum(weights(object, scale = "population")) = N_pop} (within machine precision)
#'   \item \code{weights(object, "population") = N_pop * weights(object, "probability")} (exact)
#' }
#'
#' \strong{Trimming effects}:
#' When \code{trim_cap < Inf} and trimming is active, the normalization
#' identity sum_i w_tilde_i = sum_i d_i is violated. However, this method still returns
#' weights with correct sums by using the formula:
#' \deqn{w_i = N_pop * w_tilde_i / sum_j w_tilde_j}
#'
#' @return Numeric vector of weights with length equal to number of respondents
#'
#' @references
#' Qin, J., Leung, D., & Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. \emph{Journal of the
#' American Statistical Association}, 97(457), 193-200.
#'
#' @examples
#' \dontrun{
#' res <- nmar(y_miss ~ x, data = df, engine = el_engine())
#'
#' # Probability weights (default): sum to 1
#' w_prob <- weights(res)
#' sum(w_prob) # Exactly 1.0
#'
#' # Population weights: sum to N_pop
#' w_pop <- weights(res, scale = "population")
#' sum(w_pop) # Exactly nrow(df)
#'
#' # Relationship (exact):
#' all.equal(w_pop, nrow(df) * w_prob)
#' }
#'
#' @keywords result_param
#' @export
weights.nmar_result <- function(object,
                                scale = c("probability", "population"),
                                ...) {
  scale <- match.arg(scale)

# Validate input
  if (!inherits(object, "nmar_result")) {
    stop("'object' must be of class 'nmar_result'", call. = FALSE)
  }

# Do not warn on nonconvergence here to keep downstream
# code (including tests) quiet. Convergence is reported
# in print/summary and stored in diagnostics; callers
# can inspect object$converged and diagnostics as needed.

# Extract stored unnormalized masses (single source of truth)
  info <- nmar_result_get_weights_info(object)
  w_unnorm <- info$values

  if (is.null(w_unnorm) || length(w_unnorm) == 0) {
    return(numeric(0))
  }

  w_unnorm <- as.numeric(w_unnorm)

  if (any(!is.finite(w_unnorm))) {
    stop("Non-finite weights detected in result object", call. = FALSE)
  }

# Compute sum
  sum_w <- sum(w_unnorm)

  if (sum_w <= 0) {
    stop(sprintf("Invalid weight sum: %g", sum_w), call. = FALSE)
  }

# Scale according to user request
  if (scale == "probability") {
# Paper's canonical form: p_i = w_tilde_i / sum_j w_tilde_j
# GUARANTEE: sum_i p_i = 1 (within machine precision)
    weights_out <- w_unnorm / sum_w

  } else { # scale == "population"
# Survey convention: w_i = N_pop * p_i
# GUARANTEE: sum_i w_i = N_pop (within machine precision)

    sample <- nmar_result_get_sample(object)
    N_pop <- sample$n_total

# Defensive check
    if (is.null(N_pop) || !is.finite(N_pop) || N_pop <= 0) {
      stop(sprintf(
        "Invalid or missing N_pop in result object: %s",
        if (is.null(N_pop)) "NULL" else as.character(N_pop)
      ), call. = FALSE)
    }

# CRITICAL FIX: This formula guarantees sum = N_pop even with trimming
    weights_out <- N_pop * w_unnorm / sum_w
  }

# Add informative attributes
  attr(weights_out, "scale") <- scale
  attr(weights_out, "trimmed_fraction") <- info$trimmed_fraction
  if (scale == "population") {
    sample <- nmar_result_get_sample(object)
    attr(weights_out, "N_pop") <- sample$n_total
  }

  return(weights_out)
}

#' Default formula for NMAR results
#' @description Returns the estimation formula if available.
#' @param x An `nmar_result` object.
#' @param ... Ignored.
#' @return A formula or `NULL`.
#' @keywords result_param
#' @export
formula.nmar_result <- function(x, ...) {
  x$meta$formula %||% NULL
}


#' Standard error generic
#'
#' Generic function for extracting standard errors from model objects.
#' @param object An object
#' @param ... Additional arguments passed to methods
#' @return Numeric value or vector of standard errors
#' @export
se <- function(object, ...) {
  UseMethod("se")
}

#' Extract standard error for NMAR results
#'
#' Returns the standard error of the primary mean estimate.
#' @param object An `nmar_result` or subclass.
#' @param ... Ignored.
#' @return Numeric scalar.
#' @keywords result_param
#' @exportS3Method se nmar_result
se.nmar_result <- function(object, ...) nmar_result_get_se(object)





#' Coefficient table for summary objects
#'
#' Returns a coefficients table (Estimate, Std. Error, statistic, p-value)
#' from a `summary_nmar_result*` object when response-model coefficients and a
#' variance matrix are available. If the summary does not carry response-model
#' coefficients, returns `NULL`.
#'
#' The statistic column is labelled "t value" when finite degrees of freedom
#' are available (e.g., survey designs); otherwise, it is labelled "z value".
#'
#' @param object An object of class `summary_nmar_result` (or subclass).
#' @param ... Ignored.
#' @return A data.frame with rows named by coefficient, or `NULL` if not available.
#' @keywords result_view
#' @export
coef.summary_nmar_result <- function(object, ...) {
  beta <- object$response_model
  if (is.null(beta) || length(beta) == 0) return(NULL)
  beta_names <- names(beta)
  beta_vec <- as.numeric(beta)
  se <- rep(NA_real_, length(beta_vec))
  if (!is.null(object$response_vcov) && is.matrix(object$response_vcov)) {
    se <- sqrt(diag(object$response_vcov))
  }
  stat <- beta_vec / se
  df <- object$df %||% NA_real_
  stat_label <- if (is.finite(df)) "t value" else "z value"
  pval <- if (is.finite(df)) 2 * stats::pt(-abs(stat), df = df) else 2 * stats::pnorm(-abs(stat))
  if (is.null(beta_names) || length(beta_names) != length(beta_vec) || anyNA(beta_names)) {
    beta_names <- paste0("coef", seq_along(beta_vec))
  }
  tab <- data.frame(y_hat = beta_vec, `Std. Error` = se, check.names = FALSE)
  tab[[stat_label]] <- stat
  p_label <- if (is.finite(df)) "Pr(>|t|)" else "Pr(>|z|)"
  tab[[p_label]] <- pval
  rownames(tab) <- beta_names
  tab
}

#' Confidence intervals for coefficient table (summary objects)
#'
#' Returns Wald-style confidence intervals for response-model coefficients from
#' a `summary_nmar_result*` object. Uses t-quantiles when finite degrees of
#' freedom are available, otherwise normal quantiles.
#'
#' @param object An object of class `summary_nmar_result` (or subclass).
#' @param parm A specification of which coefficients are to be given confidence intervals,
#'   either a vector of names or a vector of indices; by default, all coefficients are considered.
#' @param level The confidence level required.
#' @param ... Ignored.
#' @return A numeric matrix with columns giving lower and upper confidence limits for each parameter.
#'   Row names correspond to coefficient names. Returns `NULL` if coefficients are unavailable.
#' @keywords result_view
#' @export
confint.summary_nmar_result <- function(object, parm, level = 0.95, ...) {
  beta <- object$response_model
  if (is.null(beta) || length(beta) == 0) return(NULL)
  beta_names <- names(beta)
  beta_vec <- as.numeric(beta)
  vc <- object$response_vcov
  if (is.null(vc) || !is.matrix(vc)) {
# No variance; return NA intervals with proper shape
    idx <- seq_along(beta_vec)
    if (!missing(parm)) {
      if (is.character(parm)) idx <- match(parm, beta_names, nomatch = 0L) else idx <- as.integer(parm)
      idx <- idx[idx > 0 & idx <= length(beta_vec)]
      if (!length(idx)) return(NULL)
    }
    out <- cbind(`2.5 %` = rep(NA_real_, length(idx)), `97.5 %` = rep(NA_real_, length(idx)))
    rownames(out) <- if (!is.null(beta_names)) beta_names[idx] else paste0("coef", idx)
    return(out)
  }
  se <- sqrt(diag(vc))
# Subset if requested
  idx <- seq_along(beta_vec)
  if (!missing(parm)) {
    if (is.character(parm)) idx <- match(parm, beta_names, nomatch = 0L) else idx <- as.integer(parm)
    idx <- idx[idx > 0 & idx <= length(beta_vec)]
    if (!length(idx)) return(NULL)
  }
  alpha <- 1 - level
  df <- object$df %||% NA_real_
  crit <- if (is.finite(df)) stats::qt(1 - alpha / 2, df = df) else stats::qnorm(1 - alpha / 2)
  lo <- beta_vec[idx] - crit * se[idx]
  hi <- beta_vec[idx] + crit * se[idx]
  out <- cbind(lo, hi)
  colnames(out) <- paste0(format(100 * c(alpha / 2, 1 - alpha / 2)), " %")
  rownames(out) <- if (!is.null(beta_names)) beta_names[idx] else paste0("coef", idx)
  out
}
