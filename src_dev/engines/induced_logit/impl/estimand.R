#' Compute induced-logit estimand and diagnostics
#'
#' Shared implementation of Eq. (13) and related diagnostics for IID and survey
#' backends.
#'
#' @keywords internal
#' @noRd
il_compute_estimand <- function(spec, y_vec, mu_hat, gamma_hat_paper, glm_intercept, weights_full = NULL) {
  y_vec <- as.numeric(y_vec)
  mu_hat <- as.numeric(mu_hat)
  if (length(y_vec) != length(mu_hat)) stop("Internal error: y_vec and mu_hat length mismatch.", call. = FALSE)

  eps_hat <- as.numeric(y_vec[spec$respondent_mask] - mu_hat[spec$respondent_mask])
  if (anyNA(eps_hat) || any(!is.finite(eps_hat))) {
    stop("Residuals contain NA/Inf. Check outcome regression fit.", call. = FALSE)
  }

  glm_intercept <- as.numeric(glm_intercept)
  alpha_hat_paper <- if (length(glm_intercept) == 1L && is.finite(glm_intercept)) -glm_intercept else NA_real_

  weighted_mean <- function(x, weights = NULL) {
    x <- as.numeric(x)
    if (is.null(weights)) return(mean(x))
    weights <- as.numeric(weights)
    sum_w <- sum(weights)
    if (!is.finite(sum_w) || sum_w <= 0) stop("`weights_full` must sum to a positive number.", call. = FALSE)
    sum(weights * x) / sum_w
  }

  if (!is.null(weights_full)) {
    weights_full <- as.numeric(weights_full)
    if (length(weights_full) != length(mu_hat)) stop("`weights_full` must have the same length as mu_hat.", call. = FALSE)
    if (any(!is.finite(weights_full))) stop("`weights_full` must be finite.", call. = FALSE)
    if (any(weights_full < 0)) stop("`weights_full` must be nonnegative.", call. = FALSE)
  }

  r_vec <- as.integer(spec$respondent_mask)
  eta_hat <- weighted_mean(r_vec, weights = weights_full)
  mu_bar <- weighted_mean(mu_hat, weights = weights_full)

  w_resp <- if (is.null(weights_full)) NULL else weights_full[spec$respondent_mask]
  m2_over_m1 <- il_m2_over_m1_ratio_weighted(eps_hat = eps_hat, gamma = gamma_hat_paper, weights = w_resp)
  log_m1_hat <- il_log_m1_hat_weighted(eps_hat = eps_hat, gamma = gamma_hat_paper, weights = w_resp)
  alpha0_hat_paper <- if (is.finite(alpha_hat_paper) && is.finite(log_m1_hat)) alpha_hat_paper - log_m1_hat else NA_real_

  list(
    tau_hat = as.numeric(mu_bar + (1 - eta_hat) * m2_over_m1),
    eta_hat = as.numeric(eta_hat),
    mu_bar = as.numeric(mu_bar),
    m2_over_m1 = as.numeric(m2_over_m1),
    log_m1_hat = as.numeric(log_m1_hat),
    alpha_hat_paper = as.numeric(alpha_hat_paper),
    alpha0_hat_paper = as.numeric(alpha0_hat_paper)
  )
}
