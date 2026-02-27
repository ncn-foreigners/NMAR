#' Missingness-model fits for induced-logit
#'
#' Fits the induced logistic regression of `R` on `x1` and `mu_hat`.
#' When `standardize = TRUE`, `x1` columns (excluding the intercept) are
#' centered/scaled as a reparameterization for numerical stability, then
#' coefficients and the covariance matrix are transformed back to the original
#' scale.
#'
#' @keywords internal
#' @noRd
il_fit_resp_iid_core <- function(spec, mu_hat, standardize, glm_control) {
  il_assert_scalar_logical(standardize, name = "standardize")
  if (is.null(glm_control)) glm_control <- stats::glm.control()
  validator_assert_list(glm_control, name = "glm_control")

  r_vec <- spec$r_vec %||% as.integer(spec$respondent_mask)
  design <- il_build_response_design(spec = spec, mu_hat = mu_hat, standardize = isTRUE(standardize), weights = NULL)
  x_glm <- design$x_iid
  ident <- il_enforce_response_model_identifiability(il_response_model_diagnostics_matrix(x_glm))

  captured <- il_fit_binomial_glmfit(x = x_glm, y = r_vec, glm_control = glm_control)
  glm_fit <- captured$fit

  glm_warnings <- captured$warnings
  if (is.character(ident$warning) && nzchar(ident$warning)) {
    glm_warnings <- c(glm_warnings, ident$warning)
  }

  beta_scaled <- stats::coef(glm_fit)
  if (anyNA(beta_scaled)) {
    stop(
      "Induced-logit response model has NA coefficients (rank deficiency / non-identifiability / separation).",
      call. = FALSE
    )
  }

  mu_coef_name <- IL_COL_MU_HAT
  gamma_hat_paper <- il_extract_gamma_hat(beta_scaled, mu_coef_name = mu_coef_name)
  fitted_vals <- tryCatch(as.numeric(glm_fit$fitted.values), error = function(e) numeric(0))
  vc_scaled <- captured$vcov

  beta_out <- beta_scaled
  vc_out <- vc_scaled
  if (isTRUE(standardize)) {
    un <- il_unscale_response_model(
      beta_scaled = beta_scaled,
      vcov_scaled = vc_scaled,
      x1_recipe = design$x1_recipe,
      mu_coef_name = mu_coef_name
    )
    beta_out <- un$coefficients
    vc_out <- un$vcov
  }

  list(
    induced_glm = glm_fit,
    beta_glm = beta_out,
    vcov = vc_out,
    fitted = fitted_vals,
    gamma_hat_paper = as.numeric(gamma_hat_paper),
    x1_recipe = design$x1_recipe,
    warnings = unique(glm_warnings),
    identifiability = ident$diag
  )
}

#' Extract a model-based covariance matrix from a `glm.fit()` result.
#'
#' `stats::glm.fit()` returns a plain list with QR decomposition components but
#' has no `vcov()` method. For full-rank fits we invert the upper-triangular
#' `R` factor (via `chol2inv`) and apply the pivoting information to return a
#' covariance matrix aligned with coefficient names.
#'
#' @keywords internal
#' @noRd
il_vcov_from_glmfit <- function(fit, coef_names) {
  if (is.null(fit) || !is.list(fit)) return(NULL)
  if (is.null(coef_names) || length(coef_names) == 0L) return(NULL)
  qr <- fit$qr %||% NULL
  if (is.null(qr) || is.null(qr$qr)) return(NULL)

  p <- length(coef_names)
  rank <- as.integer(fit$rank %||% NA_integer_)
  if (!is.finite(rank) || rank < 1L) return(NULL)
  if (rank < p) return(NULL)

  R <- qr$qr[seq_len(rank), seq_len(rank), drop = FALSE]
  R[lower.tri(R)] <- 0
  cov_unscaled <- tryCatch(base::chol2inv(R), error = function(e) NULL)
  if (is.null(cov_unscaled)) return(NULL)

  pivot <- qr$pivot %||% seq_len(p)
  pivot <- as.integer(pivot)
  if (length(pivot) < p) return(NULL)

  vc <- matrix(0, nrow = p, ncol = p, dimnames = list(coef_names, coef_names))
  vc[pivot, pivot] <- cov_unscaled
  vc
}

#' Binomial GLM fit using `glm.fit()` with explicit convergence checks.
#'
#' This is used for the induced logistic regression step in the IID backend to
#' avoid formula parsing and to allow explicit control over standardization and
#' coefficient unscaling.
#'
#' @keywords internal
#' @noRd
il_fit_binomial_glmfit <- function(x, y, glm_control = stats::glm.control()) {
  if (!is.matrix(x)) stop("Internal error: `x` must be a matrix.", call. = FALSE)
  if (is.null(colnames(x))) stop("Internal error: `x` must have column names.", call. = FALSE)
  y <- as.numeric(y)
  if (length(y) != nrow(x)) stop("Internal error: `y` length must match nrow(x).", call. = FALSE)
  if (anyNA(y) || any(!y %in% c(0, 1))) stop("Internal error: `y` must be 0/1 without NA.", call. = FALSE)

  if (is.null(glm_control)) glm_control <- stats::glm.control()
  validator_assert_list(glm_control, name = "glm_control")

  captured <- il_capture_warnings(
    stats::glm.fit(x = x, y = y, family = stats::binomial(), control = glm_control)
  )
  fit <- captured$value
  if (!isTRUE(fit$converged)) {
    stop(
      "Induced-logit response model did not converge. ",
      "Try increasing iterations via `control = list(maxit = ...)` or standardize = TRUE.",
      call. = FALSE
    )
  }
  beta <- stats::coef(fit)
  vc <- il_vcov_from_glmfit(fit, coef_names = names(beta))
  list(fit = fit, vcov = vc, warnings = captured$warnings)
}

#' Build induced-logit response-model predictors
#'
#' Shared builder for the response design used by IID (`glm.fit` on a matrix) and
#' survey (`svyglm` on a design whose variables include `x1` columns and `mu_hat`).
#'
#' @keywords internal
#' @noRd
il_build_response_design <- function(spec, mu_hat, standardize = FALSE, weights = NULL) {
  mf <- spec$model_frame
  n <- nrow(mf)
  mu_col <- as.numeric(mu_hat)
  if (length(mu_col) != n) stop("Internal error: mu_hat length mismatch.", call. = FALSE)

  x1_mat_full <- spec$x1_mat_full %||% stop("Internal error: missing spec$x1_mat_full.", call. = FALSE)
  if (IL_COL_MU_HAT %in% colnames(x1_mat_full)) {
    stop("Internal error: x1 model matrix already contains '", IL_COL_MU_HAT, "'.", call. = FALSE)
  }
  x1_names <- setdiff(colnames(x1_mat_full), "(Intercept)")

  x1_recipe <- NULL
  x1_mat_scaled <- x1_mat_full[, x1_names, drop = FALSE]
  if (isTRUE(standardize) && length(x1_names) > 0L) {
    x1_scale <- il_scale_matrix(x1_mat_scaled, weights = weights)
    x1_recipe <- x1_scale$recipe
    x1_mat_scaled <- x1_scale$mat
  }

  x_iid <- cbind(`(Intercept)` = rep.int(1, n), x1_mat_scaled, mu_col)
  colnames(x_iid)[ncol(x_iid)] <- IL_COL_MU_HAT

  x1_df <- if (ncol(x1_mat_scaled) > 0L) {
    as.data.frame(x1_mat_scaled, check.names = FALSE)
  } else {
    data.frame(check.names = FALSE)
  }

  list(
    x_iid = x_iid,
    x1_df = x1_df,
    x1_recipe = x1_recipe
  )
}
