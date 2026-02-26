#' Outcome-regression fits for induced-logit
#'
#' @keywords internal
#' @noRd
new_nmar_il_mu_fit <- function(mu_rhs, coef, recipe = NULL, fit = NULL) {
  if (!inherits(mu_rhs, "formula") || length(mu_rhs) != 2L) {
    stop("Internal error: `mu_rhs` must be a one-sided formula.", call. = FALSE)
  }
  coef <- il_unquote_backticked_names(coef)
  nm <- names(coef)
  coef <- as.numeric(coef)
  names(coef) <- nm
  if (is.null(nm) || length(nm) != length(coef)) {
    stop("Internal error: mu-fit coefficients must be a named numeric vector.", call. = FALSE)
  }
  out <- list(mu_rhs = mu_rhs, coef = coef, recipe = recipe, fit = fit)
  class(out) <- "nmar_il_mu_fit"
  out
}

#' @keywords internal
#' @noRd
#' @exportS3Method predict nmar_il_mu_fit
predict.nmar_il_mu_fit <- function(object, newdata, ...) {
  if (missing(newdata)) {
    stop("`newdata` is required for predictions from internal induced-logit mu fits.", call. = FALSE)
  }
  if (!is.data.frame(newdata)) stop("`newdata` must be a data.frame.", call. = FALSE)

  mf <- stats::model.frame(object$mu_rhs, data = newdata, na.action = stats::na.pass, drop.unused.levels = TRUE)
  mm <- stats::model.matrix(object$mu_rhs, data = mf)
  if (!is.null(object$recipe)) {
    mm <- apply_nmar_scaling(mm, object$recipe)
  }

  coef <- object$coef
  mm_cols <- names(coef)
  if (!all(mm_cols %in% colnames(mm))) {
    stop(
      "Outcome regression prediction failed: model matrix columns did not match stored coefficients.",
      call. = FALSE
    )
  }
  as.numeric(mm[, mm_cols, drop = FALSE] %*% coef[mm_cols])
}

#' @keywords internal
#' @noRd
#' @exportS3Method coef nmar_il_mu_fit
coef.nmar_il_mu_fit <- function(object, ...) {
  object$coef
}

#' @keywords internal
#' @noRd
il_fit_mu_scaled <- function(spec) {
  mf <- spec$model_frame
  mu_mat_full <- spec$mu_mat_full %||% stop("Internal error: missing spec$mu_mat_full.", call. = FALSE)
  mu_scale <- il_scale_mu_matrix(mu_mat_full, has_intercept = spec$mu_intercept %||% TRUE)

  d_resp_idx <- which(spec$respondent_mask)
  y_resp <- mf[[spec$outcome]][d_resp_idx]
  x_resp <- mu_scale$mat[d_resp_idx, , drop = FALSE]

  mu_fit <- withCallingHandlers(
    stats::lm.fit(x = x_resp, y = y_resp),
    warning = function(w) invokeRestart("muffleWarning")
  )

  if (is.null(mu_fit$coefficients) || anyNA(mu_fit$coefficients)) {
    stop("Outcome regression has NA coefficients (rank deficiency).", call. = FALSE)
  }

  coef_vec <- as.numeric(mu_fit$coefficients)
  names(coef_vec) <- colnames(x_resp) %||% names(mu_fit$coefficients)
  if (is.null(names(coef_vec))) stop("Internal error: mu-fit coefficients are missing names.", call. = FALSE)

  mu_hat <- as.numeric(mu_scale$mat %*% coef_vec[colnames(mu_scale$mat)])
  if (length(mu_hat) != nrow(mf)) stop("Internal error: mu_hat length mismatch.", call. = FALSE)
  if (anyNA(mu_hat) || any(!is.finite(mu_hat))) {
    stop("Outcome regression produced NA/Inf predictions on full data.", call. = FALSE)
  }

  list(
    mu_fit = new_nmar_il_mu_fit(mu_rhs = spec$mu_rhs, coef = coef_vec, recipe = mu_scale$recipe, fit = mu_fit),
    mu_hat = mu_hat,
    mu_recipe = mu_scale$recipe
  )
}

#' @keywords internal
#' @noRd
il_fit_mu_unscaled <- function(spec) {
  mf <- spec$model_frame
  mu_mat_full <- spec$mu_mat_full %||% stop("Internal error: missing spec$mu_mat_full.", call. = FALSE)

  d_resp_idx <- which(spec$respondent_mask)
  y_resp <- mf[[spec$outcome]][d_resp_idx]
  x_resp <- mu_mat_full[d_resp_idx, , drop = FALSE]

  mu_fit <- withCallingHandlers(
    stats::lm.fit(x = x_resp, y = y_resp),
    warning = function(w) invokeRestart("muffleWarning")
  )

  if (is.null(mu_fit$coefficients) || anyNA(mu_fit$coefficients)) {
    stop("Outcome regression has NA coefficients (rank deficiency).", call. = FALSE)
  }

  coef_vec <- as.numeric(mu_fit$coefficients)
  names(coef_vec) <- colnames(x_resp) %||% names(mu_fit$coefficients)
  if (is.null(names(coef_vec))) stop("Internal error: mu-fit coefficients are missing names.", call. = FALSE)

  mu_hat <- as.numeric(mu_mat_full %*% coef_vec[colnames(mu_mat_full)])
  if (length(mu_hat) != nrow(mf)) stop("Internal error: mu_hat length mismatch.", call. = FALSE)
  if (anyNA(mu_hat) || any(!is.finite(mu_hat))) {
    stop(
      "Outcome regression produced NA/Inf predictions on full data.\n  ",
      "This can happen if the outcome model is rank-deficient or if factor levels appear only among nonrespondents.",
      call. = FALSE
    )
  }

  list(
    mu_fit = new_nmar_il_mu_fit(
      mu_rhs = spec$mu_rhs,
      coef = coef_vec[colnames(mu_mat_full)],
      recipe = NULL,
      fit = mu_fit
    ),
    mu_hat = mu_hat,
    warnings = character()
  )
}

#' @keywords internal
#' @noRd
il_fit_mu_iid_core <- function(spec, standardize) {
  il_assert_scalar_logical(standardize, name = "standardize")

  if (isTRUE(standardize)) {
    mu_res <- il_fit_mu_scaled(spec)
    return(list(
      mu_fit = mu_res$mu_fit,
      mu_hat = mu_res$mu_hat,
      mu_recipe = mu_res$mu_recipe,
      warnings = character()
    ))
  }

  mu_res <- il_fit_mu_unscaled(spec)
  list(
    mu_fit = mu_res$mu_fit,
    mu_hat = mu_res$mu_hat,
    mu_recipe = NULL,
    warnings = mu_res$warnings %||% character()
  )
}
