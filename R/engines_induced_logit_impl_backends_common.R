#' Shared induced-logit backend helpers
#'
#' @keywords internal
#' @noRd
il_unquote_backticked_names <- function(x) {
  if (is.null(x)) return(x)
  nm <- names(x)
  if (is.null(nm)) return(x)
  names(x) <- gsub("^`|`$", "", nm)
  x
}

#' @keywords internal
#' @noRd
il_unquote_backticked_dimnames <- function(mat) {
  if (!is.matrix(mat)) return(mat)
  dn <- dimnames(mat)
  if (is.null(dn)) return(mat)
  dn[[1L]] <- if (is.null(dn[[1L]])) NULL else gsub("^`|`$", "", dn[[1L]])
  dn[[2L]] <- if (is.null(dn[[2L]])) NULL else gsub("^`|`$", "", dn[[2L]])
  dimnames(mat) <- dn
  mat
}

#' Unscales standardized missingness-model coefficients and covariance back to
#' the original predictor scale.
#'
#' @keywords internal
#' @noRd
il_unscale_response_model <- function(beta_scaled, vcov_scaled, x1_recipe, mu_coef_name = IL_COL_MU_HAT) {
  beta_scaled <- il_unquote_backticked_names(beta_scaled)
  vcov_scaled <- il_unquote_backticked_dimnames(vcov_scaled)

  recipe_aug <- if (!is.null(x1_recipe)) x1_recipe else new_nmar_scaling_recipe(list())
  recipe_aug[[mu_coef_name]] <- list(mean = 0, sd = 1)

  if (is.matrix(vcov_scaled)) {
    un <- unscale_coefficients(scaled_coeffs = beta_scaled, scaled_vcov = vcov_scaled, recipe = recipe_aug)
    return(list(coefficients = un$coefficients, vcov = un$vcov))
  }

  un <- unscale_coefficients(
    scaled_coeffs = beta_scaled,
    scaled_vcov = matrix(
      0,
      nrow = length(beta_scaled),
      ncol = length(beta_scaled),
      dimnames = list(names(beta_scaled), names(beta_scaled))
    ),
    recipe = recipe_aug
  )
  list(coefficients = un$coefficients, vcov = NULL)
}

#' @keywords internal
#' @noRd
il_assert_named_coefficients <- function(coef_vec, mm_cols, message_prefix) {
  if (is.null(coef_vec) || anyNA(coef_vec)) {
    stop(paste0(message_prefix, " has NA coefficients (rank deficiency)."), call. = FALSE)
  }
  if (!all(mm_cols %in% names(coef_vec))) {
    stop(
      sprintf("Internal error: %s coefficient names do not match model matrix columns.", message_prefix),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
il_extract_gamma_hat <- function(beta_glm, mu_coef_name = IL_COL_MU_HAT) {
  if (!mu_coef_name %in% names(beta_glm)) {
    stop("Internal error: induced-logit model is missing the mu_hat coefficient.", call. = FALSE)
  }
  gamma_hat_paper <- unname(-beta_glm[mu_coef_name])
  if (!is.finite(gamma_hat_paper)) stop("Estimated gamma is not finite.", call. = FALSE)
  as.numeric(gamma_hat_paper)
}

#' @keywords internal
#' @noRd
il_require_survey <- function() {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for induced-logit with survey designs. Please install it.", call. = FALSE)
  }
}

#' @keywords internal
#' @noRd
il_get_survey_weights <- function(design, vars) {
  w_full <- tryCatch(as.numeric(stats::weights(design)), error = function(e) NULL)
  if (is.null(w_full) || length(w_full) != nrow(vars) || any(!is.finite(w_full))) {
    stop("Survey design weights must be finite and aligned with the number of observations.", call. = FALSE)
  }
  if (any(w_full < 0)) stop("Survey design weights must be nonnegative.", call. = FALSE)
  sum_w <- sum(w_full)
  if (!is.finite(sum_w) || sum_w <= 0) stop("Survey design weights must sum to a positive number.", call. = FALSE)
  list(weights = w_full, sum_weights = sum_w)
}
