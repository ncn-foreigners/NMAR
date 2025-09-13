## Shared scaling infrastructure for NMAR engines
##
## This module centralizes feature scaling and parameter unscaling so all NMAR
## estimators can share consistent, numerically stable behavior.
##
## Goals
## - Provide an engine‑agnostic API for standardizing design matrices and
##   auxiliary moments before solving.
## - Return a minimal “recipe” (per‑column mean/sd) that allows engines to
##   unscale coefficients and their covariance matrices after solving.
##
## Inputs/outputs (at a glance)
## - Inputs: `Z_un` (response model matrix incl. intercept), optional `X_un`
##   (auxiliary model matrix, no intercept), optional named `mu_x_un` (auxiliary
##   means on original scale), and a logical `standardize` flag.
## - Outputs: `Z`, `X`, `mu_x` on the scaled space, plus an
##   `nmar_scaling_recipe` used later for unscaling.
##
## Integration pattern
## - Before solving: call `validate_and_apply_nmar_scaling()` (engine-friendly) or
##   `prepare_nmar_scaling()` (low-level) to obtain scaled matrices and recipe.
## - Solve in the scaled space.
## - After solving: call `unscale_coefficients(beta_scaled, vcov_scaled, recipe)`
##   to report coefficients and vcov on the original (unscaled) scale.
## - Store `nmar_scaling_recipe` in results for diagnostics and reproducibility.
##
## Notes
## - We intentionally avoid engine‑specific assumptions; engines with non‑linear
##   parameterizations can build on top of this. The intercept column is never
##   scaled. Columns with ~0 variance are assigned sd=1 to avoid blow‑ups.
## - Engines may build the recipe on respondents‑only matrices (as EL does) or
##   on full/design matrices, depending on their estimation strategy.

#' @title Shared scaling for NMAR estimators (developer overview)
#' @name nmar_scaling_infra
#' @keywords internal
#' @noRd
NULL

new_nmar_scaling_recipe <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = "nmar_scaling_recipe")
}

validate_nmar_scaling_recipe <- function(x) {
  if (!all(vapply(x, is.list, logical(1)))) stop("Each element of a nmar_scaling_recipe must be a list.", call. = FALSE)
  has_req_fields <- vapply(x, function(comp) all(c("mean", "sd") %in% names(comp)), logical(1))
  if (!all(has_req_fields)) stop("Each element of a nmar_scaling_recipe must contain 'mean' and 'sd'.", call. = FALSE)
  x
}

#' Build a scaling recipe from one or more design matrices
#' @param ... one or more matrices with named columns
#' @param intercept_col name of the intercept column that should not be scaled
#' @keywords internal
create_nmar_scaling_recipe <- function(..., intercept_col = "(Intercept)") {
  matrices <- list(...)
  recipe <- list()
  for (mat in matrices) {
    for (col_name in colnames(mat)) {
      if (col_name == intercept_col || col_name %in% names(recipe)) next
      col_data <- mat[, col_name]
      col_mean <- mean(col_data, na.rm = TRUE)
      col_sd <- sd(col_data, na.rm = TRUE)
      if (is.na(col_sd) || col_sd < .Machine$double.eps) {
        recipe[[col_name]] <- list(mean = col_mean, sd = 1)
      } else {
        recipe[[col_name]] <- list(mean = col_mean, sd = col_sd)
      }
    }
  }
  validate_nmar_scaling_recipe(new_nmar_scaling_recipe(recipe))
}

#' Apply scaling to a matrix using a recipe
#' @keywords internal
#' @param matrix_to_scale a numeric matrix with column names present in `recipe`.
#' @param recipe an object of class `nmar_scaling_recipe`.
#' @return a matrix with each named column centered and scaled using the recipe.
apply_nmar_scaling <- function(matrix_to_scale, recipe) {
  scaled_matrix <- matrix_to_scale
  for (col_name in colnames(scaled_matrix)) {
    if (col_name %in% names(recipe)) {
      params <- recipe[[col_name]]
      scaled_matrix[, col_name] <- (scaled_matrix[, col_name] - params$mean) / params$sd
    }
  }
  scaled_matrix
}

#' Prepare scaled matrices and moments (low-level)
#' @keywords internal
#' @param Z_un response model matrix (with intercept column).
#' @param X_un auxiliary model matrix (no intercept), or NULL.
#' @param mu_x_un named numeric vector of auxiliary means on the original scale
#'   (names must match `colnames(X_un)`), or NULL.
#' @param standardize logical; apply standardization if TRUE.
#' @return a list with `Z`, `X`, `mu_x`, and `recipe`.
prepare_nmar_scaling <- function(Z_un, X_un, mu_x_un, standardize) {
  if (!standardize) {
    return(list(Z = Z_un, X = X_un, mu_x = if (is.null(mu_x_un)) numeric(0) else mu_x_un, recipe = NULL))
  }
  if (!is.null(X_un)) {
    if (!setequal(colnames(X_un), names(mu_x_un))) stop("Names of `auxiliary_means` do not match RHS variables.", call. = FALSE)
    mu_x <- mu_x_un[colnames(X_un)]
  } else {
    mu_x <- numeric(0)
  }
  recipe <- create_nmar_scaling_recipe(Z_un, if (is.null(X_un)) matrix(nrow = nrow(Z_un), ncol = 0) else X_un)
  Z <- apply_nmar_scaling(Z_un, recipe)
  X <- if (is.null(X_un)) X_un else apply_nmar_scaling(X_un, recipe)
  if (length(mu_x)) {
    mu_x <- vapply(names(mu_x), function(n) (mu_x[[n]] - recipe[[n]]$mean) / recipe[[n]]$sd, numeric(1))
  }
  list(Z = Z, X = X, mu_x = mu_x, recipe = recipe)
}

#' Validate and apply scaling (engine-friendly)
#' @keywords internal
#' @param standardize logical; apply standardization if TRUE.
#' @param has_aux logical; whether the engine uses auxiliary constraints.
#' @param response_model_matrix_unscaled response model matrix (with intercept).
#' @param auxiliary_matrix_unscaled auxiliary matrix (no intercept) or an empty matrix.
#' @param mu_x_unscaled named auxiliary means on original scale, or NULL.
#' @return a list with `nmar_scaling_recipe`, `response_model_matrix_scaled`,
#'   `auxiliary_matrix_scaled`, and `mu_x_scaled`.
validate_and_apply_nmar_scaling <- function(standardize, has_aux, response_model_matrix_unscaled,
                                            auxiliary_matrix_unscaled, mu_x_unscaled) {
  nmar_scaling_recipe <- NULL
  if (standardize) {
    if (has_aux) {
      if (!setequal(colnames(auxiliary_matrix_unscaled), names(mu_x_unscaled))) {
        stop("Names of `auxiliary_means` do not match the variables specified on the RHS of the formula.")
      }
      mu_x_unscaled <- mu_x_unscaled[colnames(auxiliary_matrix_unscaled)]
    }
    nmar_scaling_recipe <- create_nmar_scaling_recipe(response_model_matrix_unscaled, auxiliary_matrix_unscaled)
    response_model_matrix_scaled <- apply_nmar_scaling(response_model_matrix_unscaled, nmar_scaling_recipe)
    auxiliary_matrix_scaled <- apply_nmar_scaling(auxiliary_matrix_unscaled, nmar_scaling_recipe)
    mu_x_scaled <- if (has_aux) {
      vapply(names(mu_x_unscaled), function(n) (mu_x_unscaled[n] - nmar_scaling_recipe[[n]]$mean) / nmar_scaling_recipe[[n]]$sd, numeric(1))
    } else {
      numeric(0)
    }
  } else {
    response_model_matrix_scaled <- response_model_matrix_unscaled
    auxiliary_matrix_scaled <- auxiliary_matrix_unscaled
    mu_x_scaled <- if (is.null(mu_x_unscaled)) numeric(0) else mu_x_unscaled
  }
  list(
    nmar_scaling_recipe = nmar_scaling_recipe,
    response_model_matrix_scaled = response_model_matrix_scaled,
    auxiliary_matrix_scaled = auxiliary_matrix_scaled,
    mu_x_scaled = mu_x_scaled
  )
}

#' Unscale regression coefficients and covariance
#' @keywords internal
#' @param scaled_coeffs named numeric vector of coefficients estimated on the scaled space.
#' @param scaled_vcov covariance matrix of `scaled_coeffs`.
#' @param recipe `nmar_scaling_recipe` produced when scaling was applied.
#' @return a list with unscaled `coefficients` and `vcov`.
unscale_coefficients <- function(scaled_coeffs, scaled_vcov, recipe) {
  n_params <- length(scaled_coeffs)
  param_names <- names(scaled_coeffs)
  D <- diag(n_params)
  colnames(D) <- rownames(D) <- param_names
  intercept_name <- "(Intercept)"
  for (param_name in param_names) {
    if (param_name == intercept_name || !(param_name %in% names(recipe))) next
    rec <- recipe[[param_name]]
    mu <- rec$mean
    sigma <- rec$sd
    if (!is.null(sigma) && !is.na(sigma) && sigma != 0) {
      D[param_name, param_name] <- 1 / sigma
      if (intercept_name %in% param_names) D[intercept_name, param_name] <- -mu / sigma
    }
  }
  unscaled_coeffs <- drop(D %*% scaled_coeffs)
  unscaled_vcov <- D %*% scaled_vcov %*% t(D)
  list(coefficients = unscaled_coeffs, vcov = unscaled_vcov)
}
