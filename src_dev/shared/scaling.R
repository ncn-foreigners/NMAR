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
## - Engines may build the recipe on respondents-only matrices (as EL does) or
##   on the full design, and can request design-weighted scaling by supplying
##   weights through the API.

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

compute_weighted_stats <- function(values, weights = NULL) {
  values <- as.numeric(values)
  if (is.null(weights)) {
    mean_val <- mean(values, na.rm = TRUE)
    sd_val <- stats::sd(values, na.rm = TRUE)
    return(list(mean = mean_val, sd = sd_val))
  }
  if (length(weights) != length(values)) {
    stop("`weights` must have the same length as the data being scaled.", call. = FALSE)
  }
  weights <- as.numeric(weights)
  mask <- !is.na(values) & is.finite(values) & !is.na(weights)
  if (!any(mask)) {
    mean_val <- NA_real_
    sd_val <- NA_real_
    return(list(mean = mean_val, sd = sd_val))
  }
  w <- weights[mask]
  x <- values[mask]
  w_sum <- sum(w)
  if (w_sum <= 0 || all(w == 0)) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- stats::sd(x, na.rm = TRUE)
    return(list(mean = mean_val, sd = sd_val))
  }
  mean_val <- sum(w * x) / w_sum
  variance_num <- sum(w * (x - mean_val)^2)
  sd_val <- sqrt(variance_num / w_sum)
  list(mean = mean_val, sd = sd_val)
}

#' Build a scaling recipe from one or more design matrices
#' @param ... one or more matrices with named columns.
#' @param intercept_col Intercept column name that should remain unscaled.
#' @param weights Optional numeric vector of weights used to compute weighted
#'   means/standard deviations.
#' @param weight_mask Optional logical/ numeric mask applied to `weights` before
#'   computing moments (useful for respondents-only scaling).
#' @param tol_constant Numeric tolerance below which columns are treated as
#'   constant and left unscaled.
#' @param warn_on_constant Logical; emit a warning when a column is treated as
#'   constant.
#' @keywords internal
create_nmar_scaling_recipe <- function(..., intercept_col = "(Intercept)", weights = NULL,
                                       weight_mask = NULL, tol_constant = 1e-8,
                                       warn_on_constant = TRUE) {
  matrices <- list(...)
  recipe <- list()
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }
  if (!is.null(weight_mask)) {
    weight_mask <- as.logical(weight_mask)
    if (!is.null(weights)) {
      if (length(weight_mask) != length(weights)) {
        stop("`weight_mask` must have the same length as `weights`.", call. = FALSE)
      }
      weights <- weights * as.numeric(weight_mask)
    } else {
      weights <- as.numeric(weight_mask)
    }
  }
  for (mat in matrices) {
    for (col_name in colnames(mat)) {
      if (col_name == intercept_col || col_name %in% names(recipe)) next
      col_data <- mat[, col_name]
      stats <- compute_weighted_stats(col_data, weights)
      col_mean <- stats$mean
      col_sd <- stats$sd
      if (is.na(col_sd) || col_sd < tol_constant) {
        if (warn_on_constant) {
          warning(sprintf("Column '%s' is nearly constant under the scaling weights; leaving unscaled.", col_name), call. = FALSE)
        }
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
#' @param weights Optional numeric vector used for weighted scaling.
#' @param weight_mask Optional logical/numeric mask applied to `weights`.
#' @return a list with `Z`, `X`, `mu_x`, and `recipe`.
prepare_nmar_scaling <- function(Z_un, X_un, mu_x_un, standardize,
                                 weights = NULL, weight_mask = NULL) {
  if (!standardize) {
    return(list(Z = Z_un, X = X_un, mu_x = if (is.null(mu_x_un)) numeric(0) else mu_x_un, recipe = NULL))
  }
  if (!is.null(X_un)) {
    if (!setequal(colnames(X_un), names(mu_x_un))) stop("Names of `auxiliary_means` do not match RHS variables.", call. = FALSE)
    mu_x <- mu_x_un[colnames(X_un)]
  } else {
    mu_x <- numeric(0)
  }
  recipe <- create_nmar_scaling_recipe(Z_un, if (is.null(X_un)) matrix(nrow = nrow(Z_un), ncol = 0) else X_un,
    weights = weights, weight_mask = weight_mask
  )
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
#' @param weights Optional numeric vector used for weighted scaling.
#' @param weight_mask Optional logical/numeric mask applied to `weights`.
#' @return a list with `nmar_scaling_recipe`, `response_model_matrix_scaled`,
#'   `auxiliary_matrix_scaled`, and `mu_x_scaled`.
validate_and_apply_nmar_scaling <- function(standardize, has_aux, response_model_matrix_unscaled,
                                            auxiliary_matrix_unscaled, mu_x_unscaled,
                                            weights = NULL, weight_mask = NULL) {
  nmar_scaling_recipe <- NULL
  if (standardize) {
    if (has_aux) {
      if (!setequal(colnames(auxiliary_matrix_unscaled), names(mu_x_unscaled))) {
        stop("Names of `auxiliary_means` do not match the variables specified on the RHS of the formula.")
      }
      mu_x_unscaled <- mu_x_unscaled[colnames(auxiliary_matrix_unscaled)]
    }
    if (!is.null(weights) && length(weights) != nrow(response_model_matrix_unscaled)) {
      stop("`weights` must have the same length as the number of rows in the response matrix.", call. = FALSE)
    }
    if (!is.null(weight_mask) && length(weight_mask) != nrow(response_model_matrix_unscaled)) {
      stop("`weight_mask` must have the same length as the number of rows in the response matrix.", call. = FALSE)
    }
    nmar_scaling_recipe <- create_nmar_scaling_recipe(
      response_model_matrix_unscaled,
      auxiliary_matrix_unscaled,
      weights = weights,
      weight_mask = weight_mask
    )
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

#' Map unscaled coefficients to scaled space
#' @keywords internal
#' @param beta_unscaled named numeric vector of coefficients for the response
#'   model on the original scale, including an intercept named `"(Intercept)"`.
#' @param recipe `nmar_scaling_recipe` returned by scaling utilities.
#' @param columns character vector of column names (order) for the scaled design
#'   matrix (including intercept).
#' @return numeric vector of coefficients in the scaled space, ordered by
#'   `columns`.
scale_coefficients <- function(beta_unscaled, recipe, columns) {
  if (is.null(beta_unscaled) || length(beta_unscaled) == 0) {
    return(setNames(numeric(length(columns)), columns))
  }
  beta_unscaled <- beta_unscaled[intersect(names(beta_unscaled), columns)]
  out <- setNames(numeric(length(columns)), columns)
# Non-intercept: multiply by sd
  for (nm in columns) {
    if (nm == "(Intercept)") next
    if (nm %in% names(beta_unscaled)) {
      sdj <- if (!is.null(recipe) && nm %in% names(recipe)) recipe[[nm]]$sd else 1
      out[[nm]] <- beta_unscaled[[nm]] * sdj
    }
  }
# Intercept: b0_scaled = b0_unscaled + sum_j b_j_unscaled * mean_j
  b0_un <- beta_unscaled[["(Intercept)"]] %||% 0
  adj <- 0
  for (nm in columns) {
    if (nm == "(Intercept)") next
    if (nm %in% names(beta_unscaled)) {
      mj <- if (!is.null(recipe) && nm %in% names(recipe)) recipe[[nm]]$mean else 0
      adj <- adj + beta_unscaled[[nm]] * mj
    }
  }
  out[["(Intercept)"]] <- b0_un + adj
  out
}

#' Map unscaled auxiliary multipliers to scaled space
#' @keywords internal
#' @param lambda_unscaled named numeric vector of auxiliary multipliers aligned
#'   to auxiliary design columns (no intercept) on original scale.
#' @param recipe `nmar_scaling_recipe`.
#' @param columns character vector of auxiliary column names (order) for the
#'   scaled design.
#' @return numeric vector of multipliers in the scaled space.
scale_aux_multipliers <- function(lambda_unscaled, recipe, columns) {
  if (is.null(lambda_unscaled) || length(lambda_unscaled) == 0) {
    return(setNames(numeric(length(columns)), columns))
  }
  lambda_unscaled <- lambda_unscaled[intersect(names(lambda_unscaled), columns)]
  out <- setNames(numeric(length(columns)), columns)
  for (nm in columns) {
    if (nm %in% names(lambda_unscaled)) {
      sdj <- if (!is.null(recipe) && nm %in% names(recipe)) recipe[[nm]]$sd else 1
      out[[nm]] <- lambda_unscaled[[nm]] / sdj
    }
  }
  out
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
  missing_in_recipe <- setdiff(param_names[param_names != intercept_name], names(recipe))
  if (length(missing_in_recipe)) {
    stop(
      "Scaling recipe is missing entries for coefficients: ",
      paste(missing_in_recipe, collapse = ", "),
      ". Ensure all predictors were present during scaling.",
      call. = FALSE
    )
  }
  unscaled_coeffs <- drop(D %*% scaled_coeffs)
  unscaled_vcov <- D %*% scaled_vcov %*% t(D)
  list(coefficients = unscaled_coeffs, vcov = unscaled_vcov)
}
