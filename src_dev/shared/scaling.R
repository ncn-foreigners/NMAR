#' Shared scaling infrastructure for NMAR engines
#'
#' Centralized feature scaling and parameter unscaling routines used by NMAR
#' estimation engines to ensure consistent, numerically stable behavior.
#'
#' @section Goals:
#' \itemize{
#'   \item Provide an engine-agnostic API for standardizing design matrices and
#'     auxiliary moments before solving.
#'   \item Return a minimal "recipe" (per-column mean and standard deviation)
#'     for unscaling coefficients and covariance matrices after solving.
#' }
#'
#' @section Inputs/Outputs:
#' \describe{
#'   \item{Inputs}{\code{Z_un} (response model matrix with intercept),
#'     optional \code{X_un} (auxiliary model matrix, no intercept),
#'     optional named \code{mu_x_un} (auxiliary means on the original scale),
#'     and a logical \code{standardize} flag.}
#'   \item{Outputs}{Scaled matrices \code{Z}, \code{X}, and \code{mu_x}, plus an
#'     \code{nmar_scaling_recipe} used later for unscaling.}
#' }
#'
#' @section Integration pattern:
#' \enumerate{
#'   \item Before solving: call \code{validate_and_apply_nmar_scaling()} (engine-level)
#'         or \code{prepare_nmar_scaling()} (low-level) to obtain scaled matrices
#'         and recipe.
#'   \item Solve in the scaled space.
#'   \item After solving: call \code{unscale_coefficients()} to unscale coefficients
#'         and their covariance matrices.
#'   \item Store the \code{nmar_scaling_recipe} in results for diagnostics and reproducibility.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item The intercept column is never scaled.
#'   \item Columns with near-zero variance are centered but assigned
#'     \code{sd = 1} so that the corresponding parameter is not inflated by
#'     division by a very small standard deviation.
#'   \item Engines may use design-weighted scaling via the \code{weights} and
#'         \code{weight_mask} arguments.
#' }
#'
#' @name nmar_scaling_infra
#' @keywords internal
NULL

#' @keywords internal
new_nmar_scaling_recipe <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = "nmar_scaling_recipe")
}

#' @keywords internal
validate_nmar_scaling_recipe <- function(x) {
  if (!all(vapply(x, is.list, logical(1)))) stop("Each element of a nmar_scaling_recipe must be a list.", call. = FALSE)
  has_req_fields <- vapply(x, function(comp) all(c("mean", "sd") %in% names(comp)), logical(1))
  if (!all(has_req_fields)) stop("Each element of a nmar_scaling_recipe must contain 'mean' and 'sd'.", call. = FALSE)
  x
}

#' Compute (possibly weighted) mean and standard deviation
#'
#' @keywords internal
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
  mask <- is.finite(values) & is.finite(weights)
  if (!any(mask)) {
    return(list(mean = NA_real_, sd = NA_real_))
  }
  w <- weights[mask]
  x <- values[mask]
  if (any(w < 0)) stop("`weights` must be nonnegative.", call. = FALSE)
  w_sum <- sum(w)
  if (w_sum <= 0 || all(w == 0)) {
    return(list(mean = NA_real_, sd = NA_real_))
  }
  mean_val <- sum(w * x) / w_sum
# Population-style weighted standard deviation: sqrt(E[(X - mean)^2]).
  sd_val <- sqrt(sum(w * (x - mean_val)^2) / w_sum)
  list(mean = mean_val, sd = sd_val)
}

#' Build a scaling recipe from one or more design matrices
#' @param ... one or more matrices with named columns.
#' @param intercept_col Intercept column name that should remain unscaled.
#' @param weights Optional numeric vector of weights used to compute weighted
#'   means/standard deviations.
#' @param weight_mask Optional logical mask or nonnegative numeric multipliers
#'   applied to `weights` before computing moments (useful for respondents-only
#'   scaling). If `weights` is `NULL`, `weight_mask` is treated as the weights.
#' @param tol_constant Numeric tolerance below which columns are treated as
#'   constant and left unscaled.
#' @param warn_on_constant Logical; emit a warning when a column is treated as
#'   constant.
#'
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
    if (is.logical(weight_mask)) {
      if (any(is.na(weight_mask))) stop("`weight_mask` must not contain NA.", call. = FALSE)
      mask_w <- as.numeric(weight_mask)
    } else {
      mask_w <- as.numeric(weight_mask)
      if (any(!is.finite(mask_w))) stop("`weight_mask` must be finite (no NA/NaN/Inf).", call. = FALSE)
      if (any(mask_w < 0)) stop("`weight_mask` must be nonnegative.", call. = FALSE)
    }
    if (!is.null(weights)) {
      if (length(mask_w) != length(weights)) {
        stop("`weight_mask` must have the same length as `weights`.", call. = FALSE)
      }
      weights <- weights * mask_w
    } else {
      weights <- mask_w
    }
  }
  if (!is.null(weights)) {
    w_finite <- weights[is.finite(weights)]
    if (!any(w_finite > 0)) {
      stop("Scaling weights contain no positive entries after applying weight_mask.", call. = FALSE)
    }
  }
  for (mat in matrices) {
    for (col_name in colnames(mat)) {
      if (col_name == intercept_col || col_name %in% names(recipe)) next
      col_data <- mat[, col_name]
      stats <- compute_weighted_stats(col_data, weights)
      col_mean <- stats$mean
      col_sd <- stats$sd
      if (!is.finite(col_sd) || !is.finite(col_mean)) {
        if (warn_on_constant) {
          warning(sprintf("Column '%s' has undefined weighted moments; leaving as identity scale.", col_name), call. = FALSE)
        }
# Use an identity transform to avoid propagating NA values.
        recipe[[col_name]] <- list(mean = 0, sd = 1)
      } else if (col_sd < tol_constant) {
        if (warn_on_constant) {
          warning(sprintf("Column '%s' is nearly constant under the scaling weights; centering to zero column.", col_name), call. = FALSE)
        }
# Center to a zero column; coefficients on such columns are not
# identifiable and should not be inflated by dividing by a very small sd.
        recipe[[col_name]] <- list(mean = col_mean, sd = 1)
      } else {
        recipe[[col_name]] <- list(mean = col_mean, sd = col_sd)
      }
    }
  }
  validate_nmar_scaling_recipe(new_nmar_scaling_recipe(recipe))
}

#' Apply scaling to a matrix using a recipe
#' @param matrix_to_scale a numeric matrix with column names present in `recipe`.
#' @param recipe an object of class `nmar_scaling_recipe`.
#' @return a matrix with each named column centered and scaled using the recipe
#'
#' @keywords internal
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
#' @param Z_un response model matrix (with intercept column).
#' @param X_un auxiliary model matrix (no intercept), or NULL.
#' @param mu_x_un named numeric vector of auxiliary means on the original scale
#'   (names must match `colnames(X_un)`), or NULL.
#' @param standardize logical; apply standardization if TRUE.
#' @param weights Optional numeric vector used for weighted scaling.
#' @param weight_mask Optional logical mask or nonnegative numeric multipliers
#'   applied to `weights`.
#' @return a list with `Z`, `X`, `mu_x`, and `recipe`.
#'
#' @keywords internal
prepare_nmar_scaling <- function(Z_un, X_un, mu_x_un, standardize,
                                 weights = NULL, weight_mask = NULL) {
  if (!standardize) {
    return(list(Z = Z_un, X = X_un, mu_x = if (is.null(mu_x_un)) numeric(0) else mu_x_un, recipe = NULL))
  }
# Check row consistency between response and auxiliary matrices.
  if (!is.null(X_un) && nrow(Z_un) != nrow(X_un)) {
    stop("`Z_un` and `X_un` must have the same number of rows.", call. = FALSE)
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
    mu_x <- vapply(names(mu_x), function(n) {
      if (n %in% names(recipe)) {
        (mu_x[[n]] - recipe[[n]]$mean) / recipe[[n]]$sd
      } else {
        warning("No scaling info for mu_x component '", n, "'. Leaving unscaled.")
        mu_x[[n]]
      }
    }, numeric(1))
  }
  list(Z = Z, X = X, mu_x = mu_x, recipe = recipe)
}

#' Validate and apply scaling (engine-friendly)
#' @param standardize logical; apply standardization if TRUE.
#' @param has_aux logical; whether the engine uses auxiliary constraints.
#' @param response_model_matrix_unscaled response model matrix (with intercept).
#' @param aux_matrix_unscaled auxiliary matrix (no intercept) or an empty matrix.
#' @param mu_x_unscaled named auxiliary means on original scale, or NULL.
#' @param weights Optional numeric vector used for weighted scaling.
#' @param weight_mask Optional logical mask or nonnegative numeric multipliers
#'   applied to `weights`.
#' @return a list with `nmar_scaling_recipe`, `response_model_matrix_scaled`,
#'   `auxiliary_matrix_scaled`, and `mu_x_scaled`.
#'
#' @keywords internal
validate_and_apply_nmar_scaling <- function(standardize, has_aux, response_model_matrix_unscaled,
                                            aux_matrix_unscaled, mu_x_unscaled,
                                            weights = NULL, weight_mask = NULL) {
  nmar_scaling_recipe <- NULL
  if (standardize) {
# Coerce NULL auxiliaries to a 0-column matrix aligned to the response
# matrix to simplify downstream logic.
    if (is.null(aux_matrix_unscaled)) {
      aux_matrix_unscaled <- matrix(nrow = nrow(response_model_matrix_unscaled), ncol = 0,
        dimnames = list(NULL, character()))
    }
# Check row consistency between response and auxiliary matrices.
    if (nrow(response_model_matrix_unscaled) != nrow(aux_matrix_unscaled)) {
      stop("Response and auxiliary matrices must have the same number of rows.", call. = FALSE)
    }
    if (has_aux) {
# Drop any stray intercept column from the auxiliary matrix before
# matching names to auxiliary means.
      if ("(Intercept)" %in% colnames(aux_matrix_unscaled)) {
        aux_matrix_unscaled <- aux_matrix_unscaled[, setdiff(colnames(aux_matrix_unscaled), "(Intercept)"), drop = FALSE]
      }
      if (!setequal(colnames(aux_matrix_unscaled), names(mu_x_unscaled))) {
        stop("Names of `auxiliary_means` do not match the variables specified on the RHS of the formula.")
      }
      mu_x_unscaled <- mu_x_unscaled[colnames(aux_matrix_unscaled)]
    }
    if (!is.null(weights) && length(weights) != nrow(response_model_matrix_unscaled)) {
      stop("`weights` must have the same length as the number of rows in the response matrix.", call. = FALSE)
    }
    if (!is.null(weight_mask) && length(weight_mask) != nrow(response_model_matrix_unscaled)) {
      stop("`weight_mask` must have the same length as the number of rows in the response matrix.", call. = FALSE)
    }
    nmar_scaling_recipe <- create_nmar_scaling_recipe(
      response_model_matrix_unscaled,
      aux_matrix_unscaled,
      weights = weights,
      weight_mask = weight_mask
    )
    response_model_matrix_scaled <- apply_nmar_scaling(response_model_matrix_unscaled, nmar_scaling_recipe)
    auxiliary_matrix_scaled <- apply_nmar_scaling(aux_matrix_unscaled, nmar_scaling_recipe)
    mu_x_scaled <- if (has_aux && !is.null(mu_x_unscaled) && length(mu_x_unscaled) > 0) {
      vapply(names(mu_x_unscaled), function(n) (mu_x_unscaled[n] - nmar_scaling_recipe[[n]]$mean) / nmar_scaling_recipe[[n]]$sd, numeric(1))
    } else {
      numeric(0)
    }
  } else {
    response_model_matrix_scaled <- response_model_matrix_unscaled
    auxiliary_matrix_scaled <- aux_matrix_unscaled
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
#' @param beta_unscaled named numeric vector of coefficients for the response
#'   model on the original scale, including an intercept named `"(Intercept)"`.
#' @param recipe `nmar_scaling_recipe` returned by scaling utilities.
#' @param columns character vector of column names (order) for the scaled design
#'   matrix (including intercept).
#' @return numeric vector of coefficients in the scaled space, ordered by
#'   `columns`.
#'
#' @keywords internal
scale_coefficients <- function(beta_unscaled, recipe, columns) {
  if (is.null(beta_unscaled) || length(beta_unscaled) == 0) {
    return(setNames(numeric(length(columns)), columns))
  }
  beta_unscaled <- beta_unscaled[intersect(names(beta_unscaled), columns)]
  out <- setNames(numeric(length(columns)), columns)
# Non-intercept: multiply by sd.
  for (nm in columns) {
    if (nm == "(Intercept)") next
    if (nm %in% names(beta_unscaled)) {
      sdj <- if (!is.null(recipe) && nm %in% names(recipe)) recipe[[nm]]$sd else 1
      out[[nm]] <- beta_unscaled[[nm]] * sdj
    }
  }
# Intercept: b0_scaled = b0_unscaled + sum_j b_j_unscaled * mean_j.
  b0_un <- if (!is.null(beta_unscaled[["(Intercept)"]])) beta_unscaled[["(Intercept)"]] else 0
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
#' @param lambda_unscaled named numeric vector of auxiliary multipliers aligned
#'   to auxiliary design columns (no intercept) on original scale.
#' @param recipe `nmar_scaling_recipe`.
#' @param columns character vector of auxiliary column names (order) for the
#'   scaled design.
#' @return numeric vector of multipliers in the scaled space.
#'
#' @keywords internal
scale_aux_multipliers <- function(lambda_unscaled, recipe, columns) {
  if (is.null(lambda_unscaled) || length(lambda_unscaled) == 0) {
    return(setNames(numeric(length(columns)), columns))
  }
  lambda_unscaled <- lambda_unscaled[intersect(names(lambda_unscaled), columns)]
  out <- setNames(numeric(length(columns)), columns)
  for (nm in columns) {
    if (nm %in% names(lambda_unscaled)) {
# Scaling identity for centered auxiliaries: if
#   Xc_scaled = (Xc_unscaled / sd)
# and we want to preserve the linear term Xc_unscaled %*% lambda_un
# under scaling, then
#   Xc_un %*% lambda_un = Xc_scaled %*% lambda_scaled
# implies lambda_scaled = sd * lambda_un.
      sdj <- if (!is.null(recipe) && nm %in% names(recipe)) recipe[[nm]]$sd else 1
      out[[nm]] <- lambda_unscaled[[nm]] * sdj
    }
  }
  out
}

#' Unscale regression coefficients and covariance
#' @param scaled_coeffs named numeric vector of coefficients estimated on the scaled space.
#' @param scaled_vcov covariance matrix of `scaled_coeffs`.
#' @param recipe `nmar_scaling_recipe` produced when scaling was applied.
#' @return a list with unscaled `coefficients` and `vcov`.
#'
#' @keywords internal
unscale_coefficients <- function(scaled_coeffs, scaled_vcov, recipe) {
  if (is.null(recipe)) {
    return(list(coefficients = scaled_coeffs, vcov = scaled_vcov))
  }
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
    warning(
      "Scaling recipe is missing entries for coefficients: ",
      paste(missing_in_recipe, collapse = ", "),
      ". Treating as unscaled (identity).",
      call. = FALSE
    )
  }
  unscaled_coeffs <- drop(D %*% scaled_coeffs)
  unscaled_vcov <- D %*% scaled_vcov %*% t(D)
# Preserve attribute semantics: if the input covariance had no dimnames,
# do not add them on output.
  if (is.null(dimnames(scaled_vcov))) dimnames(unscaled_vcov) <- NULL
  list(coefficients = unscaled_coeffs, vcov = unscaled_vcov)
}
