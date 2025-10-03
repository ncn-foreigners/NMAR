generate_conditional_density <- function(model) {
  # Conditional density fit for f_1(y | x_1; gamma)
  #
  # Rationale and safeguards (what can go wrong):
  # - Support mismatch (lognormal/exponential): both families require y>0.
  #   If any observed respondent outcome is nonpositive, fitting
  #   log(y) ~ x_1 (lognormal) or Gamma(log) (exponential) will produce domain
  #   errors or NaNs. We therefore exclude these families from the "auto"
  #   selection unless all observed y are strictly positive. If the user
  #   explicitly requests one of these but the data violate support, we warn and
  #   revert to the normal family.
  # - Non-finite likelihood evaluations: even with a feasible family, extreme
  #   parameter estimates (e.g., sigma <= 0) or numerical issues can yield
  #   non-finite density values during EM. After selecting a family we probe the
  #   observed support; if any evaluation is non-finite we conservatively fall
  #   back to the normal family at the same design to keep the EM step defined.
  # - Design-matrix alignment: the design for f_1(y | x_1) must match the
  #   fitted coefficient vector. We check dimensions and throw an explicit error
  #   if they disagree (prevents silent recycling or misaligned multiplications).
  # - Scaling: f_1 is fit on the scaled auxiliary space when the ET solver works
  #   in scaled coordinates. For density evaluation on unobserved x_1 we apply
  #   the stored scaling recipe so the design fed to the density matches the fit.

  data_df <- data.frame(y = model$y_1, model$x_for_y_obs)
  # Use respondent slice of design weights if available; otherwise weight 1
  if (!is.null(model$design_weights) && !is.null(model$respondent_mask)) {
    w_resp_local <- model$design_weights[model$respondent_mask]
  } else {
    w_resp_local <- rep(1, length(model$y_1))
  }
  data_df$weights <- w_resp_local

  # Get covariate names (excluding y and weights)
  covar_names <- setdiff(colnames(model$x_for_y_obs), c("y", "weights"))
  n_covars <- length(covar_names)
  rhs_terms <- if (n_covars > 0) paste(covar_names, collapse = " + ") else "1"

  respondent_support_positive <- all(is.finite(model$y_1) & model$y_1 > 0)

  # Gradient and Hessian for normal distribution
  normal_gradient <- function(y, mean_val, coefs, x_vector) {
    sigma <- coefs[["sigma"]]
    residual <- y - mean_val
    d_beta <- (residual / sigma^2) * x_vector
    d_sigma <- -1 / sigma + residual^2 / sigma^3
    gradient_vector <- c(d_beta, sigma = d_sigma)
    names(gradient_vector) <- c(paste0("beta", 0:(length(x_vector) - 1)), "sigma")
    return(gradient_vector)
  }

  normal_hessian <- function(y, mean_val, coefs, x_vector) {
    sigma <- coefs[["sigma"]]
    residual <- y - mean_val
    n_params <- length(x_vector) + 1
    H <- matrix(0, nrow = n_params, ncol = n_params)
    param_names <- c(paste0("beta", 0:(length(x_vector) - 1)), "sigma")
    rownames(H) <- colnames(H) <- param_names
    for (i in 1:length(x_vector)) {
      for (j in 1:length(x_vector)) {
        H[i, j] <- -x_vector[i] * x_vector[j] / sigma^2
      }
    }
    for (i in 1:length(x_vector)) {
      H[i, n_params] <- H[n_params, i] <- -2 * residual * x_vector[i] / sigma^3
    }
    H[n_params, n_params] <- 1 / sigma^2 - 3 * residual^2 / sigma^4
    return(H)
  }

  # Supported distributions with families and link functions
  dist_list <- list(
    normal = list(
      family = gaussian(link = "identity"),
      fit = function(formula, data) {
        fit <- lm(formula, data = data, weights = weights)
        return(fit)
      },
      extra = "sigma",
      density = function(y, mean_val, coefs) {
        dnorm(y, mean = mean_val, sd = coefs[["sigma"]])
      },
      gradient = normal_gradient,
      hessian = normal_hessian
    ),
    lognormal = list(
      family = gaussian(link = "identity"),
      transform = log,
      fit = function(formula, data) {
        fit <- lm(formula, data = data, weights = weights)
        return(fit)
      },
      extra = "sigma",
      density = function(y, mean_val, coefs) {
        dlnorm(y, meanlog = mean_val, sdlog = coefs[["sigma"]])
      }
    ),
    exponential = list(
      family = Gamma(link = "log"),
      fit = function(formula, data) {
        fit <- glm(formula, data = data, family = Gamma(link = "log"), weights = weights)
        return(fit)
      },
      density = function(y, mean_val, coefs) {
        dexp(y, rate = 1 / mean_val)
      }
    )
  )

  # Choose distribution
  fit_distribution <- function(dist_name) {
    formula_str <- if (dist_name == "lognormal") {
      paste("log(y) ~", rhs_terms)
    } else {
      paste("y ~", rhs_terms)
    }
    dist_list[[dist_name]]$fit(as.formula(formula_str), data_df)
  }

  chosen_dist <- model$y_dens
  if (chosen_dist == "auto") {
    candidate_dists <- names(dist_list)
    if (!respondent_support_positive) {
      candidate_dists <- setdiff(candidate_dists, c("lognormal", "exponential"))
    }
    aics <- c()
    fits <- list()
    for (d in candidate_dists) {
      fit_try <- try(fit_distribution(d), silent = TRUE)
      if (!inherits(fit_try, "try-error")) {
        fits[[d]] <- fit_try
        aics[d] <- AIC(fit_try)
      }
    }
    if (!length(fits)) {
      stop("Unable to fit any conditional outcome density.", call. = FALSE)
    }
    chosen_dist <- names(which.min(aics))
    .model <- fits[[chosen_dist]]
  } else {
    if (chosen_dist %in% c("lognormal", "exponential") && !respondent_support_positive) {
      warning("Selected y_dens='", chosen_dist, "' requires positive observed outcomes; using normal instead.", call. = FALSE)
      chosen_dist <- "normal"
    }
    .model <- fit_distribution(chosen_dist)
  }

  # Extract coefficients and create design matrix function
  coefs <- coef(.model)
  beta_names <- names(coefs)

  # Design matrix function - only include covariates, not weights
  design_mat <- function(x) {
    as.matrix(cbind(Intercept = 1, x[covar_names]))
  }

  # Add sigma coefficient for normal and lognormal using the model's
  # weighted residual standard error (design-consistent when weights are
  # provided to lm()). Fallback to an epsilon if undefined
  if (chosen_dist %in% c("normal", "lognormal")) {
    sigma_val <- tryCatch(
      {
        s <- summary(.model)$sigma
        if (!is.finite(s) || s <= 0) NA_real_ else s
      },
      error = function(e) NA_real_
    )
    if (!is.finite(sigma_val) || is.na(sigma_val) || sigma_val <= 0) {
      sigma_val <- .Machine$double.eps
    }
    coefs <- c(coefs, sigma = sigma_val)
  }

  # Density function
  density_fun <- function(y, x) {
    x_mat <- design_mat(x)

    # Ensure coefficient vector is aligned with design matrix
    if (ncol(x_mat) != length(coefs[beta_names])) {
      stop(paste(
        "Design matrix has", ncol(x_mat), "columns but coefficients have",
        length(coefs[beta_names]), "elements"
      ))
    }

    mean_val <- x_mat %*% coefs[beta_names]
    dist_list[[chosen_dist]]$density(y, mean_val, coefs)
  }

  # Gradient and Hessian functions (only for normal)
  density_grad_fun <- if (chosen_dist == "normal") {
    function(y, x) {
      x_vector <- c(1, as.numeric(x[covar_names]))
      mean_val <- sum(x_vector * coefs[beta_names])
      normal_gradient(y, mean_val, coefs, x_vector)
    }
  } else {
    NULL
  }

  density_hess_fun <- if (chosen_dist == "normal") {
    function(y, x) {
      x_vector <- c(1, as.numeric(x[covar_names]))
      mean_val <- sum(x_vector * coefs[beta_names])
      normal_hessian(y, mean_val, coefs, x_vector)
    }
  } else {
    NULL
  }

  # Guard against densities that fail on the observed support. If detected and
  # a non-normal family was selected, refit using the normal helper instead
  if (!identical(chosen_dist, "normal")) {
    finite_eval <- tryCatch(
      {
        vals <- vapply(seq_along(model$y_1), function(idx) {
          x_row <- model$x_for_y_obs[idx, , drop = FALSE]
          if (isFALSE(model$features_are_scaled) && !is.null(model$nmar_scaling_recipe)) {
            x_row <- apply_nmar_scaling(x_row, model$nmar_scaling_recipe)
          }
          density_fun(model$y_1[idx], x_row)
        }, numeric(1))
        all(is.finite(vals))
      },
      error = function(e) FALSE
    )

    if (!finite_eval) {
      warning("Density '", chosen_dist, "' produced non-finite evaluations; reverting to normal.", call. = FALSE)
      model_fallback <- model
      model_fallback$y_dens <- "normal"
      return(generate_conditional_density(model_fallback))
    }
  }

  return(list(
    model = .model,
    density_function = density_fun,
    density_function_grad = density_grad_fun,
    density_function_hess = density_hess_fun,
    chosen_distribution = chosen_dist,
    num_of_coefs = length(coefs)
  ))
}
generate_conditional_density_matrix <- function(model) {
  # Add error handling for dimension mismatches and ensure that the covariates
  # fed to the fitted density function live on the SAME scale used during the
  # density fit. The model stores a scaling recipe, after unscaling the feature
  # matrices for reporting, we temporarily re-apply that recipe here
  tryCatch({
    n_unobs <- nrow(model$x_for_y_unobs)
    n_resp <- length(model$y_1)
    if (!n_unobs || !n_resp) {
      return(matrix(numeric(0), nrow = n_unobs, ncol = n_resp))
    }

    y_vals <- as.numeric(model$y_1)
    out <- matrix(NA_real_, nrow = n_unobs, ncol = n_resp)

    for (i in seq_len(n_unobs)) {
      x_row <- model$x_for_y_unobs[i, , drop = FALSE]
      # Re-apply scaling if the current feature matrices are on the original scale
      if (isFALSE(model$features_are_scaled) && !is.null(model$nmar_scaling_recipe)) {
        x_row <- apply_nmar_scaling(x_row, model$nmar_scaling_recipe)
      }
      # Evaluate densities for all y_j in a single vectorized call
      out[i, ] <- vapply(y_vals, function(yj) model$density_fun(yj, x_row), numeric(1))
    }
    out
  }, error = function(e) {
    stop("Error in generate_conditional_density_matrix: ", e$message)
  })
}

generate_C_matrix <- function(model) {
  # C(y_j; gamma) must sum over RESPONDENTS, not nonrespondents.
  # Paper (eq. 15): C(y; gamma) = sum_{l: delta_l=1} f1(y | x_l; gamma). In the survey case,
  # this becomes a weighted sum with design weights d_l. This term depends only
  # on y_j and gamma_hat; it must not depend on any nonrespondent covariates.
  #
  # Earlier implementation incorrectly computed C as colSums of f_matrix_nieobs,
  # which sums over nonrespondents i: sum_i f1(y_j | x_i^unobs). That is
  # inconsistent with the derivation and breaks the fractional weights.

  y_obs <- as.numeric(model$y_1)
  X_obs <- model$x_for_y_obs
  n_obs <- length(y_obs)
  if (!length(y_obs) || is.null(X_obs) || !nrow(X_obs)) {
    return(matrix(numeric(0), ncol = 1))
  }

  # Build F_resp[l, j] = f1(y_j | x_l^obs; gamma_hat), honoring the scale used during fit.
  F_resp <- matrix(NA_real_, nrow = nrow(X_obs), ncol = n_obs)
  for (l in seq_len(nrow(X_obs))) {
    x_row <- X_obs[l, , drop = FALSE]
    if (isFALSE(model$features_are_scaled) && !is.null(model$nmar_scaling_recipe)) {
      x_row <- apply_nmar_scaling(x_row, model$nmar_scaling_recipe)
    }
    F_resp[l, ] <- vapply(y_obs, function(yj) model$density_fun(yj, x_row), numeric(1))
  }

  # Include respondent sampling weights if present (Remark 2)
  if (!is.null(model$design_weights) && !is.null(model$respondent_mask)) {
    w_resp <- model$design_weights[model$respondent_mask]
  } else {
    w_resp <- rep(1, nrow(X_obs))
  }

  # C_j = sum_l w_l * f1(y_j | x_l)
  C_vec <- as.numeric(crossprod(w_resp, F_resp))
  matrix(C_vec, ncol = 1)
}
