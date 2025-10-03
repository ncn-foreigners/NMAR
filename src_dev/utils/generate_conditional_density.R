generate_conditional_density <- function(model) {

  data_df <- data.frame(y = model$y_1, model$x_for_y_obs)
  data_df$weights <- model$respondent_weights # Add weights to data frame

# Get covariate names (excluding y and weights)
  covar_names <- setdiff(colnames(model$x_for_y_obs), c("y", "weights"))
  n_covars <- length(covar_names)

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
  chosen_dist <- model$y_dens
  if (chosen_dist == "auto") {
    aics <- c()
    fits <- list()
    for (d in names(dist_list)) {
# Create formula with explicit covariates (not using .)
      formula_str <- if (d == "lognormal") {
        paste("log(y) ~", paste(covar_names, collapse = " + "))
      } else {
        paste("y ~", paste(covar_names, collapse = " + "))
      }
      fit_try <- try(
        dist_list[[d]]$fit(as.formula(formula_str), data = data_df),
        silent = TRUE
      )
      if (!inherits(fit_try, "try-error")) {
        fits[[d]] <- fit_try
        aics[d] <- AIC(fit_try)
      }
    }
    chosen_dist <- names(which.min(aics))
    .model <- fits[[chosen_dist]]
  } else {
# Create formula with explicit covariates (not using .)
    formula_str <- if (chosen_dist == "lognormal") {
      paste("log(y) ~", paste(covar_names, collapse = " + "))
    } else {
      paste("y ~", paste(covar_names, collapse = " + "))
    }
    .model <- dist_list[[chosen_dist]]$fit(as.formula(formula_str), data_df)
  }

# Extract coefficients and create design matrix function
  coefs <- coef(.model)
  beta_names <- names(coefs)

# Design matrix function - only include covariates, not weights
  design_mat <- function(x) {
    as.matrix(cbind(Intercept = 1, x[covar_names]))
  }

# Add sigma coefficient for normal and lognormal
  if (chosen_dist %in% c("normal", "lognormal")) {
    sigma_val <- sd(resid(.model))
    coefs <- c(coefs, sigma = sigma_val)
  }

# Density function
  density_fun <- function(y, x) {
    x_mat <- design_mat(x)

# Ensure coefficient vector is aligned with design matrix
    if (ncol(x_mat) != length(coefs[beta_names])) {
      stop(paste("Design matrix has", ncol(x_mat), "columns but coefficients have",
                 length(coefs[beta_names]), "elements"))
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
# Add error handling for dimension mismatches
  tryCatch({
    matrix(
      outer(1:nrow(model$x_for_y_unobs), model$y_1,
            FUN = function(i, y) model$density_fun(y, model$x_for_y_unobs[i, , drop = FALSE])),
      ncol = length(model$y_1),
      nrow = nrow(model$x_for_y_unobs)
    )
  }, error = function(e) {
    stop("Error in generate_conditional_density_matrix: ", e$message)
  })
}

generate_C_matrix <- function(model) {
  stopifnot(is.matrix(model$f_matrix_nieobs))
  col_sums <- colSums(model$f_matrix_nieobs)
  return(matrix(col_sums, ncol = 1))
}
