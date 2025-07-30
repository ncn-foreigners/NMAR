#' @importFrom bbmle mle2 coef
#' @importFrom stats as.formula dnorm dgamma sd setNames

generate_conditional_density <- function(model) {
  # Create data frame and extract covariate names
  data_df <- data.frame(y = model$y_1, model$x_for_y_obs)
  covar_names <- setdiff(colnames(data_df), "y")
  n_covars <- length(covar_names)

  if (model$y_dens == "normal") {
    # Build linear predictor string dynamically
    lin_pred_str <- "beta0"
    if (n_covars > 0) {
      beta_terms <- paste0("beta", 1:n_covars, " * ", covar_names, collapse = " + ")
      lin_pred_str <- paste(lin_pred_str, beta_terms, sep = " + ")
    }

    # Create start list with dynamic coefficients
    start_list <- list(beta0 = mean(model$y_1))
    if (n_covars > 0) {
      beta_starts <- setNames(rep(0.1, n_covars), paste0("beta", 1:n_covars))
      start_list <- c(start_list, as.list(beta_starts))
    }
    start_list$sigma <- sd(model$y_1)

    # Fit model with dynamically generated formula
    .model <- mle2(
      as.formula(paste("y ~ dnorm(mean =", lin_pred_str, ", sd = sigma)")),
      start = start_list,
      data = data_df
    )

    # Define density function
    density_fun <- function(y, x) {
      coefs <- bbmle::coef(.model)
      beta_names <- c("beta0", paste0("beta", 1:n_covars))
      beta_vec <- coefs[beta_names]

      # Create design matrix
      design_mat <- as.matrix(cbind(1, x[covar_names]))
      mean_val <- c(design_mat %*% beta_vec)

      dnorm(y, mean = mean_val, sd = coefs["sigma"])
    }

  } else if (model$y_dens == "gamma") {
    # Build linear predictor string (same as normal)
    lin_pred_str <- "beta0"
    if (n_covars > 0) {
      beta_terms <- paste0("beta", 1:n_covars, " * ", covar_names, collapse = " + ")
      lin_pred_str <- paste(lin_pred_str, beta_terms, sep = " + ")
    }

    # Create start list
    start_list <- list(beta0 = mean(model$y_1))
    if (n_covars > 0) {
      beta_starts <- setNames(rep(0.1, n_covars), paste0("beta", 1:n_covars))
      start_list <- c(start_list, as.list(beta_starts))
    }
    start_list$k <- 1

    # Create lower bounds
    lower_vec <- rep(-Inf, length(start_list))
    names(lower_vec) <- names(start_list)
    lower_vec[c("beta0", "k")] <- c(0.001, 0.001)

    # Fit model
    .model <- mle2(
      as.formula(paste("y ~ dgamma(shape = k, scale = (", lin_pred_str, ")/k)")),
      start = start_list,
      data = data_df,
      method = "L-BFGS-B",
      lower = lower_vec
    )

    # Define density function
    density_fun <- function(y, x) {
      coefs <- bbmle::coef(.model)
      beta_names <- c("beta0", paste0("beta", 1:n_covars))
      beta_vec <- coefs[beta_names]

      design_mat <- as.matrix(cbind(1, x[covar_names]))
      mean_val <- c(design_mat %*% beta_vec)
      shape <- coefs["k"]
      scale_val <- mean_val / shape

      dgamma(y, shape = shape, scale = scale_val)
    }

  } else {
    stop("Unsupported distribution type. Use 'normal' or 'gamma'.")
  }

  return(list(
    model = .model,
    density_function = density_fun
  ))
}

generate_conditional_density_matrix <- function(model) {
  matrix(
    outer(1:nrow(model$x_for_y_unobs),model$y_1 , FUN = function(i, y) model$density_fun(y, model$x_for_y_unobs[i, , drop = FALSE])),
    ncol = length(model$y_1),
    nrow = nrow(model$x_for_y_unobs)
  )
}

generate_C_matrix <- function(model) {
  stopifnot(
    is.matrix(model$f_matrix_nieobs)
  )

  col_sums <- colSums(model$f_matrix_nieobs)

  return(matrix(col_sums, ncol = 1))
}
