generate_conditional_density <- function(model) {
  data_df <- data.frame(y = model$y_1, model$x_for_y_obs)

  if (!is.null(model$respondent_weights)) {
    data_df$weights <- model$respondent_weights
  }

  covar_names <- colnames(model$x_for_y_obs)
  n_covars <- length(covar_names)

  dist_list <- list(
    normal = list(
      family = gaussian(link = "identity"),
      fit = function(formula, data, weights = NULL) {
        if (!is.null(weights)) {
          fit <- glm(formula, data = data, family = gaussian(), weights = weights)
        } else {
          fit <- glm(formula, data = data, family = gaussian())
        }
        return(fit)
      },
      extra = "sigma",
      density = function(y, mean_val, coefs) {
        dnorm(y, mean = mean_val, sd = coefs[["sigma"]])
      }
    ),
    lognormal = list(
      family = gaussian(link = "identity"),
      transform = log,
      fit = function(formula, data, weights = NULL) {
        if (!is.null(weights)) {
          fit <- lm(formula, data = data, weights = weights)
        } else {
          fit <- lm(formula, data = data)
        }
        return(fit)
      },
      extra = "sigma",
      density = function(y, mean_val, coefs) {
        dlnorm(y, meanlog = mean_val, sdlog = coefs[["sigma"]])
      }
    ),
    exponential = list(
      family = Gamma(link = "log"),
      fit = function(formula, data, weights = NULL) {
        if (!is.null(weights)) {
          fit <- glm(formula, data = data, family = Gamma(link = "log"), weights = weights)
        } else {
          fit <- glm(formula, data = data, family = Gamma(link = "log"))
        }
        return(fit)
      },
      density = function(y, mean_val, coefs) {
        dexp(y, rate = 1 / mean_val)
      }
    )
  )

  chosen_dist <- model$y_dens

  formula_str <- if (chosen_dist == "lognormal") {
    paste("log(y) ~", paste(covar_names, collapse = " + "))
  } else {
    paste("y ~", paste(covar_names, collapse = " + "))
  }

  if (!is.null(model$respondent_weights)) {
    .model <- dist_list[[chosen_dist]]$fit(as.formula(formula_str), data_df, weights = data_df$weights)
  } else {
    .model <- dist_list[[chosen_dist]]$fit(as.formula(formula_str), data_df)
  }

  coefs <- coef(.model)
  beta_names <- names(coefs)

  design_mat <- function(x) {
    as.matrix(cbind(Intercept = 1, x[covar_names]))
  }

  if (chosen_dist %in% c("normal", "lognormal")) {
    sigma_val <- sd(resid(.model))
    coefs <- c(coefs, sigma = sigma_val)
  }

  density_fun <- function(y, x) {
    x_mat <- design_mat(x)

    validator$assert_matrix_ncol(x_mat, length(coefs[beta_names]), name = "design_mat/coefs")

    mean_val <- x_mat %*% coefs[beta_names]
    dist_list[[chosen_dist]]$density(y, mean_val, coefs)
  }

  return(structure(list(
    density_model = .model,
    density_function = density_fun,
    chosen_distribution = chosen_dist,
    num_of_coefs = length(coefs),
    aic_comparison = if (exists("aics")) aics else NULL
  ), class = "nmar_density_response"))
}

generate_conditional_density_matrix <- function(model) {
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
  f_matrix_obs <- t(sapply(1:nrow(model$x_for_y_obs), function(i) {
    model$density_fun(model$y_1, model$x_for_y_obs[i, , drop = FALSE])
  }))

  C_vector <- colSums(f_matrix_obs)

  return(matrix(C_vector, ncol = 1))
}
