#' @importFrom bbmle mle2 coef
#' @importFrom stats as.formula dnorm dgamma dlnorm dweibull dexp dt dcauchy sd setNames

generate_conditional_density <- function(model) {
  data_df <- data.frame(y = model$y_1, model$x_for_y_obs)
  covar_names <- setdiff(colnames(data_df), "y")
  n_covars <- length(covar_names)

  # supported distributions
  dist_list <- list(
    normal = list(
      formula = function(lin_pred_str) {
        as.formula(paste("y ~ dnorm(mean =", lin_pred_str, ", sd = sigma)"))
      },
      start = function(y) list(beta0 = mean(y), sigma = sd(y)),
      density = function(y, mean_val, coefs) {
        dnorm(y, mean = mean_val, sd = coefs[["sigma"]])
      },
      extra = c("sigma")
    ),
    lognormal = list(
      formula = function(lin_pred_str) {
        as.formula(paste("y ~ dlnorm(meanlog =", lin_pred_str, ", sdlog = sigma)"))
      },
      start = function(y) list(beta0 = log(mean(y)), sigma = sd(log(y))),
      density = function(y, mean_val, coefs) {
        dlnorm(y, meanlog = mean_val, sdlog = coefs[["sigma"]])
      },
      extra = c("sigma")
    ),
    weibull = list(
      formula = function(lin_pred_str) {
        as.formula(paste("y ~ dweibull(shape = shape, scale =", lin_pred_str, ")"))
      },
      start = function(y) list(beta0 = mean(y), shape = 1),
      density = function(y, mean_val, coefs) {
        dweibull(y, shape = coefs[["shape"]], scale = mean_val)
      },
      extra = c("shape")
    ),
    exponential = list(
      formula = function(lin_pred_str) {
        as.formula(paste("y ~ dexp(rate = 1/(", lin_pred_str, "))"))
      },
      start = function(y) list(beta0 = mean(y)),
      density = function(y, mean_val, coefs) {
        dexp(y, rate = 1/mean_val)
      },
      extra = character(0)
    ),
    cauchy = list(
      formula = function(lin_pred_str) {
        as.formula(paste("y ~ dcauchy(location =", lin_pred_str, ", scale = gamma)"))
      },
      start = function(y) list(beta0 = median(y), gamma = 1),
      density = function(y, mean_val, coefs) {
        dcauchy(y, location = mean_val, scale = coefs[["gamma"]])
      },
      extra = c("gamma")
    )
  )


  # choose distribution
  chosen_dist <- model$y_dens
  if (chosen_dist == "auto") {
    cat('AUTO SELECTED')
    aics <- c()
    fits <- list()
    for (d in names(dist_list)) {
      lin_pred_str <- "beta0"
      if (n_covars > 0) {
        beta_terms <- paste0("beta", 1:n_covars, " * ", covar_names, collapse = " + ")
        lin_pred_str <- paste(lin_pred_str, beta_terms, sep = " + ")
      }
      start_list <- dist_list[[d]]$start(model$y_1)
      if (n_covars > 0) {
        beta_starts <- setNames(rep(0.1, n_covars), paste0("beta", 1:n_covars))
        start_list <- c(start_list, as.list(beta_starts))
      }
      lower_vec <- NULL
      if (!is.null(dist_list[[d]]$lower)) {
        lower_vec <- dist_list[[d]]$lower(start_list)
      }
      fit_try <- try(
        mle2(
          dist_list[[d]]$formula(lin_pred_str),
          start = start_list,
          data = data_df,
          method = if (is.null(lower_vec)) "BFGS" else "L-BFGS-B",
          lower = lower_vec
        ),
        silent = TRUE
      )
      if (!inherits(fit_try, "try-error")) {
        fits[[d]] <- fit_try
        aics[d] <- AIC(fit_try)
      }
    }
    print(aics)
    print(fits)
    #print aics

    chosen_dist <- names(which.min(aics))
    .model <- fits[[chosen_dist]]
  } else {
    lin_pred_str <- "beta0"
    if (n_covars > 0) {
      beta_terms <- paste0("beta", 1:n_covars, " * ", covar_names, collapse = " + ")
      lin_pred_str <- paste(lin_pred_str, beta_terms, sep = " + ")
    }
    start_list <- dist_list[[chosen_dist]]$start(model$y_1)
    if (n_covars > 0) {
      beta_starts <- setNames(rep(0.1, n_covars), paste0("beta", 1:n_covars))
      start_list <- c(start_list, as.list(beta_starts))
    }
    lower_vec <- NULL
    if (!is.null(dist_list[[chosen_dist]]$lower)) {
      lower_vec <- dist_list[[chosen_dist]]$lower(start_list)
    }
    .model <- mle2(
      dist_list[[chosen_dist]]$formula(lin_pred_str),
      start = start_list,
      data = data_df,
      method = if (is.null(lower_vec)) "BFGS" else "L-BFGS-B",
      lower = lower_vec
    )
  }

  # density function
  density_fun <- function(y, x) {
    coefs <- bbmle::coef(.model)
    beta_names <- c("beta0", paste0("beta", 1:n_covars))
    beta_vec <- coefs[beta_names]
    design_mat <- as.matrix(cbind(1, x[covar_names]))
    mean_val <- c(design_mat %*% beta_vec)
    dist_list[[chosen_dist]]$density(y, mean_val, coefs)
  }


  return(list(
    model = .model,
    density_function = density_fun,
    chosen_distribution = chosen_dist
  ))
}


generate_conditional_density_matrix <- function(model) {
  matrix(
    outer(1:nrow(model$x_for_y_unobs), model$y_1,
          FUN = function(i, y) model$density_fun(y, model$x_for_y_unobs[i, , drop = FALSE])),
    ncol = length(model$y_1),
    nrow = nrow(model$x_for_y_unobs)
  )
}

generate_C_matrix <- function(model) {
  stopifnot(is.matrix(model$f_matrix_nieobs))
  col_sums <- colSums(model$f_matrix_nieobs)
  return(matrix(col_sums, ncol = 1))
}
