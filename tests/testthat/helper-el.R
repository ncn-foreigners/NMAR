make_iid_nmar <- function(n = 200, alpha = 0.4, include_z = FALSE, seed = 1L) {
  set.seed(seed)
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 2 + 0.5 * X + Z
  p <- stats::plogis(-0.6 + alpha * scale(Y)[, 1] + if (include_z) 0.3 * Z else 0)
  R <- stats::runif(n) < p
  df <- data.frame(Y_miss = Y, X = X)
  if (include_z) df$Z <- Z
  df$Y_miss[!R] <- NA_real_
  df
}

make_el_test_data <- function(n = 100) {
  set.seed(123)
  X <- rnorm(n)
  Z <- rnorm(n)
  Y <- 1 + 0.5 * X + 0.3 * Z + rnorm(n)
  p <- plogis(-0.3 + 0.4 * scale(Y)[, 1] + 0.2 * Z)
  R <- runif(n) < p
  df <- data.frame(Y_miss = Y, X = X, Z = Z)
  df$Y_miss[!R] <- NA_real_
  df
}

make_engine <- function(variance_method = c("delta", "bootstrap", "none"),
                        family = c("logit", "probit"),
                        standardize = TRUE,
                        auxiliary_means = NULL,
                        control = list(),
                        trim_cap = Inf,
                        bootstrap_reps = 50) {
  variance_method <- match.arg(variance_method)
  family <- match.arg(family)
  NMAR::el_engine(
    variance_method = variance_method,
    family = family,
    standardize = standardize,
    auxiliary_means = auxiliary_means,
    control = control,
    trim_cap = trim_cap,
    bootstrap_reps = bootstrap_reps
  )
}
