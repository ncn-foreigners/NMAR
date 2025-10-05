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

make_engine <- function(variance_method = c("delta", "bootstrap", "none"),
                        family = c("logit", "probit"),
                        standardize = TRUE,
                        auxiliary_means = NULL,
                        solver_method = c("auto", "newton", "broyden"),
                        solver_args = list(),
                        control = list(),
                        trim_cap = Inf,
                        bootstrap_reps = 50,
                        variance_ridge = FALSE,
                        variance_pseudoinverse = FALSE,
                        suppress_warnings = TRUE) {
  variance_method <- match.arg(variance_method)
  family <- match.arg(family)
  solver_method <- match.arg(solver_method)
  NMAR::el_engine(
    variance_method = variance_method,
    family = family,
    standardize = standardize,
    auxiliary_means = auxiliary_means,
    solver_method = solver_method,
    solver_args = solver_args,
    control = control,
    trim_cap = trim_cap,
    bootstrap_reps = bootstrap_reps,
    variance_ridge = variance_ridge,
    variance_pseudoinverse = variance_pseudoinverse,
    suppress_warnings = suppress_warnings
  )
}
