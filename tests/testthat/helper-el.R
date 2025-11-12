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

prepare_el_inputs <- function(formula, data, require_na = TRUE, auxiliary_means = NULL) {
  design <- NMAR:::el_construct_design(formula, data, require_na)
  delta <- NMAR:::el_make_delta_column(data, design$outcome_var, design$respondent_mask)
  list(
    data = delta$data,
    outcome_var = design$outcome_var,
    delta_name = delta$delta_name,
    respondent_mask = design$respondent_mask,
    response_matrix = design$missingness_model_matrix,
    response_outcome = design$y_obs,
    auxiliary_matrix = design$aux_mm_full[design$respondent_mask, , drop = FALSE],
    auxiliary_matrix_full = design$aux_mm_full,
    has_aux = ncol(design$aux_mm_full) > 0
  )
}
