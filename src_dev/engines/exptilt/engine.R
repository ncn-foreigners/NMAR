#' Exponential tilting engine for NMAR estimation
#'
#' Build a configuration object for the exponential tilting (ET) estimator under
#' not-missing-at-random (NMAR) missingness. The configuration controls the
#' missingness link, outcome density model, scaling, variance computation, and
#' solver behaviour. Pass the resulting engine to [nmar()] to fit the ET
#' estimator.
#'
#' @param standardize Logical; standardise covariates prior to optimisation.
#' @param on_failure One of `"return"` or `"error"`; governs what happens if the
#'   solver fails to converge.
#' @param variance_method Variance estimator to use (`"delta"` or `"bootstrap"`).
#' @param bootstrap_reps Integer number of bootstrap replications when
#'   `variance_method = "bootstrap"`.
#' @param supress_warnings Logical; suppress variance-related warnings.
#' @param auxiliary_means Optional named numeric vector of population moments for
#'   auxiliary covariates.
#' @param control Optional list of additional solver controls.
#' @param family Missingness link function (`"logit"` or `"probit"`).
#' @param y_dens Outcome density model (`"auto"`, `"normal"`, or `"gamma"`).
#' @param min_iter Minimum number of solver iterations. (Will be migrated to
#'   `control` in a future release.)
#' @param max_iter Maximum number of solver iterations. (Will be migrated to
#'   `control` in a future release.)
#' @param optim_method Solver used for the response model (`"Newton"` or
#'   `"Broyden"`). (Will be migrated to `control` in a future release.)
#' @param tol_value Convergence tolerance for the optimisation routine. (Will be
#'   migrated to `control` in a future release.)
#'
#' @return A list of class `c("nmar_engine_exptilt", "nmar_engine")`.
#'
#' @export
exptilt_engine <- function(
    standardize = TRUE,
    on_failure = c("return", "error"),
    variance_method = c("delta", "bootstrap"),
    bootstrap_reps = 10,
    supress_warnings = FALSE,
    auxiliary_means = NULL,
    control = list(),
    family = c("logit", "probit"),
    y_dens = c("auto", "normal", "gamma"),
    min_iter = 10,  # TODO move to control
    max_iter = 100, # TODO move to control
    optim_method = c("Newton", "Broyden"), # TODO move to control
    tol_value = 1e-5 # TODO move to control
) {
  on_failure <- match.arg(on_failure)
  variance_method <- match.arg(variance_method)
  family <- match.arg(family)
  y_dens <- match.arg(y_dens)
  optim_method <- match.arg(optim_method)

  validator$assert_logical(standardize, name = "standardize")
  validator$assert_choice(on_failure, choices = c("return", "error"), name = "on_failure")
  validator$assert_positive_integer(bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  validator$assert_logical(supress_warnings, name = "supress_warnings")
  validator$assert_choice(family, choices = c("logit", "probit"), name = "family")
  validator$assert_choice(y_dens, choices = c("auto", "normal", "gamma"), name = "y_dens")
  validator$assert_choice(variance_method, choices = c("delta", "bootstrap"), name = "variance_method")
  validator$assert_positive_integer(min_iter, name = "min_iter")
  validator$assert_positive_integer(max_iter, name = "max_iter")
  validator$assert_choice(optim_method, choices = c("Newton", "Broyden"), name = "optim_method")
  validator$assert_number(tol_value, name = "tol_value", min = 0, max = Inf)

  engine <- list(
    standardize = standardize,
    on_failure = on_failure,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    supress_warnings = supress_warnings,
    auxiliary_means = auxiliary_means,
    control = control,
    prob_model_type = family,
    y_dens = y_dens,
    tol_value = tol_value,
    min_iter = min_iter,
    max_iter = max_iter,
    optim_method = optim_method
  )
  class(engine) <- c("nmar_engine_exptilt", "nmar_engine")
  engine
}
