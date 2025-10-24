#' Exponential tilting  (ET) engine for NMAR estimation
#'
#' Constructs a configuration for the exponential tilting estimator under
#' nonignorable nonresponse (NMAR)..
#' The estimator solves \eqn{S_2(\boldsymbol{\phi}, \hat{\boldsymbol{\gamma}}) = 0,} using nleqslv to apply EM algorithm.

#'
#' @param standardize logical; standardize predictors. Default \code{TRUE}.
#' @param on_failure character; \code{"return"} or \code{"error"} on solver failure
#' @param variance_method character; one of \code{"delta"}, \code{"bootstrap"}, or \code{"none"}.
#' @param bootstrap_reps integer; number of bootstrap replicates when
#'   \code{variance_method = "bootstrap"}.
#' @param supress_warnings Logical; suppress variance-related warnings.
#' @param auxiliary_means Optional named numeric vector of population moments for
#'   auxiliary covariates.
#' @param control Optional list of additional solver controls.
#' @param family character; response model family, either \code{"logit"} or \code{"probit"},
#'   or a family object created by \code{logit_family()} / \code{probit_family()}.
#' @param y_dens Outcome density model (`"auto"`, `"normal"`, `"lognormal"`, or `"exponential"`).
#' @param min_iter Minimum number of solver iterations. (Will be migrated to
#'   `control` in a future release.)
#' @param max_iter Maximum number of solver iterations. (Will be migrated to
#'   `control` in a future release.)
#' @param optim_method Solver used for the response model (`"Newton"` or
#'   `"Broyden"`). (Will be migrated to `control` in a future release.)
#' @param tol_value Convergence tolerance for the optimisation routine. (Will be
#'   migrated to `control` in a future release.)
#' @details
#' The method is a robust Propensity-Score Adjustment (PSA) approach for Not Missing at Random (NMAR).
#' It uses Maximum Likelihood Estimation (MLE), basing the likelihood on the observed part of the sample (\eqn{f(\boldsymbol{Y}_i | \delta_i = 1, \boldsymbol{X}_i)}), making it robust against outcome model misspecification.
#' The propensity score is estimated by assuming an instrumental variable*(\eqn{X_2}) that is independent of the response status given other covariates and the study variable.
#' Estimator calculates fracional imputation weights weights \eqn{w_i}
#' The final estimator is a weighted average, where the weights are the inverse of the estimated response probabilities \eqn{\hat{\pi}_i}, satisfying the estimating equation:
#' \deqn{
#' \sum_{i \in \mathcal{R}} \frac{\boldsymbol{g}(\boldsymbol{Y}_i, \boldsymbol{X}_i ; \boldsymbol{\theta})}{\hat{\pi}_i} = 0,
#' }
#' where \eqn{\mathcal{R}} is the set of observed respondents.

#' @return An engine object of class \code{c("nmar_engine_el","nmar_engine")}.
#'   This is a configuration list; it is not a fit. Pass it to \link{nmar}.
#'
#' @references
#' Minsun Kim Riddles, Jae Kwang Kim, Jongho Im
#' A Propensity-score-adjustment Method for Nonignorable Nonresponse
#' \emph{Journal of Survey Statistics and Methodology}, Volume 4, Issue 2, June 2016, Pages 215–245.
#'
#' @examples
#' \donttest{
#' generate_test_data <- function(n_rows = 500, n_cols = 1, case = 1, x_var = 0.5, eps_var = 0.9, a = 0.8, b = -0.2) {
#' # Generate X variables - fixed to match comparison
#'   X <- as.data.frame(replicate(n_cols, rnorm(n_rows, 0, sqrt(x_var))))
#'   colnames(X) <- paste0("x", 1:n_cols)
#'
#' # Generate Y - fixed coefficients to match comparison
#'   eps <- rnorm(n_rows, 0, sqrt(eps_var))
#'   if (case == 1) {
#' # Use fixed coefficient of 1 for all x variables to match: y = -1 + x1 + epsilon
#'     X$Y <- as.vector(-1 + as.matrix(X) %*% rep(1, n_cols) + eps)
#'   }
#'   else if (case == 2) {
#'     X$Y <- -2 + 0.5 * exp(as.matrix(X) %*% rep(1, n_cols)) + eps
#'   }
#'   else if (case == 3) {
#'     X$Y <- -1 + sin(2 * as.matrix(X) %*% rep(1, n_cols)) + eps
#'   }
#'   else if (case == 4) {
#'     X$Y <- -1 + 0.4 * as.matrix(X)^3 %*% rep(1, n_cols) + eps
#'   }
#'
#'   Y_original <- X$Y
#'
#' # Missingness mechanism - identical to comparison
#'   pi_obs <- 1 / (1 + exp(-(a + b * X$Y)))
#'
#' # Create missing values
#'   mask <- runif(nrow(X)) > pi_obs
#'   mask[1] <- FALSE # Ensure at least one observation is not missing
#'   X$Y[mask] <- NA
#'
#'
#'   return(list(X = X, Y_original = Y_original))
#' }
#' res_test_data <- generate_test_data(n_rows = 500, n_cols = 1, case = 1)
#' x <- res_test_data$X
#'
#' exptilt_config <- exptilt_engine(
#'   y_dens = 'normal',
#'   min_iter = 3,
#'   max_iter = 10,
#'   tol_value = 0.01,
#'   standardize = FALSE,
#'   family = 'logit', # or logit
#'   bootstrap_reps = 50
#' )
#' formula = Y ~ x1
#' res <- nmar(formula = formula, data = x, engine = exptilt_config, response_predictors = NULL)
#' summary(res)
#' }
#' @keywords engine
#' @export
exptilt_engine <- function(
    standardize = FALSE,
    on_failure = c("return", "error"),
    variance_method = c("delta", "bootstrap"),
    bootstrap_reps = 10,
    supress_warnings = FALSE,
    auxiliary_means = NULL,
    control = list(),
    family = c("logit", "probit"),
    y_dens = c("auto", "normal", "lognormal", "exponential"),
    min_iter = 10, # TODO move to control
    max_iter = 100, # TODO move to control
    optim_method = c("Newton", "Broyden"), # TODO move to control
    tol_value = 0.1 # TODO move to control
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
  validator$assert_choice(y_dens, choices = c("auto", "normal", "lognormal", "exponential"), name = "y_dens")
  validator$assert_choice(variance_method, choices = c("delta", "bootstrap"), name = "variance_method")
  validator$assert_positive_integer(min_iter, name = "min_iter")
  validator$assert_positive_integer(max_iter, name = "max_iter")
  validator$assert_choice(optim_method, choices = c("Newton", "Broyden"), name = "optim_method")
  validator$assert_number(tol_value, name = "tol_value", min = 0, max = Inf)
  if (min_iter > max_iter) {
    stop("min_iter cannot be greater than max_iter.")
  }

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
