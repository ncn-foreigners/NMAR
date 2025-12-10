#' Exponential tilting  (ET) engine for NMAR estimation
#'
#' Constructs a configuration for the exponential tilting estimator under
#' nonignorable nonresponse (NMAR).
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
#' @param control Named list of control parameters passed to \code{nleqslv::nleqslv}.
#'   Common parameters include:
#'   \itemize{
#'     \item \code{maxit}: Maximum number of iterations (default: 100)
#'     \item \code{method}: Solver method - \code{"Newton"} or \code{"Broyden"} (default: \code{"Newton"})
#'     \item \code{global}: Global strategy - \code{"dbldog"}, \code{"pwldog"}, \code{"qline"}, \code{"gline"}, \code{"hook"}, or \code{"none"} (default: \code{"dbldog"})
#'     \item \code{xtol}: Tolerance for relative error in solution (default: 1e-8)
#'     \item \code{ftol}: Tolerance for function value (default: 1e-8)
#'     \item \code{btol}: Tolerance for backtracking (default: 0.01)
#'     \item \code{allowSingular}: Allow singular Jacobians (default: \code{TRUE})
#'   }
#'   See \code{?nleqslv::nleqslv} for full details.
#' @param stopping_threshold Numeric; early stopping threshold. If the maximum absolute value
#'   of the score function falls below this threshold, the algorithm stops early (default: 1).
#' @param family character; response model family, either \code{"logit"} or \code{"probit"},
#'   or a family object created by \code{logit_family()} / \code{probit_family()}.
#' @param y_dens Outcome density model (\code{"auto"}, \code{"normal"}, \code{"lognormal"}, or \code{"exponential"}).
#' @param sample_size Integer; maximum sample size for stratified random sampling (default: 2000).
#'   When the dataset exceeds this size, a stratified random sample is drawn to optimize memory usage.
#'   The sampling preserves the ratio of respondents to non-respondents in the original data.
#' @details
#' The method is a robust Propensity-Score Adjustment (PSA) approach for Not Missing at Random (NMAR).
#' It uses Maximum Likelihood Estimation (MLE), basing the likelihood on the observed part of the sample (\eqn{f(\boldsymbol{Y}_i | \delta_i = 1, \boldsymbol{X}_i)}), making it robust against outcome model misspecification.
#' The propensity score is estimated by assuming an instrumental variable (\eqn{X_2}) that is independent of the response status given other covariates and the study variable.
#' Estimator calculates fractional imputation weights \eqn{w_i}.
#' The final estimator is a weighted average, where the weights are the inverse of the estimated response probabilities \eqn{\hat{\pi}_i}, satisfying the estimating equation:
#' \deqn{
#' \sum_{i \in \mathcal{R}} \frac{\boldsymbol{g}(\boldsymbol{Y}_i, \boldsymbol{X}_i ; \boldsymbol{\theta})}{\hat{\pi}_i} = 0,
#' }
#' where \eqn{\mathcal{R}} is the set of observed respondents.

#' @return An engine object of class \code{c("nmar_engine_exptilt","nmar_engine")}.
#'   This is a configuration list; it is not a fit. Pass it to \link{nmar}.
#'
#' @references
#' Minsun Kim Riddles, Jae Kwang Kim, Jongho Im
#' A Propensity-score-adjustment Method for Nonignorable Nonresponse
#' \emph{Journal of Survey Statistics and Methodology}, Volume 4, Issue 2, June 2016, Pages 215–245.
#'
#' @examples
#' \donttest{
#' generate_test_data <- function(
#'   n_rows = 500,
#'   n_cols = 1,
#'   case = 1,
#'   x_var = 0.5,
#'   eps_var = 0.9,
#'   a = 0.8,
#'   b = -0.2
#' ) {
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
#'   return(list(X = X, Y_original = Y_original))
#' }
#' res_test_data <- generate_test_data(n_rows = 500, n_cols = 1, case = 1)
#' x <- res_test_data$X
#'
#' exptilt_config <- exptilt_engine(
#'   y_dens = 'normal',
#'   control = list(maxit = 10),
#'   stopping_threshold = 0.1,
#'   standardize = FALSE,
#'   family = 'logit',
#'   bootstrap_reps = 5
#' )
#' formula = Y ~ x1
#' res <- nmar(formula = formula, data = x, engine = exptilt_config, trace_level = 1)
#' summary(res)
#' }
#' @keywords engine
#' @export
exptilt_engine <- function(
    standardize = FALSE,
    on_failure = c("return", "error"),
    variance_method = c("delta", "bootstrap", 'none'),
    bootstrap_reps = 10,
    supress_warnings = FALSE,
    auxiliary_means = NULL,
    control = list(),
    family = c("logit", "probit"),
    y_dens = c("normal", "lognormal", "exponential"),
    stopping_threshold = 1,
    sample_size = 2000
    ) {
  on_failure <- match.arg(on_failure)
  variance_method <- match.arg(variance_method)
  family <- match.arg(family)
  y_dens <- match.arg(y_dens)

  validator_assert_logical(standardize, name = "standardize")
  validator_assert_choice(on_failure, choices = c("return", "error"), name = "on_failure")
  validator_assert_positive_integer(bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  validator_assert_logical(supress_warnings, name = "supress_warnings")
  validator_assert_choice(family, choices = c("logit", "probit"), name = "family")
  validator_assert_choice(y_dens, choices = c("normal", "lognormal", "exponential"), name = "y_dens")
  validator_assert_choice(variance_method, choices = c("delta", "bootstrap", 'none'), name = "variance_method")
  validator_assert_number(stopping_threshold, name = "stopping_threshold", min = 0, max = Inf)
  validator_assert_list(control, name = "control")
  validator_assert_positive_integer(sample_size, name = "sample_size", is.finite = TRUE)

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
    stopping_threshold = stopping_threshold,
    sample_size = sample_size
  )

  class(engine) <- c("nmar_engine_exptilt", "nmar_engine")
  engine

}
