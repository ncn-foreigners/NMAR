#' Induced logistic regression engine for NMAR
#'
#' Constructs an engine configuration for the induced logistic regression
#' estimator of Li, Qin and Liu. Use it with
#' \code{\link{nmar}} via \code{engine = induced_logit_engine(...)}.
#'
#' @details
#' This engine supports IID \code{data.frame} inputs and \code{survey.design}
#' inputs from the suggested \pkg{survey} package. The outcome is assumed to be
#' missing for nonrespondents (\code{NA} indicates nonresponse). All covariates
#' used in either model must be fully observed for all sampled units.
#'
#' The estimator fits:
#' \enumerate{
#' \item An outcome regression on respondents to estimate
#' \eqn{E(Y \mid X=x, R=1)}.
#' \item An induced logistic regression for the response indicator \eqn{R} on
#' \eqn{x_1} and \eqn{\hat{\mu}(x)}.
#' }
#' It then estimates the population mean \eqn{\tau = E(Y)} using stabilized
#' exponential tilt moment calculations.
#'
#' \strong{Assumptions:}
#' \itemize{
#' \item This method is not doubly robust. Misspecifying either the outcome
#' regression \eqn{\mu(x;\xi)} or the response model can bias \eqn{\hat\tau}.
#' \item Step 1 is implemented as a linear regression on
#' respondents (IID: \code{stats::lm.fit}; survey: \code{survey::svyglm}).
#' \item Identification and numerical stability can be weak when
#' \eqn{\hat{\mu}(x)} is (nearly) collinear with the \eqn{x_1} design matrix. The
#' implementation fails fast on rank deficiency and extreme ill-conditioning.
#' \item The survey-design path is a weight-based pseudo-likelihood extension of
#' the IID estimator and is not derived in Li, Qin and Liu. Proceed with caution
#' for designs that require probability- or FPC-based handling beyond analysis
#' weights.
#' }
#'
#' \strong{Formula interface:} \code{y_miss ~ mu_covariates | x1_covariates}.
#' The first block specifies the outcome regression; the second specifies the
#' missingness predictors \eqn{x_1}. If \code{|} is omitted, the missingness
#' block defaults to \code{| 1}.
#'
#' \strong{Parameterization:} The paper writes the response model as
#' \deqn{\Pr(R=1 \mid x, y) = 1 / (1 + \exp\{\alpha_0 + x_1^T \beta + \gamma y\}).}
#' This implementation fits a standard logit GLM for \eqn{Pr(R=1 \mid x)} with
#' linear predictor \eqn{b_0 + b_{x_1}^T x_1 + b_{\mu} \hat{\mu}(x)}. Paper
#' parameters correspond to the negation of fitted GLM coefficients; diagnostics
#' report \code{gamma_hat_paper = -b_mu}, where \eqn{b_\mu} is the fitted GLM
#' coefficient on \eqn{\hat{\mu}(x)}.
#'
#' \strong{Survey designs:} The paper is an IID development. For
#' \code{survey.design} inputs, this engine uses design-weighted estimating
#' equations via \code{survey::svyglm}. Calibrated and post-stratified designs
#' are rejected. For other design features that go beyond analysis weights,
#' \code{survey_design_policy = "strict"} fails fast and \code{"warn"} proceeds
#' with a warning.
#'
#' \strong{Variance}: The induced logit engine supports bootstrap standard errors via
#' \code{variance_method = "bootstrap"} or can skip variance with
#' \code{variance_method = "none"}.
#'
#' Bootstrap uses no additional packages for IID resampling, and will run
#' sequentially by default. If the suggested \code{future.apply} package is
#' installed, IID bootstrap can use \code{future.apply::future_lapply()} according to
#' the user's \code{future::plan()} for parallel execution.
#' Bootstrap backend is controlled by the package option \code{nmar.bootstrap_apply}:
#' \describe{
#' \item{\code{"auto"}}{(default) Use \code{base::lapply()} unless the current
#' future plan has more than one worker, in which case use
#' \code{future.apply::future_lapply()} if available.}
#' \item{\code{"base"}}{Always use \code{base::lapply()}, even
#' if \code{future.apply} is installed.}
#' \item{\code{"future"}}{Always use \code{future.apply::future_lapply()}.}
#' }
#' For \code{survey.design} inputs, replicate-weight bootstrap requires the
#' suggested packages \code{survey} and \code{svrep}.
#'
#' \strong{Numerical stability:} When computing exponential tilt moment terms,
#' this implementation uses max-shift (log-sum-exp) stabilization to avoid
#' overflow and underflow. Setting \code{standardize = TRUE} additionally
#' centers and scales design matrices as a reparameterization.
#'
#' @param variance_method Character; one of \code{"none"} or \code{"bootstrap"}.
#' @param bootstrap_reps Integer; number of bootstrap replicates when
#' \code{variance_method = "bootstrap"}.
#' @param standardize Logical; if \code{TRUE}, internally standardizes (centers/scales)
#' the outcome-model and missingness-model design matrices for numerical stability.
#' This is implemented as a reparameterization; point estimates and fitted values
#' are invariant up to numerical tolerance.
#' @param control A list of control parameters for the missingness-model GLM fit.
#' This is forwarded to \code{stats::glm.control()}. Common entries include
#' \code{maxit} and \code{epsilon}. This affects only the missingness-model
#' estimation step.
#' @param on_failure Character; \code{"return"} (default) returns a non-converged
#' result object on error, \code{"error"} rethrows the error.
#' @param survey_design_policy Character; one of \code{"strict"} or \code{"warn"}.
#' Controls how the survey path handles design features that are outside the
#' paper and beyond the weight-only pseudo-likelihood extension used by this
#' implementation (for example PPS/FPC or multistage probabilities).
#' \code{"strict"} fails fast; \code{"warn"} proceeds with a warning.
#' @param keep_fits Logical; if \code{TRUE}, stores internal fit objects in
#' \code{result$extra$raw} for debugging. Default \code{FALSE} keeps result objects
#' lighter-weight.
#'
#' @return An engine object of class \code{c("nmar_engine_induced_logit","nmar_engine")}.
#' @references
#' Li, P., Qin, J., and Liu, Y. (2023). Instability of Inverse Probability
#' Weighting Methods and a Remedy for Nonignorable Missing Data. Biometrics,
#' 79(4), 3215-3226. \doi{10.1111/biom.13881}
#' @seealso \code{\link{nmar}}, \code{\link{summary.nmar_result}}, \code{\link{coef.nmar_result}},
#' \code{vignette("tutorial_induced_logit")}
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 0.5 * x1 - 0.25 * x2 + rnorm(n)
#'
#' # Nonignorable missingness: response depends on the possibly missing outcome
#' p <- plogis(-0.4 + 0.15 * y + 0.2 * x1)
#' r <- rbinom(n, 1, p)
#' if (all(r == 1L)) r[1] <- 0L
#' if (all(r == 0L)) r[1] <- 1L
#'
#' df <- data.frame(y_miss = ifelse(r == 1L, y, NA_real_), x1 = x1, x2 = x2)
#'
#' eng <- induced_logit_engine(variance_method = "none", standardize = TRUE, on_failure = "error")
#' fit <- nmar(y_miss ~ x1 + x2 | x1, data = df, engine = eng)
#' summary(fit)
#'
#' # Survey design usage
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   des <- survey::svydesign(ids = ~1, weights = ~1, data = df)
#'   fit_svy <- nmar(y_miss ~ x1 + x2 | x1, data = des, engine = eng)
#'   summary(fit_svy)
#' }
#' @keywords engine
#' @export
induced_logit_engine <- function(
    variance_method = c("none", "bootstrap"),
    bootstrap_reps = 500,
    standardize = FALSE,
    control = list(),
    on_failure = c("return", "error"),
    survey_design_policy = c("strict", "warn"),
    keep_fits = FALSE) {
  variance_method <- match.arg(variance_method)
  on_failure <- match.arg(on_failure)
  survey_design_policy <- match.arg(survey_design_policy)
  validator_assert_scalar_logical(standardize, name = "standardize")
  validator_assert_scalar_logical(keep_fits, name = "keep_fits")
  validator_assert_list(control, name = "control")

  engine <- list(
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    standardize = standardize,
    control = control,
    on_failure = on_failure,
    survey_design_policy = survey_design_policy,
    keep_fits = keep_fits
  )

  validate_nmar_engine_induced_logit(engine)
  new_nmar_engine_induced_logit(engine)
}

#' @keywords internal
#' @noRd
new_nmar_engine_induced_logit <- function(engine) {
  stopifnot(is.list(engine))
  class(engine) <- c("nmar_engine_induced_logit", "nmar_engine")
  engine
}

#' @keywords internal
#' @noRd
validate_nmar_engine_induced_logit <- function(engine) {
  validator_assert_list(engine, name = "engine")
  keep_fits <- if (is.null(engine$keep_fits)) FALSE else engine$keep_fits

  il_validate_induced_logit_opts(
    variance_method = engine$variance_method,
    bootstrap_reps = engine$bootstrap_reps,
    standardize = engine$standardize,
    control = engine$control %||% list(),
    on_failure = engine$on_failure,
    survey_design_policy = engine$survey_design_policy,
    keep_fits = keep_fits
  )
  engine
}
