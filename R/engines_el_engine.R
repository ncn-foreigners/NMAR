#' Empirical likelihood (EL) engine for NMAR
#'
#' Constructs an engine specification for the empirical likelihood (EL)
#' estimator of a full-data mean under nonignorable nonresponse (NMAR).
#'
#' The implementation follows Qin, Leung, and Shao (2002): the response
#' mechanism is modeled as \eqn{w(y, x; \beta) = P(R = 1 \mid Y = y, X = x)} and
#' the joint law of \eqn{(Y, X)} is represented nonparametrically by respondent
#' masses that satisfy empirical likelihood constraints. The mean is estimated
#' as a respondent weighted mean with weights proportional to
#' \eqn{\tilde w_i = a_i / D_i(\beta, W, \lambda)}, where \eqn{a_i} are base
#' weights (\eqn{a_i \equiv 1} for IID data and \eqn{a_i = d_i} for survey
#' designs) and \eqn{D_i} is the EL denominator.
#'
#' For \code{data.frame} inputs the estimator solves the Qin-Leung-Shao (QLS)
#' estimating equations for \eqn{(\beta, W, \lambda_x)} with \eqn{W} reparameterized
#' as \eqn{z = \operatorname{logit}(W)}, and profiles out the response multiplier
#' \eqn{\lambda_W} using the closed-form QLS identity (their Eq. 10). For
#' \code{survey.design} inputs the estimator uses a design-weighted analogue
#' (Chen and Sitter 1999; Wu 2005) with an explicit \eqn{\lambda_W} and an
#' additional linkage equation involving the nonrespondent design-weight total
#' \eqn{T_0}.
#'
#' Numerical stability:
#' \itemize{
#'   \item \eqn{W} is optimized on the logit scale so \eqn{0 < W < 1}.
#'   \item The response-model linear predictor is capped and EL denominators
#'     \eqn{D_i} are floored at a small positive value; the analytic Jacobian is
#'     consistent with this guard via an active-set mask.
#'   \item Optional trimming (\code{trim_cap}) is applied only after solving, to
#'     the unnormalized masses \eqn{\tilde w_i = a_i/D_i}; this changes the
#'     returned weights and therefore the point estimate.
#' }
#'
#' @param family Missingness (response) model family. Either \code{"logit"}
#'   (default) or \code{"probit"}, or a custom family object: a list with
#'   components \code{name}, \code{linkinv}, \code{mu.eta}, \code{score_eta}, and
#'   optionally \code{d2mu.deta2}. When \code{d2mu.deta2} is absent the solver
#'   uses Broyden/numeric Jacobians.
#' @param standardize logical; standardize predictors. Default \code{TRUE}.
#' @param trim_cap numeric; cap for EL weights (\code{Inf} = no trimming).
#' @param on_failure character; \code{"return"} or \code{"error"} on solver failure.
#' @param variance_method character; one of \code{"bootstrap"} or \code{"none"}.
#' @param bootstrap_reps integer; number of bootstrap replicates when
#'   \code{variance_method = "bootstrap"}.
#' @param auxiliary_means named numeric vector; population means for auxiliary
#'   design columns. Names must match the materialized model.matrix column names
#'   on the first RHS (after formula expansion), e.g., factor indicator columns
#'   created by \code{model.matrix()} or transformed terms like \code{I(X^2)}.
#'   Auxiliary intercepts are always dropped automatically, so do not supply
#'   \code{(Intercept)}. If \code{NULL} (default)
#'   and the outcome contains at least one \code{NA}, auxiliary means are estimated
#'   from the full input (including nonrespondents): IID uses unweighted column
#'   means of the auxiliary design; survey designs use the design-weighted means
#'   based on \code{weights(design)}. This corresponds to the QLS case where
#'   \eqn{\mu_x} is replaced by \eqn{\bar X} (the full-sample mean) when auxiliary
#'   variables are observed for all sampled units.
#' @param control Optional solver configuration forwarded to
#'   \code{nleqslv::nleqslv()}. Provide a single list that may include solver
#'   tolerances (e.g., \code{xtol}, \code{ftol}, \code{maxit}) and, optionally,
#'   top-level entries \code{global} and \code{xscalm} for globalization and
#'   scaling. Example:
#'   \code{control = list(maxit = 500, xtol = 1e-10, ftol = 1e-10, global = "qline", xscalm = "auto")}.
#' @param strata_augmentation logical; when \code{TRUE} (default), survey designs
#'   with an identifiable strata structure are augmented with stratum indicators
#'   and corresponding population shares in the auxiliary block (Wu-style
#'   strata augmentation). Has no effect for \code{data.frame} inputs or
#'   survey designs without strata.
#' @param n_total numeric; optional when supplying respondents-only data (no \code{NA} in the
#'   outcome). For \code{data.frame} inputs, set to the total number of sampled units
#'   before filtering to respondents. For \code{survey.design} inputs, set to the total
#'   design-weight total on the same analysis scale as \code{weights(design)}
#'   (default \code{sum(weights(design))}). If omitted and the outcome contains no NAs,
#'   the estimator errors, requesting \code{n_total}.
#' @param start list; optional starting point for the solver. Fields:
#'   \itemize{
#'     \item \code{beta}: named numeric vector of missingness-model coefficients on the
#'       original (unscaled) scale, including \code{(Intercept)}.
#'     \item \code{W} or \code{z}: starting value for population response rate (\code{0 < W < 1})
#'       or its logit (\code{z}). If both are provided, \code{z} takes precedence.
#'     \item \code{lambda}: named numeric vector of auxiliary multipliers on the original
#'       scale (names must match auxiliary design columns; no intercept). Values
#'       are mapped to the scaled space internally.
#'   }
#'
#' @return
#' A list of class \code{"nmar_engine_el"} (also inheriting from \code{"nmar_engine"})
#' containing configuration fields to be supplied to \code{nmar()}. Users rarely
#' access fields directly; instead, pass the engine to \code{nmar()} together with a
#' formula and data.
#'
#' @details
#' \strong{Formula syntax and data constraints}: \code{nmar()} accepts a
#' partitioned right-hand side \code{y_miss ~ auxiliaries | response_only}. Variables
#' left of \code{|} enter auxiliary moment constraints; variables right of \code{|}
#' enter only the response model. The outcome (LHS) is always included as a
#' response-model predictor through the evaluated LHS expression; explicit use of
#' the outcome on the RHS is rejected. The response model always includes an
#' intercept; the auxiliary block never includes an intercept.
#'
#' To include a covariate in both the auxiliary constraints and the response
#' model, repeat it on both sides, e.g. \code{y_miss ~ X | X}.
#'
#' \strong{Auxiliary means}: If \code{auxiliary_means = NULL} (default) and the
#' outcome contains at least one \code{NA}, auxiliary means are estimated from the
#' full input and used as \eqn{\bar X} in the QLS constraints. For respondents-only
#' data (no \code{NA} in the outcome), \code{n_total} must be supplied; and if the
#' auxiliary RHS is non-empty, \code{auxiliary_means} must also be supplied.
#' When \code{standardize = TRUE}, supply \code{auxiliary_means} on the original
#' data scale; the engine applies the same standardization internally.
#'
#' \strong{Survey scale}: For \code{survey.design} inputs, \code{n_total} (if
#' provided) must be on the same analysis scale as \code{weights(design)}. The
#' default is \code{sum(weights(design))}.
#'
#' \strong{Convergence and identification}: the stacked EL system can have
#' multiple solutions. Adding response-only predictors (variables to the right
#' of \code{|}) can make the problem sensitive to starting values. Inspect
#' diagnostics such as \code{jacobian_condition_number} and consider supplying
#' \code{start = list(beta = ..., W = ...)} when needed.
#'
#' \strong{Variance}: The EL engine supports bootstrap standard errors via
#' \code{variance_method = "bootstrap"} or can skip variance with
#' \code{variance_method = "none"}.
#' Set a seed for reproducible bootstrap results.
#'
#' Bootstrap requires suggested packages: for IID resampling it requires
#' \code{future.apply} (and \code{future}); for survey replicate-weight bootstrap
#' it requires \code{survey} and \code{svrep}.
#'
#' @references
#' Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical
#' Association, 97(457), 193-200. \doi{10.1198/016214502753479338}
#'
#' Chen, J., and Sitter, R. R. (1999). A pseudo empirical likelihood approach for
#' the effective use of auxiliary information in complex surveys. Statistica Sinica,
#' 9, 385-406.
#'
#' Wu, C. (2005). Algorithms and R codes for the pseudo empirical likelihood method
#' in survey sampling. Survey Methodology, 31(2), 239-243.
#' @keywords engine
#' @export
#' @seealso \code{\link{nmar}}, \code{\link{weights.nmar_result}}, \code{\link{summary.nmar_result}}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' X <- rnorm(n)
#' Y <- 2 + 0.5 * X + rnorm(n)
#' p <- plogis(-0.7 + 0.4 * scale(Y)[, 1])
#' R <- runif(n) < p
#' if (all(R)) R[1] <- FALSE
#' df <- data.frame(Y_miss = Y, X = X)
#' df$Y_miss[!R] <- NA_real_
#'
#' # Estimate auxiliary mean from full data (QLS "use Xbar" case)
#' eng <- el_engine(auxiliary_means = NULL, variance_method = "none")
#'
#' # Put X in both the auxiliary block and the response model (QLS-like)
#' fit <- nmar(Y_miss ~ X | X, data = df, engine = eng)
#' summary(fit)
#'
#' \donttest{
#' # Response-only predictors can be placed to the right of |:
#' Z <- rnorm(n)
#' df2 <- data.frame(Y_miss = Y, X = X, Z = Z)
#' df2$Y_miss[!R] <- NA_real_
#' eng2 <- el_engine(auxiliary_means = NULL, variance_method = "none")
#' fit2 <- nmar(Y_miss ~ X | Z, data = df2, engine = eng2)
#' print(fit2)
#'
#' # Survey design usage
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   des <- survey::svydesign(ids = ~1, weights = ~1, data = df)
#'   eng3 <- el_engine(auxiliary_means = NULL, variance_method = "none")
#'   fit3 <- nmar(Y_miss ~ X, data = des, engine = eng3)
#'   summary(fit3)
#' }
#' }
el_engine <- function(
    standardize = TRUE,
    trim_cap = Inf,
    on_failure = c("return", "error"),
    variance_method = c("bootstrap", "none"),
    bootstrap_reps = 500,
    auxiliary_means = NULL,
    control = list(),
    strata_augmentation = TRUE,
    n_total = NULL,
    start = NULL,
    family = c("logit", "probit")) {
  on_failure <- match.arg(on_failure)
  if (is.null(variance_method)) variance_method <- "none"
  variance_method <- match.arg(variance_method)
  if (is.character(family)) family <- match.arg(family)
  if (is.character(family)) {
    family <- switch(family,
      logit = logit_family(),
      probit = probit_family()
    )
  }

  engine <- list(
    standardize = standardize,
    trim_cap = trim_cap,
    on_failure = on_failure,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    auxiliary_means = auxiliary_means,
    control = control,
    strata_augmentation = strata_augmentation,
    n_total = n_total,
    start = start,
    family = family
  )
# Validate and coerce unsupported variance modes upfront
  validate_nmar_engine_el(engine)
  new_nmar_engine_el(engine)
}

#' Construct EL Engine Object
#' @keywords internal
new_nmar_engine_el <- function(engine) {
  stopifnot(is.list(engine))
  class(engine) <- c("nmar_engine_el", "nmar_engine")
  engine
}

#' Validate EL Engine Settings
#' @keywords internal
validate_nmar_engine_el <- function(engine) {
  validator_assert_list(engine, name = "engine")
  validator_assert_logical(engine$standardize, name = "standardize")

  validator_assert_positive_number(engine$trim_cap, name = "trim_cap", allow_infinite = TRUE)
  validator_assert_positive_integer(engine$bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  validator_assert_choice(engine$on_failure, choices = c("return", "error"), name = "on_failure")
  validator_assert_choice(engine$variance_method, choices = c("bootstrap", "none"), name = "variance_method")

  validator_assert_named_numeric(engine$auxiliary_means, name = "auxiliary_means", allow_null = TRUE)
  validator_assert_list(engine$control, name = "control")
  validator_assert_logical(engine$strata_augmentation, name = "strata_augmentation")

  if (!is.null(engine$n_total)) validator_assert_positive_number(engine$n_total, name = "n_total", allow_infinite = FALSE)

  fam <- engine$family
  if (!is.list(fam) || is.null(fam$name) || !is.function(fam$linkinv) || !is.function(fam$mu.eta) || !is.function(fam$score_eta)) {
    stop("family must be 'logit'/'probit' or a valid family object with linkinv, mu.eta, score_eta")
  }

  engine
}
