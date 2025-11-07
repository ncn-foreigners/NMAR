#' Empirical likelihood (EL) engine for NMAR
#'
#' Constructs a configuration object for the empirical likelihood estimator under
#' nonignorable nonresponse (NMAR) with optional auxiliary moment constraints.
#' The estimator solves the stacked system in \eqn{\theta = (\beta, z, \lambda_x)}
#' with \eqn{z = \operatorname{logit}(W)} using a Newton method with analytic
#' Jacobian and globalization via \link[nleqslv]{nleqslv}. Numerical safeguards
#' (bounded linear predictor, link-inverse clipping, denominator floors, and
#' stable linear algebra) improve robustness. Pass the engine to \link{nmar}
#' together with a formula and data.
#'
#' @param family character; response model family, either \code{"logit"} or \code{"probit"},
#'   or a family object created by \code{logit_family()} / \code{probit_family()}.
#' @param standardize logical; standardize predictors. Default \code{TRUE}.
#' @param trim_cap numeric; cap for EL weights (\code{Inf} = no trimming).
#' @param on_failure character; \code{"return"} or \code{"error"} on solver failure.
#' @param variance_method character; one of \code{"delta"}, \code{"bootstrap"}, or \code{"none"}.
#' @param bootstrap_reps integer; number of bootstrap replicates when
#'   \code{variance_method = "bootstrap"}.
#' @param auxiliary_means named numeric vector; population means for auxiliaries
#'   (names must match the RHS of the outcome formula). Optional.
#' @param control list; optional solver control for \code{nleqslv::nleqslv()}.
#'   Recognized fields (defaults in parentheses):
#'   \itemize{
#'     \item Top-level: \code{global} = \code{"qline"} (quadratic line search) or one of
#'       \code{"dbldog"}, \code{"pwldog"}, \code{"cline"}, \code{"gline"}, \code{"hook"}, \code{"none"};
#'       \code{xscalm} = \code{"auto"} or \code{"fixed"}
#'     \item In \code{control=}: \code{xtol}, \code{ftol}, \code{btol}, \code{maxit}, \code{trace},
#'       \code{stepmax}, \code{delta}, \code{allowSing}
#'   }
#'   Unknown names are ignored. The method is Newton with an analytic Jacobian.
#' @param n_total numeric; optional when supplying respondents-only data (no \code{NA} in the
#'   outcome). For \code{data.frame} inputs, set to the total number of sampled units
#'   before filtering to respondents. For \code{survey.design} inputs, set to the total
#'   design weight or known population total. If omitted and the outcome contains no NAs,
#'   the estimator errors, requesting \code{n_total}.
#' @param start list; optional starting point for the solver. Fields:
#'   \itemize{
#'     \item \code{beta}: named numeric vector of response-model coefficients on the
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
#' Empirical likelihood assigns masses \eqn{m_i = d_i / D_i(\theta)} to respondents
#' with base weights \eqn{d_i}. Following Qin, Leung, and Shao (2002, JASA 97:193-200),
#' the denominator is
#' \deqn{D_i(\theta) = 1 + \lambda_W\{w_i(\beta) - W\} + X_{i\cdot}^{(c)}\,\lambda_x,}
#' where \eqn{w_i(\beta) = g(\eta_i)} with link-inverse \eqn{g}, \eqn{\eta_i = x_i^\top\beta},
#' \eqn{W} is the (unknown) response rate, and \eqn{X^{(c)} = X - \mu_X} centers
#' auxiliary columns at their population means. The multiplier for the response-rate
#' equation is
#' \deqn{\lambda_W = \{N_\mathrm{pop}/\sum_i d_i - 1\}/(1 - W),}
#' which ensures scale coherence for both IID data (\eqn{d_i \equiv 1}) and survey
#' designs (\eqn{d_i} are design weights). The estimating equations impose
#' \eqn{\sum_i m_i\{w_i(\beta) - W\} = 0},
#' \eqn{\sum_i m_i\,s_i(\beta) = 0} with \eqn{s_i(\beta) = \partial \log w_i / \partial \eta_i},
#' and, when present, \eqn{\sum_i m_i X_{i\cdot}^{(c)} = 0}.
#'
#' The response-model score used in both equations and Jacobian is the derivative
#' of the Bernoulli log-likelihood with respect to the linear predictor, i.e.
#' \code{mu.eta(eta) / linkinv(eta)} (logit: \code{1 - w}; probit: Mills ratio \code{phi/Phi}).
#' We apply a consistent guarding policy (cap \eqn{\eta}, clip \eqn{w}, floor denominators
#' with an "active" mask in the Jacobian) to ensure numerical stability and to make the
#' analytic Jacobian match the piecewise-smooth equations being solved. When
#' \code{variance_method = "delta"} is requested, the estimator returns \code{NA}
#' standard errors with a message; use \code{variance_method = "bootstrap"} for SEs.
#'
#' Survey designs are handled by replacing counts in QLS (2002) with design-weighted
#' totals. In particular, the response-rate multiplier generalizes to
#' \eqn{\lambda_W = \{N_\mathrm{pop}/\sum_{i:\,R_i=1} d_i - 1\}/(1 - W)}, which reduces
#' to the QLS expression when \eqn{d_i \equiv 1}. Solver configuration uses \code{nleqslv}
#' with an analytic Jacobian and line-search globalization. Defaults are
#' \code{global = "qline"} and \code{xscalm = "auto"}; users can override via \code{control}.
#' Invalid values are coerced to these defaults with a warning.
#'
#' \strong{Formula syntax}: \code{nmar()} supports a partitioned right-hand side
#' \code{y_miss ~ aux1 + aux2 | z1 + z2}. Variables left of \code{|} are auxiliaries
#' (used in EL moment constraints); variables right of \code{|} are response-model
#' predictors only. The outcome appears on the left-hand side and is included as a
#' response predictor by default.
#'
#' \strong{Weights in results}: Calling \code{weights()} on the returned \code{nmar_result}
#' gives respondent weights on either the probability scale (sum to 1) or the population
#' scale (sum to \eqn{N_\mathrm{pop}}). The reported masses come from the empirical
#' likelihood construction \eqn{a_i/D_i(\theta)} and are normalized in \code{weights()}.
#'
#' \strong{Variance}: Analytical delta variance for EL is not implemented. Requesting
#' \code{variance_method = "delta"} is coerced to \code{"none"} with a warning. For standard
#' errors in both IID and survey settings, use \code{variance_method = "bootstrap"}.
#'
#' @references
#' Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical
#' Association, 97(457), 193-200.
#'
#' @section Progress Reporting:
#' When \code{variance_method = "bootstrap"}, progress reporting is available via the
#' \code{progressr} package. To enable it:
#'
#' \preformatted{
#' library(progressr)
#' library(future)
#'
#' # Enable progress reporting
#' handlers(global = TRUE)
#' handlers("txtprogressbar")  # or "progress", "cli", etc.
#'
#' # Set parallel backend (optional)
#' plan(multisession, workers = 4)
#'
#' # Always set seed for reproducibility
#' set.seed(123)
#'
#' # Run with progress bar
#' result <- nmar(Y ~ X, data = df,
#'                engine = el_engine(variance_method = "bootstrap",
#'                                   bootstrap_reps = 500))
#'
#' # Reset to sequential
#' plan(sequential)
#' }
#'
#' To disable progress in simulations or batch jobs:
#'
#' \code{handlers("void")  # Silent}
#'
#' If progressr is not installed or no handlers are set, bootstrap runs silently
#' (default behavior). Progress reporting works with all future backends and does
#' not affect reproducibility.
#'
#' @return An engine object of class \code{c("nmar_engine_el","nmar_engine")}.
#'   This is a configuration list; it is not a fit. Pass it to \link{nmar}.
#'
#' @seealso \link{nmar}; see the vignette "Empirical Likelihood Theory for NMAR" for derivations.
#'
#' @references
#' Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. \emph{Journal of the American
#' Statistical Association}, 97(457), 193-200.
#'
#' Wu, C., and Sitter, R. R. (2001). A model-calibration approach to using
#' complete auxiliary information from survey data. \emph{Journal of the American
#' Statistical Association}, 96(453), 185-193. Related to design-based
#' calibration; our EL approach balances auxiliary moments through empirical
#' likelihood constraints rather than calibration adjustments to weights.
#'
#'
#' @keywords engine
#' @export
#' @seealso [nmar()], [weights.nmar_result()], [summary.nmar_result]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200
#' X <- rnorm(n)
#' Z <- rnorm(n)
#' Y <- 2 + 0.5 * X + Z
#' p <- plogis(-0.7 + 0.4 * scale(Y)[, 1])
#' R <- runif(n) < p
#' df <- data.frame(Y_miss = Y, X = X)
#' df$Y_miss[!R] <- NA_real_
#' eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
#' fit <- nmar(Y_miss ~ X, data = df, engine = eng)
#' summary(fit)
#'
#' # Response-only predictors can be placed to the right of `|`:
#' df2 <- data.frame(Y_miss = Y, X = X, Z = Z)
#' df2$Y_miss[!R] <- NA_real_
#' eng2 <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
#' fit2 <- nmar(Y_miss ~ X | Z, data = df2, engine = eng2)
#' print(fit2)
#'
#' # Survey design usage
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   des <- survey::svydesign(ids = ~1, weights = ~1, data = df)
#'   eng3 <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
#'   fit3 <- nmar(Y_miss ~ X, data = des, engine = eng3)
#'   summary(fit3)
#' }
#' }
el_engine <- function(
    standardize = TRUE,
    trim_cap = Inf,
    on_failure = c("return", "error"),
    variance_method = c("delta", "bootstrap", "none"),
    bootstrap_reps = 500,
    auxiliary_means = NULL,
    control = list(),
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
    n_total = n_total,
    start = start,
    family = family
  )
# Validate and coerce unsupported variance modes upfront
  validate_nmar_engine_el(engine)
  if (identical(engine$variance_method, "delta")) {
    warning("Empirical likelihood 'delta' variance is not implemented; using variance_method='none'.", call. = FALSE)
    engine$variance_method <- "none"
  }
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
  validator$assert_list(engine, name = "engine")
  validator$assert_logical(engine$standardize, name = "standardize")

  validator$assert_positive_number(engine$trim_cap, name = "trim_cap", allow_infinite = TRUE)
  validator$assert_positive_integer(engine$bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  validator$assert_choice(engine$on_failure, choices = c("return", "error"), name = "on_failure")
  validator$assert_choice(engine$variance_method, choices = c("delta", "bootstrap", "none"), name = "variance_method")

  validator$assert_named_numeric(engine$auxiliary_means, name = "auxiliary_means", allow_null = TRUE)
  validator$assert_list(engine$control, name = "control")

  if (!is.null(engine$n_total)) validator$assert_positive_number(engine$n_total, name = "n_total", allow_infinite = FALSE)

  fam <- engine$family
  if (!is.list(fam) || is.null(fam$name) || !is.function(fam$linkinv) || !is.function(fam$mu.eta) || !is.function(fam$score_eta)) {
    stop("family must be 'logit'/'probit' or a valid family object with linkinv, mu.eta, score_eta")
  }

  engine
}
