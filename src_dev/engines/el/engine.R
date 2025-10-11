#' Empirical likelihood (EL) engine for NMAR
#'
#' Constructs a configuration for the empirical likelihood estimator under
#' nonignorable nonresponse (NMAR), following Qin, Leung and Shao (2002).
#' Pass the returned object to `nmar()` together with a formula and data.
#'
#' @param family Response model family; a string ("logit" or "probit") or a
#'   family object created by `logit_family()` / `probit_family()`.
#' @details
#' The EL estimator uses weights that satisfy estimating equations for the
#' response mechanism and optional auxiliary moment constraints. The response
#' model score is the derivative of the Bernoulli log‑likelihood with respect
#' to the linear predictor, i.e. `mu.eta(eta)/linkinv(eta)`, valid for both
#' logit and probit links. Response predictors do not need to coincide with
#' auxiliary predictors; only auxiliaries require known population moments.
#' See the theory note and Qin, Leung and Shao (2002) for derivations.
#'
#' @param standardize Logical; standardize predictors. Default TRUE.
#' @param trim_cap Numeric; cap for EL weights (Inf = no trimming).
#' @param on_failure Character; "return" or "error" on solver failure.
#' @param variance_method Character; one of "delta", "bootstrap", or "none".
#' @param variance_jacobian Character; "auto", "analytic", or "numeric".
#' @param solver_jacobian Character; "auto", "analytic", or "none".
#' @param solver_method Character; policy for the root solver:
#'   - "auto" (default): Newton with analytic Jacobian, perturbation restarts, then Broyden fallback.
#'   - "newton": Newton only (no Broyden fallback).
#'   - "broyden": Broyden only (secant updates, no analytic Jacobian).
#' @param variance_pseudoinverse Logical; allow pseudo-inverse for variance.
#' @param variance_ridge Logical or numeric; if TRUE, apply adaptive ridge to
#'   stabilize Jacobian inversion; if numeric, treated as ridge epsilon.
#' @param bootstrap_reps Integer; reps when variance_method = "bootstrap".
#' @param suppress_warnings Logical; suppress variance method warnings.
#' @param auxiliary_means Named numeric vector; population means for auxiliaries
#'   (names must match RHS of outcome formula). Optional.
#' @param control List; optional solver control for `nleqslv::nleqslv(control=...)`.
#'   Recognized fields include `xtol`, `ftol`, `btol`, `maxit`, `trace`, `stepmax`,
#'   and `delta`.
#' @param solver_args List; optional top‑level `nleqslv` arguments such as
#'   `global` (e.g., "dbldog", "qline") and `xscalm` (e.g., "auto", "fixed").
#'   These are passed at the top level of the `nleqslv` call and take precedence
#'   over similarly named entries provided in `control`.
#' @param n_total Optional when supplying respondents-only data (no NA in the
#'   outcome). For `data.frame` inputs, set to the total number of sampled
#'   units before filtering to respondents. For `survey.design` inputs, set to
#'   the total design weight or known population total. If omitted and the
#'   outcome contains no NAs, the estimator errors with an instruction to
#'   provide `n_total`.
#'
#' @return An engine object of class `c('nmar_engine_el','nmar_engine')`.
#'   This is a configuration list; it is not a fit. Pass it to `nmar()`.
#'
#' @seealso `nmar()`; see the vignette "Empirical Likelihood Theory for NMAR" for derivations.
#' @references
#' Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American
#' Statistical Association, 97(457), 193–200.
#' @keywords engine
#' @export
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
#' eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "delta")
#' fit <- nmar(Y_miss ~ X, data = df, engine = eng)
#' summary(fit)
#' }
el_engine <- function(
    standardize = TRUE,
    trim_cap = Inf,
    on_failure = c("return", "error"),
    variance_method = c("delta", "bootstrap", "none"),
    variance_jacobian = c("auto", "analytic", "numeric"),
    solver_jacobian = c("auto", "analytic", "none"),
    solver_method = c("auto", "newton", "broyden"),
    solver_args = list(),
    variance_pseudoinverse = FALSE,
    variance_ridge = FALSE,
    bootstrap_reps = 500,
    suppress_warnings = FALSE,
    auxiliary_means = NULL,
    control = list(),
    n_total = NULL,
    family = c("logit", "probit")) {
  on_failure <- match.arg(on_failure)
  if (is.null(variance_method)) variance_method <- "none"
  variance_method <- match.arg(variance_method)
  variance_jacobian <- match.arg(variance_jacobian)
  solver_jacobian <- match.arg(solver_jacobian)
  solver_method <- match.arg(solver_method)
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
    variance_jacobian = variance_jacobian,
    solver_jacobian = solver_jacobian,
    solver_method = solver_method,
    solver_args = solver_args,
    variance_pseudoinverse = variance_pseudoinverse,
    variance_ridge = variance_ridge,
    bootstrap_reps = bootstrap_reps,
    suppress_warnings = suppress_warnings,
    auxiliary_means = auxiliary_means,
    control = control,
    n_total = n_total,
    family = family
  )
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
  stopifnot(is.list(engine))
# Logical flags
  if (!is.logical(engine$standardize) || length(engine$standardize) != 1) stop("standardize must be a single logical")
  if (!is.logical(engine$variance_pseudoinverse) || length(engine$variance_pseudoinverse) != 1) stop("variance_pseudoinverse must be a single logical")
  if (!is.logical(engine$suppress_warnings) || length(engine$suppress_warnings) != 1) stop("suppress_warnings must be a single logical")
# Numeric settings
  if (!is.numeric(engine$trim_cap) || length(engine$trim_cap) != 1 || (!is.infinite(engine$trim_cap) && engine$trim_cap <= 0)) stop("trim_cap must be a positive number or Inf")
  if (!is.numeric(engine$bootstrap_reps) || length(engine$bootstrap_reps) != 1 || engine$bootstrap_reps < 1) stop("bootstrap_reps must be a positive integer")
# Factors handled via match.arg upstream; here we assert presence
  if (!engine$on_failure %in% c("return", "error")) stop("on_failure must be 'return' or 'error'")
  if (!engine$variance_method %in% c("delta", "bootstrap", "none")) stop("variance_method must be 'delta', 'bootstrap', or 'none'")
  if (!engine$variance_jacobian %in% c("auto", "analytic", "numeric")) stop("variance_jacobian must be one of 'auto','analytic','numeric'")
  if (!engine$solver_jacobian %in% c("auto", "analytic", "none")) stop("solver_jacobian must be one of 'auto','analytic','none'")
# Auxiliary means: NULL or named numeric
  if (!is.null(engine$auxiliary_means)) {
    if (!is.numeric(engine$auxiliary_means) || is.null(names(engine$auxiliary_means)) || anyNA(names(engine$auxiliary_means))) {
      stop("auxiliary_means must be a named numeric vector or NULL")
    }
  }
# Control should be a list (passed through)
  if (!is.list(engine$control)) stop("control must be a list")
  if (!is.list(engine$solver_args)) stop("solver_args must be a list of top-level nleqslv args (e.g., global, xscalm)")
# Optional n_total must be a positive integer-like scalar if provided
  if (!is.null(engine$n_total)) {
    if (!is.numeric(engine$n_total) || length(engine$n_total) != 1 || !is.finite(engine$n_total) || engine$n_total <= 0) {
      stop("n_total must be a positive number when provided")
    }
  }
# Shallow validation of known top-level args
  if (!is.null(engine$solver_args$global) && !(is.character(engine$solver_args$global) && length(engine$solver_args$global) == 1)) {
    stop("solver_args$global must be a single character value (e.g., 'dbldog', 'qline').")
  }
  if (!is.null(engine$solver_args$xscalm) && !(is.character(engine$solver_args$xscalm) && length(engine$solver_args$xscalm) == 1)) {
    stop("solver_args$xscalm must be a single character value (e.g., 'auto', 'fixed').")
  }
  fam <- engine$family
  if (!is.list(fam) || is.null(fam$name) || !is.function(fam$linkinv) || !is.function(fam$mu.eta) || !is.function(fam$score_eta)) {
    stop("family must be 'logit'/'probit' or a valid family object with linkinv, mu.eta, score_eta")
  }
# variance_ridge can be logical or positive numeric epsilon
  if (!(is.logical(engine$variance_ridge) || (is.numeric(engine$variance_ridge) && is.finite(engine$variance_ridge) && engine$variance_ridge > 0))) {
    stop("variance_ridge must be logical or a positive numeric epsilon")
  }
  engine
}
