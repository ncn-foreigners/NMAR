#' @title Not Missing at Random (NMAR) Estimation
#'
#' @description High-level interface for estimating population means and related
#'   statistics under Not Missing At Random (NMAR) nonresponse. Uses a unified
#'   formula interface and an engine object to fit NMAR estimators on data
#'   frames or survey designs.
#'
#' @param formula A two-sided formula of the form `y_miss ~ aux1 + aux2 | z1 + z2`.
#'   The left-hand side is the outcome (with `NA` values indicating nonresponse).
#'   The right-hand side is split by `|` into two parts:
#'   - left of `|`: auxiliary variables that enter moment constraints;
#'   - right of `|`: missingness (response) model predictors that enter the
#'     missingness model only.
#'   If `|` is omitted, only auxiliary variables are used. The outcome variable
#'   is implicitly included in the missingness model via the evaluated left-hand
#'   side.
#' @param data A `data.frame` or a `survey.design` containing the variables
#'   referenced by `formula`.
#' @param engine An NMAR engine configuration object, typically created by
#'   [el_engine()], [exptilt_engine()], or [exptilt_nonparam_engine()]. This
#'   object defines the estimation method and its tuning parameters and must
#'   inherit from class `"nmar_engine"`.
#' @param trace_level Integer 0–3; controls verbosity during estimation
#'   (default 0):
#'   \itemize{
#'     \item 0: no output (silent mode);
#'     \item 1: major steps only (initialization, convergence, final results);
#'     \item 2: iteration summaries and key diagnostics;
#'     \item 3: full diagnostic output.
#'   }
#'
#' @return An object of class `"nmar_result"` with an engine-specific subclass
#'   (for example `"nmar_result_el"`). Methods such as [summary()],
#'   [se()], [confint()], [generics::tidy()], [generics::glance()],
#'   [weights()], [coef()], and [fitted()] provide access to estimates,
#'   standard errors, weights, and diagnostics.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n)
#' y_true <- 0.5 + 0.3 * x1 + rnorm(n, sd = 0.3)
#' resp <- rbinom(n, 1, plogis(2 + 0.1 * y_true))
#' if (all(resp == 1)) resp[sample.int(n, 1)] <- 0L
#' y_obs <- ifelse(resp == 1, y_true, NA_real_)
#'
#' # Empirical likelihood engine
#' df_el <- data.frame(Y_miss = y_obs, X = x1)
#' eng_el <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
#' fit_el <- nmar(Y_miss ~ X, data = df_el, engine = eng_el)
#' summary(fit_el)
#'
#' \donttest{
#' # Exponential tilting engine on the same data
#' dat_et <- data.frame(y = y_obs, x1 = x1)
#' eng_et <- exptilt_engine(
#'   y_dens = "normal",
#'   family = "logit",
#'   variance_method = "none",
#' )
#' fit_et <- nmar(y ~ x1, data = dat_et, engine = eng_et)
#' summary(fit_et)
#'
#' # Survey design example (same outcome, random weights)
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   w <- runif(n, 0.5, 2)
#'   des <- survey::svydesign(ids = ~1, weights = ~w,
#'                            data = data.frame(y = y_obs, x1 = x1))
#'   eng_svy <- el_engine(auxiliary_means = NULL, variance_method = "none")
#'   fit_svy <- nmar(y ~ x1, data = des, engine = eng_svy)
#'   summary(fit_svy)
#' }
#'
#' # Bootstrap variance usage
#' if (requireNamespace("future.apply", quietly = TRUE)) {
#'   eng_boot <- el_engine(
#'     auxiliary_means = c(X = 0),
#'     variance_method = "bootstrap",
#'     bootstrap_reps = 50
#'   )
#'   fit_boot <- nmar(Y_miss ~ X, data = df_el, engine = eng_boot)
#'   se(fit_boot)
#' }
#' }
#' @keywords nmar
#' @export
nmar <- function(formula, data, engine, trace_level = 0) {
  stopifnot(inherits(engine, "nmar_engine"))

  validator_assert_choice(trace_level, choices = 0:3, name = "trace_level")
  validate_data(data)

# Dispatch to engine
  res <- run_engine(engine, formula, data, trace_level)

# Ensure meta$call records the outer nmar() call with the actual formula object
  nmar_call <- match.call()
  nmar_call$formula <- formula
  if (is.null(res$meta) || !is.list(res$meta)) res$meta <- list()
  res$meta$call <- nmar_call
  res
}

run_engine <- function(engine, formula, data, trace_level) {
  UseMethod("run_engine")
}
