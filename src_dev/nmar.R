#' @title Not Missing at Random (NMAR) Estimation
#'
#' @description High-level interface for NMAR estimation.
#'
#' \code{nmar()} validates basic inputs and dispatches to an engine (for example
#' \code{\link{el_engine}}). The engine controls the estimation method and
#' interprets \code{formula}; see the engine documentation for model-specific
#' requirements.
#'
#' @param formula A two-sided formula. Many engines support a partitioned
#'   right-hand side via \code{|}, for example
#'   \code{y_miss ~ block1_vars | block2_vars}. The meaning of these blocks is
#'   engine-specific (see the engine documentation). In the common
#'   "missing values indicate nonresponse" workflow, the left-hand side is the
#'   outcome with \code{NA} values for nonrespondents.
#' @param data A \code{data.frame} or a \code{survey.design} containing the
#'   variables referenced by \code{formula}.
#' @param engine An NMAR engine configuration object, typically created by
#'   \code{\link{el_engine}}, \code{\link{exptilt_engine}}, or
#'   \code{\link{exptilt_nonparam_engine}}. This object defines the estimation
#'   method and tuning parameters and must inherit from class
#'   \code{"nmar_engine"}.
#' @param trace_level Integer 0-3; controls verbosity during estimation
#'   (default \code{0}):
#'   \itemize{
#'     \item 0: no output (silent mode);
#'     \item 1: major steps only (initialization, convergence, final results);
#'     \item 2: iteration summaries and key diagnostics;
#'     \item 3: full diagnostic output.
#'   }
#'
#' @return An object of class \code{"nmar_result"} with an engine-specific subclass
#'   (for example \code{"nmar_result_el"}). Use \code{summary()},
#'   \code{\link{se}}, \code{confint()}, \code{weights()}, \code{coef()},
#'   \code{fitted()}, and \code{generics::tidy()} / \code{generics::glance()} to
#'   access estimates, standard errors, weights, and diagnostics.
#'
#' @seealso \code{\link{el_engine}}, \code{\link{exptilt_engine}},
#'   \code{\link{exptilt_nonparam_engine}}, \code{\link{summary.nmar_result}},
#'   \code{\link{weights.nmar_result}}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n)
#' z1 <- rnorm(n)
#' y_true <- 0.5 + 0.3 * x1 + 0.2 * z1 + rnorm(n, sd = 0.3)
#' resp <- rbinom(n, 1, plogis(2 + 0.1 * y_true + 0.1 * z1))
#' if (all(resp == 1)) resp[sample.int(n, 1)] <- 0L
#' y_obs <- ifelse(resp == 1, y_true, NA_real_)
#'
#' # Empirical likelihood engine
#' df_el <- data.frame(Y_miss = y_obs, X = x1, Z = z1)
#' eng_el <- el_engine(variance_method = "none")
#' fit_el <- nmar(Y_miss ~ X | Z, data = df_el, engine = eng_el)
#' summary(fit_el)
#'
#' \donttest{
#' # Exponential tilting engine (illustrative)
#' dat_et <- data.frame(y = y_obs, x2 = z1, x1 = x1)
#' eng_et <- exptilt_engine(
#'   y_dens = "normal",
#'   family = "logit",
#'   variance_method = "none"
#' )
#' fit_et <- nmar(y ~ x2 | x1, data = dat_et, engine = eng_et)
#' summary(fit_et)
#'
#' # Survey design example (same outcome, random weights)
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   w <- runif(n, 0.5, 2)
#'   des <- survey::svydesign(ids = ~1, weights = ~w,
#'                            data = data.frame(Y_miss = y_obs, X = x1, Z = z1))
#'   eng_svy <- el_engine(variance_method = "none")
#'   fit_svy <- nmar(Y_miss ~ X | Z, data = des, engine = eng_svy)
#'   summary(fit_svy)
#' }
#'
#' # Bootstrap variance usage
#' if (requireNamespace("future.apply", quietly = TRUE)) {
#'   set.seed(2)
#'   eng_boot <- el_engine(
#'     variance_method = "bootstrap",
#'     bootstrap_reps = 20
#'   )
#'   fit_boot <- nmar(Y_miss ~ X | Z, data = df_el, engine = eng_boot)
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
