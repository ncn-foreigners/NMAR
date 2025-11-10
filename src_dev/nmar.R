#' @title Not Missing at Random (NMAR) Estimation
#'
#' @description Provides a unified interface for Not Missing at Random (NMAR) estimation.
#'   This function orchestrates the estimation process by validating inputs and
#'   dispatching to the appropriate engine based on the provided `engine` object.
#'   It ensures that all necessary data and model specifications are correctly
#'   formatted before computation begins.
#'
#' @param formula A two-sided formula of the form `y_miss ~ aux1 + aux2 | z1 + z2`.
#'   The left-hand side is the outcome (with `NA` values indicating nonresponse).
#'   The right-hand side is split by `|` into two parts:
#'   - left of `|`: auxiliary variables (enter moment constraints);
#'   - right of `|`: response-model predictors (enter the missingness model only).
#'   If `|` is omitted, only auxiliary variables are used for both parsing and printing.
#'   The outcome variable is implicitly included in the response model.
#' @param data A data frame or `survey.design` containing the variables referenced by the
#'   formula.
#' @param engine An engine configuration object, typically created by an
#'   engine constructor function like `exptilt()`. This object defines the
#'   specific NMAR estimation method and its parameters. It must inherit from
#'   class `nmar_engine`.
#' @param trace_level Integer 0-3; controls verbosity level during estimation (default: 1):
#'   \itemize{
#'     \item 0: No output (silent mode)
#'     \item 1: Major steps only (initialization, convergence, final results)
#'     \item 2: Moderate detail (iteration summaries, key diagnostics)
#'     \item 3: Full detail (all diagnostics, intermediate values)
#'   }
#'
#' @return An object containing the estimation results, whose structure will be
#'   specific to the `engine` used. This might include estimated parameters,
#'   convergence information, and other relevant output from the chosen NMAR method.
#' @keywords nmar
#' @export
nmar <- function(formula, data, engine, trace_level = 0) {
  stopifnot(inherits(engine, "nmar_engine"))

  validator$assert_choice(trace_level, choices = 0:3, name = "trace_level")

  spec <- parse_nmar_spec(
    formula = formula,
    data = data,
    env = parent.frame()
  )

  traits <- engine_traits(engine)
# Activate respondents-only relaxation only when the engine supplies
# the required extra information (currently: n_total for EL).
  if (isTRUE(traits$allow_respondents_only)) {
    has_n_total <- !is.null(engine$n_total)
    traits$allow_respondents_only <- isTRUE(has_n_total)
  }
  validate_nmar_args(spec, traits)

# Wrap the validated spec and engine traits into a task object so every
# engine sees the same downstream interface.
  task <- new_nmar_task(spec, traits)

# Pass trace_level to the engine
  task$trace_level <- trace_level
# Preserve the outer nmar() call for result metadata
  task$nmar_call <- match.call()

  run_engine(engine, task)
}

run_engine <- function(engine, task) {
  UseMethod("run_engine")
}
