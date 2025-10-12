#' @title Not Missing at Random (NMAR) Estimation
#'
#' @description Provides a unified interface for Not Missing at Random (NMAR) estimation.
#'   This function orchestrates the estimation process by validating inputs and
#'   dispatching to the appropriate engine based on the provided `engine` object.
#'   It ensures that all necessary data and model specifications are correctly
#'   formatted before computation begins.
#'
#' @param formula A two-sided formula of the form `y_miss ~ x1 + x2 + ...` specifying the
#'   outcome (with `NA` values indicating nonresponse) and auxiliary variables used to
#'   stabilise estimation.
#' @param data A data frame or `survey.design` containing the variables referenced by the
#'   formula.
#' @param engine An engine configuration object, typically created by an
#'   engine constructor function like `exptilt()`. This object defines the
#'   specific NMAR estimation method and its parameters. It must inherit from
#'   class `nmar_engine`.
#' @param response_predictors Optional character vector naming additional predictors for the
#'   response (missingness) model. These variables need not appear on the right-hand side of
#'   `formula` and are interpreted by the chosen engine.
#'
#' @return An object containing the estimation results, whose structure will be
#'   specific to the `engine` used. This might include estimated parameters,
#'   convergence information, and other relevant output from the chosen NMAR method.
#' @keywords nmar
#' @export
nmar <- function(formula, data, engine, response_predictors = NULL) {
  stopifnot(inherits(engine, "nmar_engine"))

  spec <- parse_nmar_spec(
    formula = formula,
    data = data,
    response_predictors = response_predictors,
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

  run_engine(engine, task)
}

run_engine <- function(engine, task) {
  UseMethod("run_engine")
}
