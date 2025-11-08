#' @title Not Missing at Random (NMAR) Estimation
#'
#' @description Provides a unified interface for Not Missing at Random (NMAR) estimation.
#'   This function orchestrates the estimation process by validating inputs and
#'   dispatching to the appropriate engine based on the provided `engine` object.
#'   It ensures that all necessary data and model specifications are correctly
#'   formatted before computation begins.
#'
#' @param formula A two-sided formula specifying the model structure.
#'
#'   Use the partitioned form `y_miss ~ aux | response` where:
#'   - The left-hand side is the outcome (with `NA` indicating nonresponse);
#'   - The part left of `|` lists auxiliary variables (used in auxiliary moment
#'     constraints and in printing/reporting);
#'   - The part right of `|` lists response-model predictors (used only in the
#'     missingness model).
#'
#'   Both RHS partitions support standard R formula features — transforms
#'   (e.g., `I(Z^2)`), interactions (e.g., `X1:X2`, `X1*X2`), and parentheses.
#'   The original formula environment is preserved and evaluated via
#'   `model.matrix()`, so transforms behave exactly as in base modeling.
#'
#'   If `|` is omitted, the RHS is interpreted as auxiliaries only and the
#'   response model uses only the outcome (no additional response-only
#'   predictors). The outcome variable is always included implicitly in the
#'   response model.
#'
#'   Engine-specific rules (e.g., whether the outcome may appear on the response
#'   side, or whether auxiliary and response sets may overlap) are enforced via
#'   [`engine_traits()`]. For example, the empirical likelihood engine permits
#'   overlap and the outcome on the RHS; the exponential-tilting engine does not.
#' @param data A data frame or `survey.design` containing the variables the
#'   formula references. For survey designs, the formula is evaluated against
#'   `design$variables`.
#' @param engine An engine configuration object created by a constructor such as
#'   [el_engine()] or [exptilt_engine()]. This object defines the specific NMAR
#'   estimation method and its parameters. It must inherit from class
#'   `nmar_engine`.
#' @param trace_level Integer 0-3; controls verbosity level during estimation (default: 0):
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
#'
#' @details
#'   \strong{Auxiliary means naming}: For engines that use auxiliary moment
#'   constraints (e.g., [el_engine()]), when supplying `auxiliary_means`, the
#'   names must match the column names of the auxiliary model matrix generated
#'   from the left-of-`|` RHS after removing the intercept. If transforms or
#'   interactions are used, names should follow `model.matrix()` conventions —
#'   for example, `I(X^2)` for squared terms, and `X1:X2` for interaction terms.
#'   Mismatched or missing names are dropped, and if no names match, auxiliary
#'   constraints are disabled.
#'
#'   \strong{Factor variables}: Categorical variables (factors or character
#'   vectors) are automatically converted to dummy variables via `model.matrix()`.
#'   The first level (alphabetically) is used as the reference category, and
#'   coefficients represent differences from the reference. For example, if
#'   `region` is a factor with levels `c("North", "South", "West")`, the formula
#'   `Y ~ region` creates dummy variables `regionSouth` and `regionWest`, with
#'   coefficients representing effects relative to North (the reference).
#'
#'   \strong{Formula transformations}: All standard R formula transformations work
#'   correctly. Transformations such as `I(X^2)`, `log(X)`, `poly(X, 2)` on the
#'   right-hand side are evaluated first, then the resulting transformed variables
#'   are standardized (if `standardize = TRUE` in the engine). Coefficients are
#'   automatically unscaled to the original parameter space after estimation.
#'
#'   \strong{Outcome transformations (IMPORTANT)}: Transformations applied to the
#'   left-hand side outcome variable estimate the mean of the \emph{transformed}
#'   outcome, NOT the mean of the original outcome. Formally, `nmar(g(Y) ~ X, ...)`
#'   estimates E[g(Y)], not g(E[Y]).
#'
#'   \itemize{
#'     \item `nmar(log(Y) ~ X, ...)` estimates E[log(Y)], the expected log-outcome
#'     \item `nmar(sqrt(Y) ~ X, ...)` estimates E[sqrt(Y)], the expected square root
#'     \item `nmar(I(Y^2) ~ X, ...)` estimates E[Y^2], the second moment
#'   }
#'
#'   To obtain estimates on the original scale, you must back-transform the point
#'   estimate. However, note that for nonlinear transformations, the naive
#'   back-transform is \emph{biased} due to Jensen's inequality. For example, if
#'   you estimate E[log(Y)] = 2.5 with SE = 0.3, then exp(2.5) does NOT equal
#'   E[Y] (it underestimates the mean due to convexity). Bias-corrected
#'   back-transformations for common cases (log-normal, Box-Cox) are beyond the
#'   scope of this function; consult the literature or compute bootstrap estimates
#'   on the original scale directly (by specifying `Y ~ X` rather than `log(Y) ~ X`).
#'
#'   \strong{Environment preservation}: The formula's original environment is
#'   preserved throughout the estimation process. Any user-defined functions
#'   used in transforms must be visible in that environment.
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
  validate_nmar_args(spec, traits)

# Wrap the validated spec and engine traits into a task object so every
# engine sees the same downstream interface.
  task <- new_nmar_task(spec, traits, trace_level = trace_level)

  result <- run_engine(engine, task)
  if (!inherits(result, "nmar_result")) {
    stop("`run_engine()` must return an object inheriting from 'nmar_result'.", call. = FALSE)
  }
  primary_class <- setdiff(class(result), "nmar_result")
  target_class <- if (length(primary_class)) primary_class[[1]] else "nmar_result"
  result <- validate_nmar_result(result, target_class)
  result
}

run_engine <- function(engine, task) {
  UseMethod("run_engine")
}
