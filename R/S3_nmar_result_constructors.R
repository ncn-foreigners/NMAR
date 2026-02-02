#' Construct for result objects
#'
#' Builds an `nmar_result` list using the shared schema and validates it.
#'
#' @details
#' Engine-level constructors should call this helper with named arguments rather
#' than assembling result lists by hand. At minimum, engines should supply
#' \code{estimate} (numeric scalar) and \code{converged} (logical). All other
#' fields are optional:
#' \itemize{
#' \item \code{estimate_name}: label for the primary estimand (defaults to
#' \code{NA_character_} if omitted).
#' \item \code{se}: standard error for the primary estimand (defaults to
#' \code{NA_real_} when not available).
#' \item \code{model}, \code{weights_info}, \code{sample}, \code{inference},
#' \code{diagnostics}, \code{meta}, \code{extra}: lists that may be partially
#' pecified or \code{NULL}; \code{validate_nmar_result()} will back-fill
#' missing subfields with safe defaults.
#' item \code{class}: engine-specific result subclass name, e.g.
#' \code{"nmar_result_el"}; it is combined with the parent class
#' \code{"nmar_result"}.
#' }
#'
#' Calling \code{new_nmar_result()} ensures that every engine returns objects
#' that satisfy the shared schema and are immediately compatible with parent
#' S3 methods such as \code{vcov()}, \code{confint()}, \code{tidy()},
#' \code{glance()}, and \code{weights()}.
#'
#' @keywords internal
new_nmar_result <- function(...) {
  dots <- list(...)

  class_name <- dots$class %||% "nmar_result"

  result <- list(
    y_hat = dots$estimate,
    estimate_name = dots$estimate_name,
    se = dots$se,
    converged = dots$converged,
    model = dots$model,
    weights_info = dots$weights_info,
    sample = dots$sample,
    inference = dots$inference,
    diagnostics = dots$diagnostics,
    meta = dots$meta,
    extra = dots$extra
  )

  class(result) <- c(class_name, "nmar_result")
  validate_nmar_result(result, class_name)
}

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
