#' Engine trait declarations and defaults
#'
#' Canonical defaults and S3 generic for engine traits. Traits describe how an
#' engine expects input validation to behave. Users can call `engine_traits()`
#' on an engine object to introspect behaviour such as whether the outcome may
#' appear on the response RHS explicitly or whether multiple outcomes are
#' supported.
#'
#' Important: The NMAR response model always includes the outcome implicitly by
#' design. The `allow_outcome_in_missingness` trait governs whether the outcome
#' may be specified explicitly on the response RHS (e.g., `Y ~ aux | Y + Z`).
#'
#' @keywords internal
NMAR_DEFAULT_TRAITS <- list(
  allow_outcome_in_missingness = FALSE,
  allow_covariate_overlap = FALSE,
  requires_single_outcome = TRUE,
  allow_respondents_only = FALSE,
  drop_auxiliary_intercept = TRUE
)

#' Engine trait S3 generic
#'
#' @param engine An object inheriting from class `nmar_engine`.
#' @return A named list of trait flags.
#' @keywords engine_view
#' @export
engine_traits <- function(engine) {
  UseMethod("engine_traits")
}

#' @keywords engine_view
#' @export
engine_traits.default <- function(engine) {
  NMAR_DEFAULT_TRAITS
}

#' @keywords engine_view
#' @export
engine_traits.nmar_engine <- function(engine) {
  engine_traits.default(engine)
}
