#' Validate nmar_result
#'
#' Ensures both the child class and the parent schema are satisfied. The
#' validator also back-fills defaults so downstream code can rely on the
#' presence of optional components without defensive checks.
#'
#' @details
#' This helper is the single authority on the `nmar_result` schema. It expects
#' a list that already carries class \code{c(class_name, "nmar_result")} and
#' at least a primary estimate stored in \code{y_hat}. All other components are
#' optional. When they are \code{NULL} or missing, the validator supplies safe
#' defaults:
#' \itemize{
#' \item Core scalars: \code{se} (numeric, default \code{NA_real_}),
#' \code{estimate_name} (character, default \code{NA_character_}),
#' \code{converged} (logical, default \code{NA}).
#' \item \code{model}: list with \code{coefficients} and \code{vcov}, both
#' defaulting to \code{NULL}.
#' \item \code{weights_info}: list with \code{values} (default \code{NULL}) and
#' \code{trimmed_fraction} (default \code{NA_real_}).
#' \item \code{sample}: list with \code{n_total}, \code{n_respondents},
#' \code{is_survey}, and \code{design}, defaulted to missing/empty values.
#' \item \code{inference}: list with \code{variance_method}, \code{df}, and
#' \code{message}, all defaulted to missing values.
#' \item \code{diagnostics}, \code{meta}, and \code{extra}: defaulted to empty
#' lists, with \code{meta} carrying \code{engine_name}, \code{call}, and
#' \code{formula} when unset.
#' }
#'
#' Engine constructors should normally call \code{new_nmar_result()} rather than
#' invoking this function directly. \code{new_nmar_result()} attaches classes and
#' funnels all objects through \code{validate_nmar_result()} so downstream S3
#' methods can assume a consistent structure.
#'
#' @keywords internal
validate_nmar_result <- function(x, class_name) {
  stopifnot(is.list(x), inherits(x, class_name), inherits(x, "nmar_result"))

  if (is.null(x$y_hat)) stop("`y_hat` must be supplied by the result object.")
  validator_assert_scalar_numeric(x$y_hat, name = "y_hat", allow_na = TRUE, finite = FALSE)

  if (is.null(x$se)) x$se <- NA_real_
  validator_assert_scalar_numeric(x$se, name = "se", allow_na = TRUE, finite = FALSE)

  if (is.null(x$estimate_name)) x$estimate_name <- NA_character_
  validator_assert_scalar_character(x$estimate_name, name = "estimate_name", allow_na = TRUE, non_empty = FALSE)

  if (is.null(x$converged)) x$converged <- NA
  validator_assert_scalar_logical(x$converged, name = "converged")

  if (is.null(x$model) || !is.list(x$model)) x$model <- list()
  if (is.null(x$model$coefficients)) x$model$coefficients <- NULL
  if (is.null(x$model$vcov)) x$model$vcov <- NULL

  if (is.null(x$weights_info) || !is.list(x$weights_info)) x$weights_info <- list()
  if (is.null(x$weights_info$values)) x$weights_info$values <- NULL
  if (is.null(x$weights_info$trimmed_fraction)) x$weights_info$trimmed_fraction <- NA_real_

  sample_defaults <- list(n_total = NA_integer_, n_respondents = NA_integer_, is_survey = FALSE, design = NULL)
  if (is.null(x$sample) || !is.list(x$sample)) x$sample <- list()
  for (nm in names(sample_defaults)) {
    if (is.null(x$sample[[nm]])) x$sample[[nm]] <- sample_defaults[[nm]]
  }

  inference_defaults <- list(
    variance_method = NA_character_,
    df = NA_real_,
    message = NA_character_
  )
  if (is.null(x$inference) || !is.list(x$inference)) x$inference <- list()
  for (nm in names(inference_defaults)) {
    if (is.null(x$inference[[nm]])) x$inference[[nm]] <- inference_defaults[[nm]]
  }

  if (is.null(x$diagnostics) || !is.list(x$diagnostics)) x$diagnostics <- list()
  if (is.null(x$meta) || !is.list(x$meta)) x$meta <- list()
  meta_defaults <- list(engine_name = NA_character_, call = NULL, formula = NULL)
  for (nm in names(meta_defaults)) {
    if (is.null(x$meta[[nm]])) x$meta[[nm]] <- meta_defaults[[nm]]
  }

  if (is.null(x$extra) || !is.list(x$extra)) x$extra <- list()

  x
}
