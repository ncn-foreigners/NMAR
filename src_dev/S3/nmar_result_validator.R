#' Validate nmar_result structure
#'
#' Ensures both the child class and the parent schema are satisfied. The
#' validator also back-fills defaults so downstream code can rely on the
#' presence of optional components without defensive checks.
#'
#' @keywords internal
validate_nmar_result <- function(x, class_name) {
  stopifnot(is.list(x), inherits(x, class_name), inherits(x, "nmar_result"))

# Core scalars
  if (is.null(x$estimate)) stop("`estimate` must be supplied by the result object.")
  if (length(x$estimate) != 1 || !is.numeric(x$estimate)) stop("`estimate` must be a numeric scalar.")

  if (is.null(x$se)) x$se <- NA_real_
  if (length(x$se) != 1 || !is.numeric(x$se)) stop("`se` must be a numeric scalar (NA allowed).")

  if (is.null(x$estimate_name)) x$estimate_name <- NA_character_
  if (length(x$estimate_name) != 1) stop("`estimate_name` must be a character scalar (NA allowed).")

  if (is.null(x$converged)) x$converged <- NA
  if (length(x$converged) != 1 || !is.logical(x$converged)) stop("`converged` must be a logical scalar.")

# Model
  if (is.null(x$model) || !is.list(x$model)) x$model <- list()
  if (is.null(x$model$coefficients)) x$model$coefficients <- NULL
  if (is.null(x$model$vcov)) x$model$vcov <- NULL

# Weights
  if (is.null(x$weights_info) || !is.list(x$weights_info)) x$weights_info <- list()
  if (is.null(x$weights_info$values)) x$weights_info$values <- NULL
  if (is.null(x$weights_info$trimmed_fraction)) x$weights_info$trimmed_fraction <- NA_real_

# Sample metadata
  sample_defaults <- list(n_total = NA_integer_, n_respondents = NA_integer_, is_survey = FALSE, design = NULL)
  if (is.null(x$sample) || !is.list(x$sample)) x$sample <- list()
  for (nm in names(sample_defaults)) {
    if (is.null(x$sample[[nm]])) x$sample[[nm]] <- sample_defaults[[nm]]
  }

# Inference metadata
  inference_defaults <- list(
    variance_method = NA_character_,
    df = NA_real_,
    message = NA_character_,
    used_pseudoinverse = FALSE,
    used_ridge = FALSE
  )
  if (is.null(x$inference) || !is.list(x$inference)) x$inference <- list()
  for (nm in names(inference_defaults)) {
    if (is.null(x$inference[[nm]])) x$inference[[nm]] <- inference_defaults[[nm]]
  }

# Diagnostics / metadata
  if (is.null(x$diagnostics) || !is.list(x$diagnostics)) x$diagnostics <- list()
  if (is.null(x$meta) || !is.list(x$meta)) x$meta <- list()
  meta_defaults <- list(engine_name = NA_character_, call = NULL, formula = NULL)
  for (nm in names(meta_defaults)) {
    if (is.null(x$meta[[nm]])) x$meta[[nm]] <- meta_defaults[[nm]]
  }

  if (is.null(x$extra) || !is.list(x$extra)) x$extra <- list()

  x
}
