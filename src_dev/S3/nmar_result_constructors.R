#' Construct Result Object (parent helper)
#'
#' Builds a normalized `nmar_result` list using the modern schema.
#' Engines must pass named fields; no legacy positional signature is supported.
#'
#' @keywords internal
new_nmar_result <- function(...) {
  dots <- list(...)

  y_hat <- dots$estimate
  se <- dots$se
  estimate_name <- dots$estimate_name %||% NA_character_
  converged <- dots$converged
  model <- dots$model %||% list()
  weights_info <- dots$weights_info %||% list()
  sample <- dots$sample %||% list()
  inference <- dots$inference %||% list()
  diagnostics <- dots$diagnostics %||% list()
  meta <- dots$meta %||% list()
  extra <- dots$extra %||% list()
  class_name <- dots$class %||% "nmar_result"

# Normalize components
  if (!is.list(model)) model <- list()
  if (is.null(model$coefficients)) model$coefficients <- NULL
  if (is.null(model$vcov)) model$vcov <- NULL

  if (!is.list(weights_info)) weights_info <- list()
  if (is.null(weights_info$values)) weights_info$values <- NULL
  if (is.null(weights_info$trimmed_fraction)) weights_info$trimmed_fraction <- NA_real_

  sample_defaults <- list(n_total = NA_integer_, n_respondents = NA_integer_, is_survey = FALSE, design = NULL)
  if (!is.list(sample)) sample <- list()
  for (nm in names(sample_defaults)) if (is.null(sample[[nm]])) sample[[nm]] <- sample_defaults[[nm]]

  inference_defaults <- list(
    variance_method = NA_character_,
    df = NA_real_,
    message = NA_character_
  )
  if (!is.list(inference)) inference <- list()
  for (nm in names(inference_defaults)) if (is.null(inference[[nm]])) inference[[nm]] <- inference_defaults[[nm]]

  meta_defaults <- list(engine_name = NA_character_, call = NULL, formula = NULL)
  if (!is.list(meta)) meta <- list()
  for (nm in names(meta_defaults)) if (is.null(meta[[nm]])) meta[[nm]] <- meta_defaults[[nm]]

  if (!is.list(diagnostics)) diagnostics <- list()
  if (!is.list(extra)) extra <- list()

  result <- list(
    y_hat = y_hat,
    estimate_name = estimate_name,
    se = se,
    converged = converged,
    model = model,
    weights_info = weights_info,
    sample = sample,
    inference = inference,
    diagnostics = diagnostics,
    meta = meta,
    extra = extra
  )
  structure(result, class = c(class_name, "nmar_result"))
}

# Fallback definition for the `%||%` helper used across the S3 stack
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
