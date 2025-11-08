#' Construct Result Object (parent helper)
#'
#' Builds a normalized `nmar_result` list using the modern schema.
#' Engines must pass named fields; no legacy positional signature is supported.
#'
#' All fields are guaranteed to exist with proper defaults after construction,
#' so downstream code can access fields without defensive NULL checks.
#'
#' @keywords internal
new_nmar_result <- function(...) {
  dots <- list(...)

# Core scalars - all must have defaults
  y_hat <- dots$estimate # Required, will be validated
  se <- dots$se %||% NA_real_
  estimate_name <- dots$estimate_name %||% NA_character_
  converged <- dots$converged %||% FALSE

# Model component - ensure list with required sub-fields
  model <- dots$model %||% list()
  if (!is.list(model)) model <- list()
  model$coefficients <- model$coefficients %||% NULL
  model$vcov <- model$vcov %||% NULL

# Weights component - ensure list with required sub-fields
  weights_info <- dots$weights_info %||% list()
  if (!is.list(weights_info)) weights_info <- list()
  weights_info$values <- weights_info$values %||% NULL
  weights_info$trimmed_fraction <- weights_info$trimmed_fraction %||% NA_real_

# Sample metadata - ensure all required fields exist
  sample <- dots$sample %||% list()
  if (!is.list(sample)) sample <- list()
  sample$n_total <- sample$n_total %||% NA_integer_
  sample$n_respondents <- sample$n_respondents %||% NA_integer_
  sample$is_survey <- sample$is_survey %||% FALSE
  sample$design <- sample$design %||% NULL

# Inference metadata - ensure all required fields exist
  inference <- dots$inference %||% list()
  if (!is.list(inference)) inference <- list()
  inference$variance_method <- inference$variance_method %||% NA_character_
  inference$df <- inference$df %||% NA_real_
  inference$message <- inference$message %||% NA_character_

# Meta information - ensure all required fields exist
  meta <- dots$meta %||% list()
  if (!is.list(meta)) meta <- list()
  meta$engine_name <- meta$engine_name %||% NA_character_
  meta$call <- meta$call %||% NULL
  meta$formula <- meta$formula %||% NULL

# Diagnostics and extra - ensure lists exist (content is engine-specific)
  diagnostics <- dots$diagnostics %||% list()
  if (!is.list(diagnostics)) diagnostics <- list()

  extra <- dots$extra %||% list()
  if (!is.list(extra)) extra <- list()

# Class name
  class_name <- dots$class %||% "nmar_result"

# Construct result object with all fields guaranteed to exist
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
