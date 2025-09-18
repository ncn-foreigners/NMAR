#' Construct Result Object (parent helper)
#'
#' Supports both the legacy positional signature (`y_hat`, `se`, ...)
#' and the new schema-based signature introduced for estimator-agnostic
#' results. Engines will migrate to the new form, but the function keeps
#' the compatibility layer so intermediate refactors do not break.
#'
#' @keywords internal
new_nmar_result <- function(...) {
  dots <- list(...)

  # Detect legacy usage: unnamed first argument assumed to be `y_hat`
  legacy_names <- c("y_hat", "se", "weights", "coefficients", "vcov", "converged", "class")
  legacy_call <-
    (length(dots) >= 6) &&
      (
        is.null(names(dots)) ||
          all(names(dots) %in% legacy_names)
      ) &&
      is.null(dots$estimate)

  if (legacy_call) {
    y_hat <- dots[[1]]
    se <- dots[[2]]
    weights_vec <- dots[[3]]
    coefficients <- dots[[4]]
    vcov <- dots[[5]]
    converged <- dots[[6]]
    class_name <- dots$class %||% dots[[7]] %||% "nmar_result"

    estimate <- y_hat
    estimate_name <- NA_character_
    std_error <- se
    model <- list(coefficients = coefficients, vcov = vcov)
    weights_info <- list(values = weights_vec, trimmed_fraction = NA_real_)
    sample <- list(n_total = NA_integer_, n_respondents = NA_integer_, is_survey = FALSE, design = NULL)
    inference <- list(
      variance_method = NA_character_,
      df = NA_real_,
      message = NA_character_,
      used_pseudoinverse = FALSE,
      used_ridge = FALSE
    )
    diagnostics <- list()
    meta <- list(engine_name = NA_character_, call = NULL, formula = NULL)
    extra <- list()
  } else {
    estimate <- dots$estimate
    std_error <- dots$std_error
    estimate_name <- dots$estimate_name %||% NA_character_
    converged <- dots$converged
    model <- dots$model %||% list()
    weights_info <- dots$weights_info %||% dots$weights %||% list()
    sample <- dots$sample %||% list()
    inference <- dots$inference %||% list()
    diagnostics <- dots$diagnostics %||% list()
    meta <- dots$meta %||% list()
    extra <- dots$extra %||% list()
    class_name <- dots$class %||% "nmar_result"

    # Coerce plain vectors into the list structure.
    if (!is.list(weights_info) || is.null(names(weights_info))) {
      weights_info <- list(values = weights_info, trimmed_fraction = NA_real_)
    }
    if (is.null(weights_info$values) && !is.null(dots$weights) && !is.list(dots$weights)) {
      weights_info$values <- dots$weights
    }
  }

  # Normalise optional components
  if (is.null(model$coefficients)) model$coefficients <- NULL
  if (is.null(model$vcov)) model$vcov <- NULL

  if (is.null(weights_info$values)) weights_info$values <- NULL
  if (is.null(weights_info$trimmed_fraction)) weights_info$trimmed_fraction <- NA_real_

  sample_defaults <- list(n_total = NA_integer_, n_respondents = NA_integer_, is_survey = FALSE, design = NULL)
  for (nm in names(sample_defaults)) {
    if (is.null(sample[[nm]])) sample[[nm]] <- sample_defaults[[nm]]
  }

  inference_defaults <- list(
    variance_method = NA_character_,
    df = NA_real_,
    message = NA_character_,
    used_pseudoinverse = FALSE,
    used_ridge = FALSE
  )
  for (nm in names(inference_defaults)) {
    if (is.null(inference[[nm]])) inference[[nm]] <- inference_defaults[[nm]]
  }

  meta_defaults <- list(engine_name = NA_character_, call = NULL, formula = NULL)
  for (nm in names(meta_defaults)) {
    if (is.null(meta[[nm]])) meta[[nm]] <- meta_defaults[[nm]]
  }

  diagnostics <- diagnostics %||% list()
  extra <- extra %||% list()

  # Build result
  result <- list(
    estimate = estimate,
    estimate_name = estimate_name,
    std_error = std_error,
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
