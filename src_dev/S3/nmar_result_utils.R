#' Internal helpers for nmar_result objects
#'
#' @keywords internal
nmar_result_get_estimate <- function(x) {
  x$estimate %||% NA_real_
}

#' @keywords internal
nmar_result_get_std_error <- function(x) {
  x$std_error %||% NA_real_
}

#' @keywords internal
nmar_result_get_estimate_name <- function(x) {
  x$estimate_name %||%
    (if (!is.null(names(x$estimate)) && length(x$estimate) == 1) names(x$estimate) else NULL) %||%
    "estimand"
}

#' @keywords internal
nmar_result_get_sample <- function(x) {
  sample <- x$sample %||% list()
  sample$n_total <- sample$n_total %||% NA_integer_
  sample$n_respondents <- sample$n_respondents %||% NA_integer_
  sample$is_survey <- sample$is_survey %||% FALSE
  sample$design <- sample$design %||% NULL
  sample
}

#' @keywords internal
nmar_result_get_inference <- function(x) {
  inference <- x$inference %||% list()
  inference$variance_method <- inference$variance_method %||% NA_character_
  inference$df <- inference$df %||% NA_real_
  inference$message <- inference$message %||% NA_character_
  inference$used_pseudoinverse <- inference$used_pseudoinverse %||% FALSE
  inference$used_ridge <- inference$used_ridge %||% FALSE
  inference
}

#' @keywords internal
nmar_result_get_weights_info <- function(x) {
  weights_info <- x$weights_info %||% list(values = NULL, trimmed_fraction = NA_real_)
  weights_info$values <- weights_info$values %||% NULL
  weights_info$trimmed_fraction <- weights_info$trimmed_fraction %||% NA_real_
  weights_info
}

#' @keywords internal
nmar_result_get_diagnostics <- function(x) {
  x$diagnostics %||% list()
}

#' @keywords internal
nmar_result_get_model <- function(x) {
  model <- x$model %||% list()
  model$coefficients <- model$coefficients %||% NULL
  model$vcov <- model$vcov %||% NULL
  model
}
