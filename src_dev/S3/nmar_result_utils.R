#' Internal helpers for nmar_result objects
#'
#' @keywords internal
nmar_result_get_estimate <- function(x) {
  x$estimate %||% x$y_hat %||% NA_real_
}

#' @keywords internal
nmar_result_get_std_error <- function(x) {
  x$std_error %||% x$se %||% NA_real_
}

#' @keywords internal
nmar_result_get_estimate_name <- function(x) {
  x$estimate_name %||% x$data_info$outcome_var %||%
    (if (!is.null(names(x$estimate)) && length(x$estimate) == 1) names(x$estimate) else NULL) %||%
    "estimand"
}

#' @keywords internal
nmar_result_get_sample <- function(x) {
  sample <- x$sample %||% list()
  sample$n_total <- sample$n_total %||% x$data_info$nobs %||% NA_integer_
  sample$n_respondents <- sample$n_respondents %||% x$data_info$nobs_resp %||% NA_integer_
  sample$is_survey <- sample$is_survey %||% isTRUE(x$data_info$is_survey)
  sample$design <- sample$design %||% x$data_info$design %||% NULL
  sample
}

#' @keywords internal
nmar_result_get_inference <- function(x) {
  inference <- x$inference %||% list()
  inference$variance_method <- inference$variance_method %||% x$data_info$variance_method %||% NA_character_
  inference$df <- inference$df %||% NA_real_
  inference$message <- inference$message %||% NA_character_
  inference$used_pseudoinverse <- inference$used_pseudoinverse %||% isTRUE(x$diagnostics$used_pseudoinverse)
  inference$used_ridge <- inference$used_ridge %||% isTRUE(x$diagnostics$used_ridge)
  inference
}

#' @keywords internal
nmar_result_get_weights_info <- function(x) {
  weights_info <- x$weights_info %||% list(values = x$weights, trimmed_fraction = NA_real_)
  weights_info$values <- weights_info$values %||% x$weights
  weights_info$trimmed_fraction <- weights_info$trimmed_fraction %||% attr(x$weights, "trimmed_fraction") %||% NA_real_
  weights_info
}

#' @keywords internal
nmar_result_get_diagnostics <- function(x) {
  x$diagnostics %||% list()
}

#' @keywords internal
nmar_result_get_model <- function(x) {
  model <- x$model %||% list()
  model$coefficients <- model$coefficients %||% x$coefficients %||% NULL
  model$vcov <- model$vcov %||% x$vcov %||% NULL
  model
}
