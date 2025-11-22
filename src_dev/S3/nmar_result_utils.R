#' Internal helpers for nmar_result objects
#'
#' @keywords internal
nmar_result_get_estimate <- function(x) {
  x$y_hat %||% NA_real_
}



#' @keywords internal
nmar_result_get_estimate_name <- function(x) {
  x$estimate_name %||%
    (if (!is.null(names(x$y_hat)) && length(x$y_hat) == 1) names(x$y_hat) else NULL) %||%
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

#' @keywords internal
nmar_result_get_se <- function(x) {
  x$se %||% NA_real_
}

#' Resolve global digits setting for printing
#' @keywords internal
nmar_get_digits <- function() {
  d <- getOption("nmar.digits", 6L)
  if (!is.numeric(d) || length(d) != 1L || is.na(d) || d < 0) return(6L)
  as.integer(d)
}

#' Format a number with fixed decimal places using nmar.digits
#' @keywords internal
nmar_fmt_num <- function(x, digits = nmar_get_digits()) {
  x <- as.numeric(x)
  if (length(x) == 0) return(character(0))
  fin <- is.finite(x)
  out <- character(length(x))
  out[fin] <- sprintf(paste0("%0.", digits, "f"), x[fin])
  out[!fin] <- "NA"
  out
}

#' Format an abridged call line for printing
#'
#' Builds a concise one-line summary of the original call without
#' materializing large objects (e.g., full data frames). Intended for
#' use by print/summary methods.
#'
#' Uses option `nmar.show_call` (default TRUE). Width can be tuned via
#' option `nmar.call_width` (default 120), but the formatter aims to keep
#' the line compact regardless of width.
#'
#' @keywords internal
nmar_format_call_line <- function(x) {
  if (!isTRUE(getOption("nmar.show_call", TRUE))) return(NULL)
  meta <- x$meta %||% list()
  sample <- nmar_result_get_sample(x)
# Formula
  fml <- meta$formula %||% NULL
  fml_str <- if (!is.null(fml) && inherits(fml, "formula")) {
    paste(deparse(fml, width.cutoff = max(60L, getOption("nmar.call_width", 120L))), collapse = " ")
  } else {
    "<formula>"
  }
# Data descriptor (avoid printing the object itself). Display N using a
# rounded/cleaned representation to avoid showing floating-point noise.
  N_raw <- sample$n_total
  if (is.finite(N_raw)) {
# Treat N as integer when it is very close to one, using a relative
# tolerance to avoid spurious decimals from floating-point arithmetic.
    if (abs(N_raw - round(N_raw)) < 1e-6 * max(1, abs(N_raw))) {
      n_str <- paste0("N=", as.integer(round(N_raw)))
    } else {
      formatted <- formatC(N_raw, digits = 6, format = "fg", flag = "#")
      formatted <- sub("^\\s+", "", formatted)
      n_str <- paste0("N=", formatted)
    }
  } else {
    n_str <- "N=?"
  }
  data_desc <- if (isTRUE(sample$is_survey)) paste0("<survey.design: ", n_str, ">") else paste0("<data.frame: ", n_str, ">")
# Engine label
  eng <- meta$engine_name %||% "nmar_engine"
  sprintf("Call: nmar(%s, data = %s, engine = %s)", fml_str, data_desc, eng)
}
