#' Data/metadata bundle for NMAR results
#'
#' Standardizes commonly used metadata so that engines report a uniform shape
#' and downstream S3 methods (tidy/glance/confint) can rely on it.
#'
#' @keywords internal
new_nmar_data_info <- function(x = list()) {
  stopifnot(is.list(x))
  x <- validate_nmar_data_info(x)
  structure(x, class = "nmar_data_info")
}

#' @keywords internal
validate_nmar_data_info <- function(x) {
  ensure <- function(name, value) if (is.null(x[[name]])) x[[name]] <<- value

  ensure("outcome_var", NA_character_)
  ensure("response_var", NA_character_)
  ensure("formula", NULL)
  ensure("nobs", NA_integer_)
  ensure("nobs_resp", NA_integer_)
  ensure("is_survey", FALSE)
  # `design` can be NULL or a survey.design; we do not enforce a class here
  ensure("design", NULL)
  ensure("variance_method", NA_character_)

  x
}
