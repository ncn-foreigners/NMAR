#' Validate nmar_result structure
#'
#' Ensures both the child class and the parent schema are satisfied.
#' This function is a PURE validator - it only checks/asserts structure
#' and types, but does NOT modify the object. All default-filling is
#' handled by the constructor (new_nmar_result).
#'
#' @keywords internal
validate_nmar_result <- function(x, class_name) {
  stopifnot(is.list(x), inherits(x, class_name), inherits(x, "nmar_result"))

# Core scalars - assert existence and validate types
  stopifnot("`y_hat` must exist" = !is.null(x$y_hat))
  validator$assert_scalar_numeric(x$y_hat, name = "y_hat", allow_na = TRUE, finite = FALSE)

  stopifnot("`se` must exist" = !is.null(x$se))
  validator$assert_scalar_numeric(x$se, name = "se", allow_na = TRUE, finite = FALSE)

  stopifnot("`estimate_name` must exist" = !is.null(x$estimate_name))
  validator$assert_scalar_character(x$estimate_name, name = "estimate_name", allow_na = TRUE, non_empty = FALSE)

  stopifnot("`converged` must exist" = !is.null(x$converged))
  validator$assert_scalar_logical(x$converged, name = "converged")

# Component lists - assert existence and type
  stopifnot("`model` must be a list" = is.list(x$model))
  stopifnot("`weights_info` must be a list" = is.list(x$weights_info))
  stopifnot("`sample` must be a list" = is.list(x$sample))
  stopifnot("`inference` must be a list" = is.list(x$inference))
  stopifnot("`diagnostics` must be a list" = is.list(x$diagnostics))
  stopifnot("`meta` must be a list" = is.list(x$meta))
  stopifnot("`extra` must be a list" = is.list(x$extra))

# Required sub-fields in component lists
  stopifnot("`sample$n_total` must exist" = !is.null(x$sample$n_total))
  stopifnot("`sample$n_respondents` must exist" = !is.null(x$sample$n_respondents))
  stopifnot("`sample$is_survey` must exist" = !is.null(x$sample$is_survey))

  stopifnot("`inference$variance_method` must exist" = !is.null(x$inference$variance_method))
  stopifnot("`inference$df` must exist" = !is.null(x$inference$df))

  stopifnot("`meta$engine_name` must exist" = !is.null(x$meta$engine_name))

# Return object unchanged (pure validator)
  x
}
