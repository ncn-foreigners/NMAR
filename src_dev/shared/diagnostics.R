#' Diagnostics bundle for NMAR results
#'
#' Provides a lightweight, standardized container for diagnostic fields that
#' engines can populate and downstream S3 can consume consistently.
#'
#' @keywords internal
new_nmar_diagnostics <- function(x = list()) {
  stopifnot(is.list(x))
  x <- validate_nmar_diagnostics(x)
  structure(x, class = "nmar_diagnostics")
}

#' @keywords internal
validate_nmar_diagnostics <- function(x) {
  # Ensure required names exist with reasonable defaults
  ensure <- function(name, value) if (is.null(x[[name]])) x[[name]] <<- value

  # Solver / convergence
  ensure("convergence_code", NA_integer_)
  ensure("message", NA_character_)
  ensure("solver_jacobian", NA_character_)
  ensure("solver_method", NA_character_)
  ensure("solver_iterations", NA_integer_)
  ensure("reparam_W", NA_character_)

  # Equation diagnostics
  ensure("max_equation_residual", NA_real_)
  ensure("min_denominator", NA_real_)
  ensure("fraction_small_denominators", NA_real_)
  ensure("constraint_sum_W", NA_real_)
  ensure("constraint_sum_aux", NA_real_)

  # Trimming
  ensure("trimmed_fraction", NA_real_)

  # Jacobian diagnostics (variance path)
  ensure("jacobian_source", NA_character_)
  ensure("jacobian_condition_number", NA_real_)
  ensure("jacobian_rel_diff", NA_real_)
  ensure("jacobian_auto_rule", NA_character_)

  # Inversion diagnostics (variance path)
  ensure("invert_rule", NA_character_)
  ensure("used_pseudoinverse", FALSE)
  ensure("used_ridge", FALSE)
  ensure("vcov_message", NA_character_)

  x
}
