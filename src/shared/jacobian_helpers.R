#' Choose Jacobian source (analytic vs numeric) and compare
#' @keywords internal
choose_jacobian <- function(analytic_fun, numeric_fun, at) {
  A_analytic <- NULL
  A_numeric <- NULL
  source <- NULL
  kappa_val <- NA_real_
  rel_diff <- NA_real_
  if (!is.null(analytic_fun)) {
    A_analytic <- tryCatch(analytic_fun(at), error = function(e) NULL)
  }
  A_numeric <- tryCatch(numeric_fun(at), error = function(e) NULL)
  if (!is.null(A_analytic) && !is.null(A_numeric)) {
    denom <- max(1e-8, norm(A_numeric, type = "F"))
    rel_diff <- tryCatch(norm(A_analytic - A_numeric, type = "F") / denom, error = function(e) NA_real_)
  }
  if (!is.null(A_analytic)) {
    source <- "analytic"
    A <- A_analytic
  } else {
    source <- "numeric"
    A <- A_numeric
  }
  kappa_val <- tryCatch(if (!is.null(A)) kappa(A) else NA_real_, error = function(e) NA_real_)
  list(A = A, A_analytic = A_analytic, A_numeric = A_numeric, source = source, kappa = kappa_val, rel_diff = rel_diff)
}
