#' Nonparametric exponential tilting engine
#'
#' Build a configuration for the nonparametric exponential-tilting EM estimator
#' used in NMAR problems with grouped outcomes and refusal counts. Pass the
#' resulting engine to [nmar()] together with the appropriate formula and data.
#'
#' @param refusal_col Column name in `data` containing refusal counts.
#' @param max_iter Maximum number of EM iterations.
#' @param tol_value Convergence tolerance for the EM updates.
#'
#' @return A list of class `c("nmar_engine_exptilt_nonparam", "nmar_engine")`.
#' @keywords engine
#' @export
exptilt_nonparam_engine <- function(
    refusal_col = '',
    max_iter = 100,
    tol_value = 1e-6) {

  validator$assert_positive_integer(max_iter, name = "max_iter")
  validator$assert_number(tol_value, name = "tol_value", min = 0, max = Inf)

  engine <- list(
    refusal_col = refusal_col,
    max_iter = max_iter,
    tol_value = tol_value
  )
  class(engine) <- c("nmar_engine_exptilt_nonparam", "nmar_engine")
  engine
}
