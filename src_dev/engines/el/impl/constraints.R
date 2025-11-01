#' Constraint summaries for EL diagnostics
#' @keywords internal
constraint_summaries <- function(w_i_hat, W_hat, mass_untrim, X_centered) {
# Report the raw constraint sums used in the estimating equations.
# With mass_untrim = respondent_weights / Di, equations enforce (QLS):
#   sum mass_untrim * (w_i - W) = 0  [W equation]
#   and, if auxiliaries present, t(X_centered) %*% mass_untrim = 0  [aux constraints].
  out <- list(constraint_sum_W = sum(mass_untrim * (w_i_hat - W_hat)))
  if (!is.null(X_centered) && ncol(X_centered) > 0) {
    out$constraint_sum_aux <- as.numeric(crossprod(mass_untrim, X_centered))
    names(out$constraint_sum_aux) <- colnames(X_centered)
  } else {
    out$constraint_sum_aux <- numeric(0)
  }
  out
}
