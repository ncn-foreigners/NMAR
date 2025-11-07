#' EL internal math utilities (engine-local)
#' @name el_utils_math
#' @keywords internal
#' @noRd
NULL

#' Compute lambda_W from C_const and W
#' @param C_const numeric scalar: (N_pop / sum(d_resp) - 1)
#' @param W numeric scalar in (0,1)
#' @keywords internal
el_lambda_W <- function(C_const, W) {
  W <- min(max(as.numeric(W)[1], 1e-12), 1 - 1e-12)
  C_const / (1 - W)
}

#' Build denominator and floor pack
#' @param lambda_W numeric scalar
#' @param W numeric scalar in (0,1)
#' @param Xc_lambda numeric vector (X_centered \%*\% lambda_x) or 0
#' @param p_i numeric vector of response probabilities
#' @param floor numeric scalar > 0, denominator floor
#' @return list with denom, active, inv, inv_sq
#' @keywords internal
el_denominator <- function(lambda_W, W, Xc_lambda, p_i, floor) {
  denom <- 1 + lambda_W * (p_i - W)
  if (length(Xc_lambda) > 1L || (is.numeric(Xc_lambda) && Xc_lambda[1] != 0)) {
    denom <- denom + as.numeric(Xc_lambda)
  }
  active <- as.numeric(denom > floor)
  denom_guard <- pmax(denom, floor)
  inv <- 1 / denom_guard
  list(denom = denom_guard, active = active, inv = inv, inv_sq = inv * inv)
}

#' EL masses and probabilities from denominators
#' @param weights numeric respondent base weights (d_i)
#' @param denom numeric denominators Di after floor guard
#' @param floor numeric small positive guard (unused in core logic here, kept for API symmetry)
#' @param trim_cap numeric cap (>0) or Inf
#' @return list with mass_untrim, mass_trimmed, prob_mass, trimmed_fraction
#' @keywords internal
el_masses <- function(weights, denom, floor, trim_cap) {
  mass_untrim <- as.numeric(weights) / as.numeric(denom)
# enforce nonnegativity softly then trim
  nn <- enforce_nonneg_weights(mass_untrim)
  if (!nn$ok) stop(nn$message, call. = FALSE)
  mass_untrim <- nn$weights
  trim_res <- trim_weights(mass_untrim, cap = trim_cap)
  mass_trim <- trim_res$weights
  total <- sum(mass_trim)
  prob_mass <- if (total > 0) mass_trim / total else rep(NA_real_, length(mass_trim))
  list(mass_untrim = mass_untrim,
       mass_trimmed = mass_trim,
       prob_mass = prob_mass,
       trimmed_fraction = trim_res$trimmed_fraction)
}

#' Mean from probability masses
#' @keywords internal
el_mean <- function(prob_mass, y) {
  sum(prob_mass * as.numeric(y))
}
