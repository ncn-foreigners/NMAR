#' Survey backend: outcome-model fits
#'
#' @keywords internal
#' @noRd
il_fit_mu_survey_core <- function(design, spec, vars, w_full, standardize) {
  il_assert_scalar_logical(standardize, name = "standardize")

  if (!isTRUE(standardize)) return(il_survey_fit_mu_unscaled(design = design, spec = spec))
  il_survey_fit_mu_scaled(design = design, spec = spec, vars = vars, w_full = w_full)
}

#' @keywords internal
#' @noRd
il_survey_fit_mu_unscaled <- function(design, spec) {
  mu_mat_full <- spec$mu_mat_full %||% stop("Internal error: missing spec$mu_mat_full.", call. = FALSE)
  design_resp <- subset(design, spec$respondent_mask)
  mu_formula <- il_reformulate(
    response = spec$outcome,
    term_labels = spec$mu_term_labels %||% character(0),
    intercept = spec$mu_intercept %||% TRUE
  )
  captured <- il_capture_warnings(survey::svyglm(mu_formula, design = design_resp))
  mu_fit <- captured$value
  mu_coef <- il_unquote_backticked_names(stats::coef(mu_fit))
  mm_cols <- colnames(mu_mat_full)
  il_assert_named_coefficients(mu_coef, mm_cols, "Outcome regression")
  mu_hat <- as.numeric(mu_mat_full %*% mu_coef[mm_cols])
  list(
    mu_fit = new_nmar_il_mu_fit(mu_rhs = spec$mu_rhs, coef = mu_coef[mm_cols], recipe = NULL, fit = mu_fit),
    mu_hat = mu_hat,
    mu_recipe = NULL,
    warnings = captured$warnings
  )
}

#' @keywords internal
#' @noRd
il_survey_fit_mu_scaled <- function(design, spec, vars, w_full) {
  mu_mat_full <- spec$mu_mat_full %||% stop("Internal error: missing spec$mu_mat_full.", call. = FALSE)
  respondent_mask <- spec$respondent_mask
  design_resp <- subset(design, respondent_mask)
  mu_scale <- il_scale_mu_matrix(
    mu_mat_full,
    has_intercept = spec$mu_intercept %||% TRUE,
    weights = w_full,
    weight_mask = respondent_mask
  )
  mu_mat_scaled <- mu_scale$mat
  mu_cols <- setdiff(colnames(mu_mat_scaled), "(Intercept)")

  y_resp <- vars[[spec$outcome]][respondent_mask]
  d_mu_resp <- data.frame(y_resp, check.names = FALSE)
  names(d_mu_resp) <- spec$outcome
  if (length(mu_cols) > 0L) {
    d_mu_resp <- cbind(d_mu_resp, as.data.frame(mu_mat_scaled[respondent_mask, mu_cols, drop = FALSE], check.names = FALSE))
  }
  design_resp_scaled <- design_resp
  design_resp_scaled$variables <- d_mu_resp
  mu_terms <- if (length(mu_cols) > 0L) vapply(mu_cols, il_quote_name_for_formula, character(1)) else character(0)
  mu_formula_scaled <- il_reformulate(
    response = spec$outcome,
    term_labels = mu_terms,
    intercept = spec$mu_intercept %||% TRUE
  )

  captured <- il_capture_warnings(survey::svyglm(mu_formula_scaled, design = design_resp_scaled))
  mu_fit <- captured$value
  mu_coef <- il_unquote_backticked_names(stats::coef(mu_fit))
  mm_cols <- colnames(mu_mat_scaled)
  il_assert_named_coefficients(mu_coef, mm_cols, "Outcome regression")
  mu_hat <- as.numeric(mu_mat_scaled %*% mu_coef[mm_cols])
  list(
    mu_fit = new_nmar_il_mu_fit(mu_rhs = spec$mu_rhs, coef = mu_coef[mm_cols], recipe = mu_scale$recipe, fit = mu_fit),
    mu_hat = mu_hat,
    mu_recipe = mu_scale$recipe,
    warnings = captured$warnings
  )
}
