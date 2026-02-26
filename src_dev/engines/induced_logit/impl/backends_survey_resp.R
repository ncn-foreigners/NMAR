#' Survey backend: missingness-model fits
#'
#' @keywords internal
#' @noRd
il_fit_resp_survey_core <- function(design, spec, mu_hat, w_full, standardize, glm_control) {
  il_assert_scalar_logical(standardize, name = "standardize")

  if (is.null(glm_control)) glm_control <- stats::glm.control()
  validator_assert_list(glm_control, name = "glm_control")

  r_vec <- spec$r_vec %||% as.integer(spec$respondent_mask)
  mu_col <- as.numeric(mu_hat)
  mu_coef_name <- IL_COL_MU_HAT

  if (isTRUE(standardize) && is.null(w_full)) {
    stop("Internal error: `w_full` is required when standardize = TRUE.", call. = FALSE)
  }

# Use cached model.matrix columns so the survey and IID paths use identical predictors
# This avoids re-parsing formulas and reduces the chance of silent divergence
  design_out <- il_build_response_design(
    spec = spec,
    mu_hat = mu_col,
    standardize = isTRUE(standardize),
    weights = if (isTRUE(standardize)) w_full else NULL
  )
  x1_df <- design_out$x1_df
  x1_recipe <- design_out$x1_recipe

  d_glm <- cbind(
    data.frame(setNames(list(r_vec, mu_col), c(IL_COL_R, IL_COL_MU_HAT)), check.names = FALSE),
    x1_df
  )
  design_glm <- design
  design_glm$variables <- d_glm

  rhs_terms <- c(
    if (ncol(x1_df) > 0) vapply(names(x1_df), il_quote_name_for_formula, character(1)) else character(0),
    il_quote_name_for_formula(IL_COL_MU_HAT)
  )
  resp_formula <- il_reformulate(response = IL_COL_R, term_labels = rhs_terms, intercept = TRUE)
  resp_mm <- stats::model.matrix(resp_formula, data = d_glm)

  captured <- il_capture_warnings(
    do.call(
      survey::svyglm,
      list(
        formula = resp_formula,
        design = design_glm,
        family = stats::quasibinomial(),
        control = glm_control
      )
    )
  )
  resp_fit <- captured$value
  beta_raw <- il_unquote_backticked_names(stats::coef(resp_fit))
  vc_raw <- tryCatch(il_unquote_backticked_dimnames(stats::vcov(resp_fit)), error = function(e) NULL)
  gamma_hat_paper <- il_extract_gamma_hat(beta_raw, mu_coef_name = mu_coef_name)
  fitted_vals <- tryCatch(as.numeric(stats::fitted(resp_fit)), error = function(e) numeric(0))

  if (is.null(beta_raw) || anyNA(beta_raw)) {
    stop(
      "Induced-logit response model has NA coefficients (rank deficiency / non-identifiability / separation).",
      call. = FALSE
    )
  }

  ident <- il_enforce_response_model_identifiability(il_response_model_diagnostics_matrix(resp_mm))
  warnings <- captured$warnings
  if (is.character(ident$warning) && nzchar(ident$warning)) warnings <- c(warnings, ident$warning)

  beta_out <- beta_raw
  vc_out <- vc_raw
  if (isTRUE(standardize)) {
    un <- il_unscale_response_model(
      beta_scaled = beta_raw,
      vcov_scaled = vc_raw,
      x1_recipe = x1_recipe,
      mu_coef_name = mu_coef_name
    )
    beta_out <- un$coefficients
    vc_out <- un$vcov
  }

  list(
    induced_glm = resp_fit,
    beta_glm = beta_out,
    vcov = vc_out,
    fitted = fitted_vals,
    gamma_hat_paper = as.numeric(gamma_hat_paper),
    x1_recipe = x1_recipe,
    warnings = unique(warnings),
    identifiability = ident$diag
  )
}
