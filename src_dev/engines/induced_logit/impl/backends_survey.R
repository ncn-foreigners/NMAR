#' Survey backend
#'
#' @keywords internal
#' @noRd
il_backend_survey <- function(design) {
  il_require_survey()
  vars <- design$variables
  w_info <- il_get_survey_weights(design, vars)
  w_full <- w_info$weights
  sum_w <- w_info$sum_weights

  list(
    name = "survey",
    design = design,
    vars = vars,
    weights_full = w_full,
    sum_weights = sum_w,
    fit_mu = function(spec, standardize) {
      il_fit_mu_survey_core(design = design, spec = spec, vars = vars, w_full = w_full, standardize = standardize)
    },
    fit_resp = function(spec, mu_hat, standardize, glm_control = stats::glm.control()) {
      il_fit_resp_survey_core(
        design = design,
        spec = spec,
        mu_hat = mu_hat,
        w_full = if (isTRUE(standardize)) w_full else NULL,
        standardize = standardize,
        glm_control = glm_control
      )
    }
  )
}
