#' IID backend
#'
#' @keywords internal
#' @noRd
il_backend_iid <- function(data) {
  list(
    name = "iid",
    vars = data,
    weights_full = NULL,
    fit_mu = function(spec, standardize) {
      il_fit_mu_iid_core(spec = spec, standardize = standardize)
    },
    fit_resp = function(spec, mu_hat, standardize, glm_control = stats::glm.control()) {
      res <- il_fit_resp_iid_core(
        spec = spec,
        mu_hat = mu_hat,
        standardize = standardize,
        glm_control = glm_control
      )
      list(
        induced_glm = res$induced_glm,
        beta_glm = il_unquote_backticked_names(res$beta_glm),
        vcov = res$vcov,
        fitted = res$fitted,
        gamma_hat_paper = res$gamma_hat_paper,
        x1_recipe = res$x1_recipe %||% NULL,
        warnings = res$warnings %||% character(),
        identifiability = res$identifiability %||% list()
      )
    }
  )
}
