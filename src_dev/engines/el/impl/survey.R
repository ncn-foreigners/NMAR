#' Empirical likelihood for survey designs (NMAR)
#' @description Method for `el()` when `data` is a `survey.design`. Uses design‑based
#'   covariance for variance estimation when `variance_method = 'delta'`.
#' @param data A `survey.design` created with [survey::svydesign()].
#' @param formula Two-sided formula: NA-valued outcome on LHS; auxiliaries on RHS.
#' @param response_predictors Optional character vector for the response model RHS. These may include variables not on
#'   the RHS of the outcome formula; such variables will enter only the response model (no auxiliary constraint).
#' @param auxiliary_means Named numeric vector of population means for auxiliaries.
#' @param standardize Logical; standardize predictors.
#' @param trim_cap Numeric; cap for EL weights (Inf = no trimming).
#' @param control List; solver control.
#' @param on_failure Character; "return" or "error" on solver failure.
#' @param variance_method Character; "delta" or "bootstrap".
#' @param variance_jacobian Character; "auto", "analytic", or "numeric".
#' @param solver_jacobian Character; "auto", "analytic", or "none".
#' @param variance_pseudoinverse Logical; allow pseudo-inverse for variance.
#' @param bootstrap_reps Integer; reps when `variance_method = "bootstrap"`.
#' @param suppress_warnings Logical; suppress variance method warnings.
#' @param ... Passed to solver.
#' @param variance_ridge Logical or numeric; if TRUE, apply adaptive ridge in
#'   Jacobian inversion; if numeric, treated as ridge epsilon.
#' @details
#' Implements the empirical likelihood estimator of Qin, Leung and Shao (2002)
#' with design weights and design‑based covariance of total scores for delta
#' variance. Bootstrap variance via replicate weights is also supported.
#' @references
#' Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American
#' Statistical Association, 97(457), 193–200.
#' @return `c('nmar_result_el','nmar_result')`.
#' @export
el.survey.design <- function(data, formula, response_predictors = NULL,
                             auxiliary_means = NULL, standardize = TRUE,
                             trim_cap = Inf, control = list(),
                             on_failure = c("return", "error"),
                             variance_method = c("delta", "bootstrap"),
                             variance_jacobian = c("auto", "analytic", "numeric"),
                             solver_jacobian = c("auto", "analytic", "none"),
                             variance_pseudoinverse = FALSE, variance_ridge = FALSE,
                             bootstrap_reps = 500, suppress_warnings = FALSE, ...) {
  cl <- match.call()
  on_failure <- match.arg(on_failure)
  variance_method <- match.arg(variance_method)
  variance_jacobian <- match.arg(variance_jacobian)
  solver_jacobian <- match.arg(solver_jacobian)

  design <- data

  parsed_inputs <- prepare_el_inputs(formula, design$variables, response_predictors)
  design$variables <- parsed_inputs$data
  internal_formula <- parsed_inputs$formula_list
  response_var <- all.vars(internal_formula$response)[1]
  observed_mask <- design$variables[[response_var]] == 1
  observed_indices <- which(observed_mask)
  resp_design <- subset(design, observed_mask)
  respondent_weights <- weights(resp_design)
  N_pop <- sum(weights(design))

  compute_score_covariance_func_survey <- function(U_matrix_resp, full_design) {
    U_full <- matrix(0, nrow = nrow(full_design$variables), ncol = ncol(U_matrix_resp))
    score_variables <- paste0("..score_", seq_len(ncol(U_matrix_resp)))
    colnames(U_full) <- score_variables
    U_full[observed_indices, ] <- U_matrix_resp
    design_scores <- full_design
    design_scores$variables <- cbind(design_scores$variables, as.data.frame(U_full))
    score_formula <- stats::as.formula(paste("~", paste(score_variables, collapse = " + ")))
    score_totals <- survey::svytotal(score_formula, design_scores, na.rm = TRUE)
    stats::vcov(score_totals)
  }

  user_args <- list(
    formula = formula, response_predictors = response_predictors,
    auxiliary_means = auxiliary_means, standardize = standardize,
    trim_cap = trim_cap, control = control,
    suppress_warnings = suppress_warnings, ...
  )

  core_results <- el_estimator_core(
    full_data = design,
    respondent_data = resp_design$variables,
    respondent_weights = respondent_weights, N_pop = N_pop,
    internal_formula = internal_formula, auxiliary_means = auxiliary_means,
    standardize = standardize, trim_cap = trim_cap, control = control,
    compute_score_variance_func = compute_score_covariance_func_survey, on_failure = on_failure,
    variance_method = variance_method, bootstrap_reps = bootstrap_reps,
    variance_jacobian = variance_jacobian, solver_jacobian = solver_jacobian,
    variance_pseudoinverse = variance_pseudoinverse, variance_ridge = variance_ridge, user_args = user_args, ...
  )

  sample_info <- list(
    outcome_var = all.vars(internal_formula$outcome)[1],
    response_var = response_var,
    formula = formula,
    nobs = nrow(design$variables),
    nobs_resp = length(observed_indices),
    is_survey = TRUE,
    design = design,
    variance_method = variance_method
  )

  if (!core_results$converged) {
    diag_list <- core_results$diagnostics
    if (is.null(diag_list)) diag_list <- list()
    msg <- diag_list$message
    if (is.null(msg)) msg <- NA_character_
    result <- new_nmar_result(
      estimate = NA_real_,
      estimate_name = sample_info$outcome_var,
      std_error = NA_real_,
      converged = FALSE,
      model = list(coefficients = NULL, vcov = NULL),
      weights_info = list(values = numeric(0), trimmed_fraction = NA_real_),
      sample = list(n_total = sample_info$nobs, n_respondents = sample_info$nobs_resp, is_survey = TRUE, design = design),
      inference = list(variance_method = variance_method, df = NA_real_, message = msg, used_pseudoinverse = FALSE, used_ridge = FALSE),
      diagnostics = diag_list,
      meta = list(engine_name = "empirical_likelihood", call = cl, formula = formula),
      extra = list(nmar_scaling_recipe = core_results$nmar_scaling_recipe),
      class = "nmar_result_el"
    )
    return(validate_nmar_result(result, 'nmar_result_el'))
  }

  result <- new_nmar_result_el(
    y_hat = core_results$y_hat, se = core_results$se, weights = core_results$weights,
    coefficients = core_results$coefficients, vcov = core_results$vcov,
    converged = core_results$converged, diagnostics = core_results$diagnostics,
    data_info = sample_info,
    nmar_scaling_recipe = core_results$nmar_scaling_recipe,
    fitted_values = core_results$fitted_values, call = cl
  )

  validate_nmar_result(result, 'nmar_result_el')
}
