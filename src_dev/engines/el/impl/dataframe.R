#' Empirical likelihood for data frames (NMAR)
#' @description Internal method dispatched by `el()` when `data` is a `data.frame`.
#'   Returns `c('nmar_result_el','nmar_result')` with estimate, standard error, weights,
#'   coefficients, diagnostics and metadata.
#' @param data A `data.frame` where the outcome column contains `NA` for nonrespondents.
#' @param formula Two-sided formula `Y_miss ~ auxiliaries`.
#' @param response_predictors Optional character vector naming predictors for the response (missingness) model.
#'   These variables enter only the response model (no auxiliary constraint) and do not require population means.
#' @param auxiliary_means Named numeric vector of population means for auxiliary variables (names must match RHS of outcome formula).
#' @param standardize Logical; whether to standardize predictors prior to estimation.
#' @param trim_cap Numeric; cap for EL weights (`Inf` = no trimming).
#' @param control List; optional solver control parameters for `nleqslv(control=...)`.
#' @param solver_args List; optional top-level `nleqslv` args (e.g., `global`, `xscalm`).
#' @param on_failure Character; one of `"return"` or `"error"` on solver failure.
#' @param variance_method Character; one of `"delta"`, `"bootstrap"`, or `"none"`.
#' @param variance_jacobian Character; one of `"auto"`, `"analytic"`, or `"numeric"`.
#' @param solver_jacobian Character; one of `"auto"`, `"analytic"`, or `"none"`.
#' @param variance_pseudoinverse Logical; allow pseudo-inverse for variance when needed.
#' @param bootstrap_reps Integer; number of bootstrap reps if `variance_method = "bootstrap"`.
#' @param suppress_warnings Logical; suppress variance-related warnings.
#' @param variance_ridge Logical or numeric; if TRUE, apply adaptive ridge in
#'   Jacobian inversion; if numeric, treated as ridge epsilon.
#' @param ... Additional arguments passed to the solver.
#' @details
#' Implements the empirical likelihood estimator of Qin, Leung and Shao (2002)
#' for IID data. The response‑model score is the derivative of the Bernoulli
#' log‑likelihood with respect to the linear predictor; auxiliary moment constraints
#' are optional.
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193–200.
#' @keywords internal
el.data.frame <- function(data, formula, response_predictors = NULL,
                          auxiliary_means = NULL, standardize = TRUE,
                          trim_cap = Inf, control = list(),
                          on_failure = c("return", "error"),
                          variance_method = c("delta", "bootstrap", "none"),
                          variance_jacobian = c("auto", "analytic", "numeric"),
                          solver_jacobian = c("auto", "analytic", "none"),
                          solver_method = c("auto", "newton", "broyden"),
                          solver_args = list(),
                          variance_pseudoinverse = FALSE, variance_ridge = FALSE,
                          bootstrap_reps = 500, suppress_warnings = FALSE, ...) {
  cl <- match.call()
  on_failure <- match.arg(on_failure)
  if (is.null(variance_method)) variance_method <- "none"
  variance_method <- match.arg(variance_method)
  variance_jacobian <- match.arg(variance_jacobian)
  solver_jacobian <- match.arg(solver_jacobian)
  solver_method <- match.arg(solver_method)

  parsed_inputs <- prepare_el_inputs(formula, data, response_predictors)
  estimation_data <- parsed_inputs$data
  internal_formula <- parsed_inputs$formula_list
  response_var <- all.vars(internal_formula$response)[1]
  observed_indices <- which(estimation_data[[response_var]] == 1)

  respondent_weights <- rep(1, length(observed_indices))
  N_pop <- nrow(estimation_data)

  compute_score_covariance_func_df <- function(U_matrix_resp, full_df) {
    U_full <- matrix(0, nrow = nrow(full_df), ncol = ncol(U_matrix_resp))
    colnames(U_full) <- colnames(U_matrix_resp)
    U_full[observed_indices, ] <- U_matrix_resp
    crossprod(U_full)
  }

  user_args <- list(
    formula = formula, response_predictors = response_predictors,
    auxiliary_means = auxiliary_means, standardize = standardize,
    trim_cap = trim_cap, control = control,
    suppress_warnings = suppress_warnings, ...
  )

  core_results <- el_estimator_core(
    full_data = estimation_data,
    respondent_data = estimation_data[observed_indices, ],
    respondent_weights = respondent_weights, N_pop = N_pop,
    internal_formula = internal_formula, auxiliary_means = auxiliary_means,
    standardize = standardize, trim_cap = trim_cap, control = control,
    compute_score_variance_func = compute_score_covariance_func_df, on_failure = on_failure,
    variance_method = variance_method, bootstrap_reps = bootstrap_reps,
    variance_jacobian = variance_jacobian, solver_jacobian = solver_jacobian,
    solver_method = solver_method, solver_args = solver_args,
    variance_pseudoinverse = variance_pseudoinverse, variance_ridge = variance_ridge, user_args = user_args, ...
  )

  sample_info <- list(
    outcome_var = all.vars(internal_formula$outcome)[1],
    response_var = response_var,
    formula = formula,
    nobs = nrow(estimation_data),
    nobs_resp = length(observed_indices),
    is_survey = FALSE,
    design = NULL,
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
      se = NA_real_,
      converged = FALSE,
      model = list(coefficients = NULL, vcov = NULL),
      weights_info = list(values = numeric(0), trimmed_fraction = NA_real_),
      sample = list(n_total = sample_info$nobs, n_respondents = sample_info$nobs_resp, is_survey = FALSE, design = NULL),
      inference = list(variance_method = variance_method, df = NA_real_, message = msg, used_pseudoinverse = FALSE, used_ridge = FALSE),
      diagnostics = diag_list,
      meta = list(engine_name = "empirical_likelihood", call = cl, formula = formula),
      extra = list(nmar_scaling_recipe = core_results$nmar_scaling_recipe),
      class = "nmar_result_el"
    )
    return(validate_nmar_result(result, "nmar_result_el"))
  }

  result <- new_nmar_result_el(
    y_hat = core_results$y_hat, se = core_results$se, weights = core_results$weights,
    coefficients = core_results$coefficients, vcov = core_results$vcov,
    converged = core_results$converged, diagnostics = core_results$diagnostics,
    data_info = sample_info,
    nmar_scaling_recipe = core_results$nmar_scaling_recipe,
    fitted_values = core_results$fitted_values, call = cl
  )

  validate_nmar_result(result, "nmar_result_el")
}
