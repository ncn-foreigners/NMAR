#' Variance helper utilities
#'
#' Centralizes Jacobian selection for variance and the delta-variance pipeline
#' to keep the core estimator compact and consistent.
#'
#' Tuning options (set via options()):
#'  - nmar.var_auto_rel_diff_thr: relative Frobenius threshold between analytic
#'    and numeric Jacobians for the "auto" selector (default 1e-3).
#'  - nmar.var_auto_kappa_ratio_thr: condition-number ratio gate to switch from
#'    analytic to numeric in "auto" mode (default 10).
#'  - nmar.var_auto_kappa_ridge_thr: singularity threshold to enable ridge
#'    stabilization in variance inversion (default 1e12).
#'  - nmar.var_auto_kappa_pinv_thr: singularity threshold to enable pseudo-
#'    inverse for variance inversion (default 1e14).
#'  - nmar.var_ridge_base: base for adaptive ridge epsilon (default 1e-6).
#' @references Qin, J., Leung, D., and Shao, J. (2002). Estimation with survey data under
#' nonignorable nonresponse or informative sampling. Journal of the American Statistical Association, 97(457), 193–200.
#' @keywords internal

#' Select Jacobian for variance with diagnostics
#'
#' Applies a consistent policy across EL to choose between analytic and numeric
#' Jacobians for variance, recording diagnostics to surface to users.
#'
#' @param equation_system_func Function mapping parameter vector to equations.
#' @param analytical_jac_func Analytic Jacobian function; may be NULL.
#' @param estimates Numeric vector at the solution.
#' @param variance_jacobian Character; one of "auto", "analytic", or "numeric".
#' @return List with entries: `A_matrix_var`, `A_source`, `A_condition`,
#'   `A_diff_norm`, and `jacobian_auto_rule`.
el_select_variance_jacobian <- function(equation_system_func,
                                        analytical_jac_func,
                                        estimates,
                                        variance_jacobian = c("auto", "analytic", "numeric")) {
  variance_jacobian <- match.arg(variance_jacobian)
  jac_choice <- choose_jacobian(
    analytic_fun = analytical_jac_func,
    numeric_fun = function(x) numDeriv::jacobian(func = equation_system_func, x = x),
    at = estimates
  )
  A_matrix_analytic <- jac_choice$A_analytic
  A_matrix_numeric <- jac_choice$A_numeric
  A_diff_norm <- jac_choice$rel_diff
  A_condition_numeric <- tryCatch(if (!is.null(A_matrix_numeric)) kappa(A_matrix_numeric) else NA_real_, error = function(e) NA_real_)
  A_condition_analytic <- tryCatch(if (!is.null(A_matrix_analytic)) kappa(A_matrix_analytic) else NA_real_, error = function(e) NA_real_)
  jacobian_auto_rule <- "default"
  if (variance_jacobian == "numeric") {
    A_matrix_var <- A_matrix_numeric
    A_condition <- A_condition_numeric
    A_source <- "numeric"
  } else if (variance_jacobian == "analytic") {
    A_matrix_var <- A_matrix_analytic
    A_condition <- A_condition_analytic
    A_source <- "analytic"
  } else {
# auto: prefer analytic when available, unless quality gates suggest numeric
    REL_THR <- getOption("nmar.var_auto_rel_diff_thr", 1e-3)
    KAPPA_RATIO_THR <- getOption("nmar.var_auto_kappa_ratio_thr", 10)
    if (!is.null(A_matrix_analytic)) {
# Start with analytic as default
      A_matrix_var <- A_matrix_analytic
      A_condition <- A_condition_analytic
      A_source <- "analytic"
# If numeric also available, evaluate quality gates
      if (!is.null(A_matrix_numeric)) {
# Gate 1: large relative difference between analytic and numeric
        if (is.finite(A_diff_norm) && A_diff_norm > REL_THR) {
          A_matrix_var <- A_matrix_numeric
          A_condition <- A_condition_numeric
          A_source <- "numeric"
          jacobian_auto_rule <- "rel_diff_high"
        } else {
# Gate 2: analytic much worse conditioned than numeric
          if (is.finite(A_condition_analytic) && is.finite(A_condition_numeric) &&
              A_condition_analytic > KAPPA_RATIO_THR * A_condition_numeric) {
            A_matrix_var <- A_matrix_numeric
            A_condition <- A_condition_numeric
            A_source <- "numeric"
            jacobian_auto_rule <- "kappa_ratio_high"
          }
        }
      }
    } else if (!is.null(A_matrix_numeric)) {
      A_matrix_var <- A_matrix_numeric
      A_condition <- A_condition_numeric
      A_source <- "numeric"
    } else {
# Fallback to whatever choose_jacobian provided
      A_matrix_var <- jac_choice$A
      A_source <- jac_choice$source
      A_condition <- jac_choice$kappa
    }
  }
  list(
    A_matrix_var = A_matrix_var,
    A_source = A_source,
    A_condition = A_condition,
    A_diff_norm = A_diff_norm,
    jacobian_auto_rule = jacobian_auto_rule
  )
}

#' Delta variance computation wrapper
#'
#' Encapsulates Jacobian selection, stability heuristics (ridge/pinv), and
#' analytic vs numeric gradient selection for the mean functional, returning
#' standard errors and diagnostics in a single place.
#'
#' @noRd
#' @keywords internal
#' @return List with entries: `se_y_hat`, `vcov_unscaled`, `vcov_result`, and
#'   `diag` (a list of diagnostics: A_source/condition/diff/auto_rule, used_ridge,
#'   used_pseudoinverse, invert_rule, and variance_auto_rule).
el_variance_delta <- function(equation_system_func,
                              analytical_jac_func,
                              estimates,
                              variance_jacobian,
                              family,
                              response_model_matrix_scaled,
                              response_model_matrix_unscaled,
                              auxiliary_matrix_scaled,
                              mu_x_scaled,
                              eta_i_hat,
                              w_i_hat,
                              W_hat,
                              denominator_hat,
                              lambda_W_hat,
                              full_data,
                              compute_score_variance_func,
                              respondent_weights,
                              N_pop,
                              n_resp_weighted,
                              trim_cap,
                              outcome_vec,
                              K_beta,
                              standardize,
                              nmar_scaling_recipe,
                              variance_ridge,
                              variance_pseudoinverse) {
  sel <- el_select_variance_jacobian(
    equation_system_func = equation_system_func,
    analytical_jac_func = analytical_jac_func,
    estimates = estimates,
    variance_jacobian = variance_jacobian
  )
  A_matrix_var <- sel$A_matrix_var
  A_source <- sel$A_source
  A_condition <- sel$A_condition
  A_diff_norm <- sel$A_diff_norm
  jacobian_auto_rule <- sel$jacobian_auto_rule

  se_y_hat <- NA_real_
  vcov_unscaled <- matrix(NA_real_, K_beta, K_beta, dimnames = list(colnames(response_model_matrix_unscaled), colnames(response_model_matrix_unscaled)))
  auto_var_rule <- NA_character_

  vcov_try <- tryCatch({
    A_matrix <- if (!is.null(A_matrix_var)) A_matrix_var else if (!is.null(analytical_jac_func)) analytical_jac_func(estimates) else numDeriv::jacobian(func = equation_system_func, x = estimates)
# Heuristic variance stabilization: if user did not request ridge/pinv, enable when diagnostics indicate fragility
    var_ridge_eff <- variance_ridge
    var_pinv_eff <- variance_pseudoinverse
    if (identical(variance_ridge, FALSE) && !isTRUE(variance_pseudoinverse)) {
      kappa_ridge_thr <- getOption("nmar.var_auto_kappa_ridge_thr", 1e12)
      kappa_pinv_thr <- getOption("nmar.var_auto_kappa_pinv_thr", 1e14)
      if (is.finite(A_condition) && A_condition > kappa_ridge_thr) {
        var_ridge_eff <- TRUE
        auto_var_rule <- "ridge_due_to_kappa"
      }
      if (is.finite(A_diff_norm) && A_diff_norm > getOption("nmar.var_auto_rel_diff_thr", 1e-3) && !isTRUE(var_ridge_eff)) {
        var_ridge_eff <- TRUE
        auto_var_rule <- "ridge_due_to_rel_diff"
      }
      if (is.finite(A_condition) && A_condition > kappa_pinv_thr) {
        var_pinv_eff <- TRUE
        auto_var_rule <- "pinv_due_to_kappa"
      }
    }
    res <- el_compute_delta_variance(
      A_matrix = A_matrix,
      family = family,
      response_model_matrix_scaled = response_model_matrix_scaled,
      auxiliary_matrix_scaled = auxiliary_matrix_scaled,
      mu_x_scaled = mu_x_scaled,
      eta_i_hat = eta_i_hat,
      w_i_hat = w_i_hat,
      W_hat = W_hat,
      denominator_hat = denominator_hat,
      lambda_W_hat = lambda_W_hat,
      full_data = full_data,
      compute_score_variance_func = compute_score_variance_func,
      respondent_weights = respondent_weights,
      N_pop = N_pop,
      n_resp_weighted = n_resp_weighted,
      trim_cap = trim_cap,
      outcome_vec = outcome_vec,
      estimates = estimates,
      variance_ridge = var_ridge_eff,
      variance_pseudoinverse = var_pinv_eff
    )
    list(result = res, message = "Calculation successful", var_ridge_eff = var_ridge_eff, var_pinv_eff = var_pinv_eff)
  }, error = function(e) list(result = NULL, message = e$message))

  vcov_result <- vcov_try$result
  vcov_message <- vcov_try$message
  var_ridge_eff <- vcov_try$var_ridge_eff
  var_pinv_eff <- vcov_try$var_pinv_eff

# Fallback: if numeric source failed and analytic is available, retry with analytic A
  need_var_fb <- (is.null(vcov_result) || !is.finite(as.numeric(vcov_result$var_y_hat))) && (!is.null(analytical_jac_func)) && (variance_jacobian == "numeric" || (variance_jacobian == "auto" && A_source == "numeric"))
  if (need_var_fb) {
    vcov_try_fb <- tryCatch({
      A_matrix_alt <- analytical_jac_func(estimates)
      res <- el_compute_delta_variance(
        A_matrix = A_matrix_alt,
        family = family,
        response_model_matrix_scaled = response_model_matrix_scaled,
        auxiliary_matrix_scaled = auxiliary_matrix_scaled,
        mu_x_scaled = mu_x_scaled,
        eta_i_hat = eta_i_hat,
        w_i_hat = w_i_hat,
        W_hat = W_hat,
        denominator_hat = denominator_hat,
        lambda_W_hat = lambda_W_hat,
        full_data = full_data,
        compute_score_variance_func = compute_score_variance_func,
        respondent_weights = respondent_weights,
        N_pop = N_pop,
        n_resp_weighted = n_resp_weighted,
        trim_cap = trim_cap,
        outcome_vec = outcome_vec,
        estimates = estimates,
        variance_ridge = variance_ridge,
        variance_pseudoinverse = variance_pseudoinverse
      )
      list(result = res, message = "Calculation successful")
    }, error = function(e) list(result = vcov_result, message = vcov_message))
    vcov_result <- vcov_try_fb$result
    vcov_message <- vcov_try_fb$message
  }

  if (!is.null(vcov_result)) {
    var_y_hat_val <- as.numeric(vcov_result$var_y_hat)
    if (is.finite(var_y_hat_val)) se_y_hat <- sqrt(pmax(var_y_hat_val, 0))
    vcov_beta_scaled <- vcov_result$vcov_matrix_sandwich_scaled[1:K_beta, 1:K_beta, drop = FALSE]
    vcov_unscaled <- if (standardize) unscale_coefficients(estimates[1:K_beta], vcov_beta_scaled, nmar_scaling_recipe)$vcov else vcov_beta_scaled
    used_pseudoinverse <- isTRUE(vcov_result$used_pseudoinverse)
    used_ridge <- isTRUE(vcov_result$used_ridge)
    invert_rule <- if (!is.null(vcov_result$invert_rule)) vcov_result$invert_rule else NA_character_
  } else {
    used_pseudoinverse <- FALSE
    used_ridge <- FALSE
    invert_rule <- NA_character_
  }
  if (!is.finite(se_y_hat)) se_y_hat <- NA_real_

  list(
    se_y_hat = se_y_hat,
    vcov_unscaled = vcov_unscaled,
    vcov_result = vcov_result,
    vcov_message = vcov_message,
    diag = list(
      A_condition = A_condition,
      A_source = A_source,
      A_diff_norm = A_diff_norm,
      jacobian_auto_rule = jacobian_auto_rule,
      used_pseudoinverse = used_pseudoinverse,
      used_ridge = used_ridge,
      invert_rule = invert_rule,
      variance_auto_rule = auto_var_rule
    )
  )
}

#' Analytic gradient of the mean functional g(θ)
#'
#' Computes the analytic gradient of the respondent-weighted mean functional
#'   g = sum_i(pi_i y_i)/sum_i(pi_i) under smooth conditions (no trimming),
#'   for the EL reparameterization (β, z = logit(W), λ_x). This mirrors the
#'   guarded denominators used in post-solution weight construction.
#'
#' @keywords internal
el_grad_g_analytic <- function(family,
                               X_beta,
                               Xc,
                               mu_x_scaled,
                               respondent_weights,
                               eta_i_hat,
                               w_i_hat,
                               W_hat,
                               denominator_hat,
                               lambda_W_hat,
                               outcome_vec,
                               n_resp_weighted,
                               N_pop) {
  K_beta <- ncol(X_beta)
  K_aux <- if (is.null(Xc) || ncol(Xc) == 0) 0 else ncol(Xc)
  denom <- pmax(as.numeric(denominator_hat), 1e-8)
  pi <- respondent_weights / denom
  B_sum <- sum(pi)
  if (!is.finite(B_sum) || B_sum <= 0) return(NULL)
  g_val <- sum(pi * outcome_vec) / B_sum
  mu_eta <- family$mu.eta(eta_i_hat)
# d denom / d beta = lambda_W * mu_eta * X
  dden_deta <- lambda_W_hat * mu_eta
  w_resp <- respondent_weights
  factor <- (outcome_vec - g_val) * (-w_resp) / (denom^2)
  grad_beta <- as.numeric(colSums(X_beta * as.numeric(dden_deta * factor))) / B_sum
# Grad wrt z (W): d denom / d z = d lambda_W / d z * (w_i - W) - lambda_W * dW/dz
  C_const <- (N_pop / n_resp_weighted) - 1
  dW_dz <- W_hat * (1 - W_hat)
  dlam_dz <- C_const * W_hat / (1 - W_hat) # = lambda_W * W/(1-W)
  dden_dz <- dlam_dz * (w_i_hat - W_hat) - lambda_W_hat * dW_dz
  grad_z <- sum((outcome_vec - g_val) * (-w_resp * dden_dz / (denom^2))) / B_sum
# Grad wrt lambda_x: d denom / d lambda_x = Xc
  if (K_aux > 0) {
    grad_lambda <- as.numeric(colSums(Xc * as.numeric((outcome_vec - g_val) * (-w_resp) / (denom^2)))) / B_sum
  } else {
    grad_lambda <- numeric(0)
  }
  c(grad_beta, grad_z, grad_lambda)
}
