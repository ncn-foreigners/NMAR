#' @keywords internal
new_nmar_result_exptilt <- function(estimate, se, coefficients, vcov, model,
                                    converged = TRUE, weights = NULL,
                                    variance_message = NA_character_) {
  outcome_name <- model$outcome_label %||% model$col_y %||% NA_character_
  n_total <- if (!is.null(model$x)) nrow(model$x) else NA_integer_
  n_resp <- if (!is.null(model$x) && !is.null(model$col_y) && model$col_y %in% colnames(model$x)) {
    sum(!is.na(model$x[, model$col_y]))
  } else {
    NA_integer_
  }

  diagnostics <- list(
    loss_value = model$loss_value,
    iterations = model$iterations,
    variance_method = model$variance_method,
    bootstrap_reps = model$bootstrap_reps,
    control = model$control,
    stopping_threshold = model$stopping_threshold,
    sampling_performed = model$sampling_performed %||% FALSE,
    sample_size = model$sample_size %||% NA_integer_,
    original_n_total = model$original_n_total %||% NA_integer_,
    original_n_resp = model$original_n_resp %||% NA_integer_,
    original_n_nonresp = model$original_n_nonresp %||% NA_integer_
  )

# The unified data-frame/survey path populates model$is_survey, we also
# honor designs that arrive directly via the survey method so result
# metadata reflects the original call
  is_survey <- isTRUE(model$is_survey) || inherits(model$design, "survey.design")
  sample <- list(
    n_total = n_total,
    n_respondents = n_resp,
    is_survey = is_survey,
    design = if (is_survey) model$design else NULL
  )

  inference <- list(
    variance_method = model$variance_method %||% NA_character_,
    df = NA_real_,
    message = variance_message,
    used_pseudoinverse = FALSE,
    used_ridge = FALSE
  )

  meta <- list(
    engine_name = "exponential_tilting",
    call = model$call %||% NULL,
    formula = model$user_formula %||% model$formula %||% NULL
  )

  coeffs_vec <- coefficients %||% NULL
  vcov_mat <- vcov %||% NULL

  result <- new_nmar_result(
    estimate = estimate,
    estimate_name = outcome_name,
    se = se,
    converged = converged,
    model = list(
      coefficients = coeffs_vec,
      vcov = vcov_mat
    ),
    weights_info = list(values = weights, trimmed_fraction = NA_real_),
    sample = sample,
    inference = inference,
    diagnostics = diagnostics,
    meta = meta,
    extra = list(
      bootstrap_reps = model$bootstrap_reps,
      variance_method = model$variance_method,
      loss_value = model$loss_value,
      iterations = model$iterations,
      raw = list(model = model) # allow post-hoc diagnostics (e.g. score checks)
    ),
    class = "nmar_result_exptilt"
  )

  result
}
