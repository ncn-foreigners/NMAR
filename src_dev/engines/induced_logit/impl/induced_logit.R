#' Induced logistic regression estimator
#'
#' Internal entry point used by `run_engine.nmar_engine_induced_logit()`.
#' S3 generic for `data.frame` and `survey.design` inputs.
#'
#' @keywords internal
#' @noRd
induced_logit <- function(data, formula,
                          variance_method = c("none", "bootstrap"),
                          bootstrap_reps = 500,
                          standardize = FALSE,
                          control = list(),
                          on_failure = c("return", "error"),
                          survey_design_policy = c("strict", "warn"),
                          keep_fits = FALSE,
                          ...) {
  UseMethod("induced_logit", data)
}

#' @keywords internal
#' @noRd
#' @exportS3Method induced_logit data.frame
induced_logit.data.frame <- function(data, formula,
                                     variance_method = c("none", "bootstrap"),
                                     bootstrap_reps = 500,
                                     standardize = FALSE,
                                     control = list(),
                                     on_failure = c("return", "error"),
                                     survey_design_policy = c("strict", "warn"),
                                     keep_fits = FALSE,
                                     ...) {
  cl <- match.call()
  variance_method <- match.arg(variance_method, choices = c("none", "bootstrap"))
  on_failure <- match.arg(on_failure, choices = c("return", "error"))
  survey_design_policy <- match.arg(survey_design_policy, choices = c("strict", "warn"))
  opts <- il_validate_induced_logit_opts(
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    standardize = standardize,
    control = control,
    on_failure = on_failure,
    survey_design_policy = survey_design_policy,
    keep_fits = keep_fits
  )

  tryCatch(
    {
      spec <- induced_logit_prepare_inputs(formula = formula, data = data)
      fit <- il_fit_from_backend(
        spec = spec,
        backend = il_backend_iid(data),
        standardize = opts$standardize,
        control = opts$control
      )
      tau_hat <- as.numeric(fit$point$tau_hat)
      if (!is.finite(tau_hat)) stop("Primary estimate is not finite.", call. = FALSE)

      variance_out <- il_compute_variance(
        data = data,
        formula = formula,
        point_estimate = tau_hat,
        variance_method = opts$variance_method,
        bootstrap_reps = opts$bootstrap_reps,
        standardize = opts$standardize,
        control = opts$control,
        survey_design_policy = opts$survey_design_policy,
        resample_guard = function(indices, data) {
          any(spec$respondent_mask[indices]) && any(!spec$respondent_mask[indices])
        }
      )

      il_build_success_result_iid(
        spec = spec,
        fit = fit,
        variance_out = variance_out,
        formula = formula,
        call = cl,
        variance_method = opts$variance_method,
        survey_design_policy = opts$survey_design_policy,
        keep_fits = opts$keep_fits
      )
    },
    error = function(e) {
      if (identical(opts$on_failure, "error")) stop(e$message, call. = FALSE)
      il_build_failure_result(
        message = e$message,
        formula = formula,
        call = cl,
        variance_method = opts$variance_method,
        survey_design_policy = opts$survey_design_policy,
        n_total = nrow(data),
        is_survey = FALSE,
        data_for_respondents = data,
        design = NULL,
        df = NA_real_
      )
    }
  )
}

#' @keywords internal
#' @noRd
#' @exportS3Method induced_logit survey.design
induced_logit.survey.design <- function(data, formula,
                                        variance_method = c("none", "bootstrap"),
                                        bootstrap_reps = 500,
                                        standardize = FALSE,
                                        control = list(),
                                        on_failure = c("return", "error"),
                                        survey_design_policy = c("strict", "warn"),
                                        keep_fits = FALSE,
                                        ...) {
  cl <- match.call()
  variance_method <- match.arg(variance_method, choices = c("none", "bootstrap"))
  on_failure <- match.arg(on_failure, choices = c("return", "error"))
  survey_design_policy <- match.arg(survey_design_policy, choices = c("strict", "warn"))
  opts <- il_validate_induced_logit_opts(
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    standardize = standardize,
    control = control,
    on_failure = on_failure,
    survey_design_policy = survey_design_policy,
    keep_fits = keep_fits
  )

  design <- data
  il_validate_survey_entry(design)

  tryCatch(
    {
      survey_assumptions <- il_validate_supported_survey_design(
        design = design,
        survey_design_policy = opts$survey_design_policy,
        context = "induced-logit survey path"
      )

      spec <- induced_logit_prepare_inputs(formula = formula, data = design$variables)
      backend <- il_backend_survey(design)
      fit <- il_fit_from_backend(spec = spec, backend = backend, standardize = opts$standardize, control = opts$control)
      tau_hat <- as.numeric(fit$point$tau_hat)
      if (!is.finite(tau_hat)) stop("Primary estimate is not finite.", call. = FALSE)

      variance_out <- il_compute_variance(
        data = design,
        formula = formula,
        point_estimate = tau_hat,
        variance_method = opts$variance_method,
        bootstrap_reps = opts$bootstrap_reps,
        standardize = opts$standardize,
        control = opts$control,
        survey_design_policy = opts$survey_design_policy
      )

      il_build_success_result_survey(
        spec = spec,
        fit = fit,
        design = design,
        backend = backend,
        survey_assumptions = survey_assumptions,
        variance_out = variance_out,
        formula = formula,
        call = cl,
        variance_method = opts$variance_method,
        survey_design_policy = opts$survey_design_policy,
        keep_fits = opts$keep_fits
      )
    },
    error = function(e) {
      if (identical(opts$on_failure, "error")) stop(e$message, call. = FALSE)
      il_build_failure_result(
        message = e$message,
        formula = formula,
        call = cl,
        variance_method = opts$variance_method,
        survey_design_policy = opts$survey_design_policy,
        n_total = il_survey_weight_total(design),
        is_survey = TRUE,
        data_for_respondents = design$variables,
        design = design,
        df = il_survey_df(design)
      )
    }
  )
}

#' @keywords internal
#' @noRd
il_validate_induced_logit_opts <- function(variance_method,
                                           bootstrap_reps,
                                           standardize,
                                           control,
                                           on_failure,
                                           survey_design_policy,
                                           keep_fits) {
  validator_assert_choice(variance_method, choices = c("none", "bootstrap"), name = "variance_method")
  validator_assert_choice(on_failure, choices = c("return", "error"), name = "on_failure")
  validator_assert_choice(survey_design_policy, choices = c("strict", "warn"), name = "survey_design_policy")

  validator_assert_positive_integer(bootstrap_reps, name = "bootstrap_reps", is.finite = TRUE)
  if (identical(variance_method, "bootstrap") && bootstrap_reps < 2) {
    stop("`bootstrap_reps` must be at least 2 when variance_method='bootstrap'.", call. = FALSE)
  }

  validator_assert_scalar_logical(standardize, name = "standardize")
  validator_assert_scalar_logical(keep_fits, name = "keep_fits")
  validator_assert_list(control, name = "control")

  list(
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    standardize = standardize,
    control = control,
    on_failure = on_failure,
    survey_design_policy = survey_design_policy,
    keep_fits = keep_fits
  )
}

#' @keywords internal
#' @noRd
#' @exportS3Method induced_logit default
induced_logit.default <- function(data, formula,
                                  variance_method = c("none", "bootstrap"),
                                  bootstrap_reps = 500,
                                  standardize = FALSE,
                                  control = list(),
                                  on_failure = c("return", "error"),
                                  survey_design_policy = c("strict", "warn"),
                                  keep_fits = FALSE,
                                  ...) {
  stop("induced-logit currently supports `data.frame` and `survey.design` inputs only.", call. = FALSE)
}
