#' Empirical likelihood for data frames
#'
#' @keywords internal
el.data.frame <- function(data, formula,
                          auxiliary_means = NULL, standardize = TRUE,
                          trim_cap = Inf, control = list(),
                          on_failure = c("return", "error"),
                          variance_method = c("bootstrap", "none"),
                          bootstrap_reps = 500,
                          n_total = NULL, start = NULL, trace_level = 0,
                          family = logit_family(), ...) {
  cl <- match.call()
  on_failure <- match.arg(on_failure)
  if (is.null(variance_method)) variance_method <- "none"
  variance_method <- match.arg(variance_method)

  inputs <- el_prepare_inputs(
    formula = formula,
    data = data,
    weights = NULL,
    n_total = n_total,
    design_object = NULL
  )

  respondents_only <- isTRUE(all(inputs$respondent_mask))
  has_aux <- is.matrix(inputs$aux_design_full) && ncol(inputs$aux_design_full) > 0L
  if (respondents_only && is.null(n_total)) {
    stop("Respondents-only data detected (no NAs in outcome), but 'n_total' was not provided.", call. = FALSE)
  }
  if (respondents_only && has_aux && is.null(auxiliary_means)) {
    stop(
      "Respondents-only data with auxiliary constraints requires auxiliary_means. Provide population auxiliary means via auxiliary_means=.",
      call. = FALSE
    )
  }

  auxiliary_summary <- el_resolve_auxiliaries(
    aux_design_full = inputs$aux_design_full,
    respondent_mask = inputs$respondent_mask,
    auxiliary_means = auxiliary_means,
    weights_full = NULL
  )

  core_results <- el_estimator_core(
    missingness_design = inputs$missingness_design,
    aux_matrix = auxiliary_summary$auxiliary_design,
    aux_means = auxiliary_summary$means,
    respondent_weights = inputs$respondent_weights,
    analysis_data = inputs$analysis_data,
    outcome_expr = inputs$outcome_expr,
    N_pop = inputs$N_pop,
    formula = formula,
    standardize = standardize,
    trim_cap = trim_cap,
    control = control,
    on_failure = on_failure,
    family = family,
    variance_method = variance_method,
    bootstrap_reps = bootstrap_reps,
    start = start,
    trace_level = trace_level,
    auxiliary_means = auxiliary_means
  )

  if (is.list(core_results$diagnostics)) {
    core_results$diagnostics$auxiliary_means <- auxiliary_summary$means
    core_results$diagnostics$auxiliary_matrix <- auxiliary_summary$auxiliary_design
  }

  inputs$variance_method <- variance_method
  el_build_result(core_results, inputs, cl, formula)
}
