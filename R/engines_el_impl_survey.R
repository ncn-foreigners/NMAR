#' Empirical likelihood estimator for survey designs
#'
#' @keywords internal
el.survey.design <- function(data, formula,
                             auxiliary_means = NULL, standardize = TRUE,
                             strata_augmentation = TRUE,
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
  design <- data
  survey_ctx <- el_get_design_context(design)
  weights_initial <- as.numeric(weights(design))
  el_validate_survey_weights(weights_initial, n_obs = nrow(design$variables))
  strata_factor <- el_extract_strata_factor(design)
  respondent_weights_full <- weights_initial

  inputs <- tryCatch(
    el_prepare_inputs(
      formula = formula,
      data = design$variables,
      weights = respondent_weights_full,
      n_total = n_total,
      design_object = design
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      msg2 <- paste0(msg, sprintf("\nSurvey design info: ids = %s, strata = %s", survey_ctx$ids, survey_ctx$strata))
      stop(msg2, call. = FALSE)
    }
  )

  n_resp_weighted <- sum(inputs$respondent_weights)
  if (!is.finite(n_resp_weighted) || n_resp_weighted <= 0) {
    stop("Respondent weights must sum to a positive number for EL estimation.", call. = FALSE)
  }
  if (!is.finite(inputs$N_pop) || inputs$N_pop <= 0) {
    stop("`n_total`/`N_pop` must be a single positive number.", call. = FALSE)
  }
  if (inputs$N_pop + 1e-8 < n_resp_weighted) {
    stop(
      sprintf(
        "`n_total`/`N_pop` must be >= sum(respondent weights). Got N_pop = %.6f, sum(d_i) = %.6f.",
        inputs$N_pop, n_resp_weighted
      ),
      call. = FALSE
    )
  }

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

  auxiliary_means_eff <- auxiliary_means
  if (isTRUE(strata_augmentation) && !is.null(strata_factor) && length(unique(strata_factor)) > 1L) {
    if (isTRUE(all(inputs$respondent_mask))) {
      warning(
        "Wu-style strata augmentation with respondents-only survey data uses respondent weights ",
        "to approximate stratum shares; this is generally not recommended. ",
        "Consider strata_augmentation = FALSE or encoding known stratum totals via auxiliary_means.",
        call. = FALSE
      )
    }
    aug <- el_augment_strata_aux(
      aux_design_full = inputs$aux_design_full,
      strata_factor = strata_factor %||% NULL,
      weights_full = weights_initial,
      N_pop = inputs$N_pop,
      auxiliary_means = auxiliary_means_eff
    )
    inputs$aux_design_full <- aug$mat
    auxiliary_means_eff <- aug$means
  }

  X_full <- inputs$aux_design_full
  if (is.matrix(X_full) && ncol(X_full) > 0L) {
    X_resp <- X_full[inputs$respondent_mask, , drop = FALSE]
    if (nrow(X_resp) > 0L) {
      qr_rank <- tryCatch(qr(X_resp)$rank, error = function(e) NA_integer_)
      if (is.finite(qr_rank) && qr_rank < ncol(X_resp)) {
        warning(
          "Auxiliary + strata design appears rank-deficient (rank = ", qr_rank,
          " < ", ncol(X_resp), "). ",
          "Redundant constraints (for example, auxiliaries that are functions of strata) ",
          "can cause instability in the empirical likelihood solver. ",
          "Consider removing redundant auxiliaries or disabling strata_augmentation.",
          call. = FALSE
        )
      }
    }
  }

  auxiliary_summary <- el_resolve_auxiliaries(
    aux_design_full = inputs$aux_design_full,
    respondent_mask = inputs$respondent_mask,
    auxiliary_means = auxiliary_means_eff,
    weights_full = weights_initial
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
    auxiliary_means = auxiliary_means_eff
  )

  if (is.list(core_results$diagnostics)) {
    core_results$diagnostics$auxiliary_means <- auxiliary_summary$means
    core_results$diagnostics$auxiliary_matrix <- auxiliary_summary$auxiliary_design
  }

  inputs$variance_method <- variance_method
  el_build_result(core_results, inputs, cl, formula)
}

el_get_design_context <- function(design) {
  ctx <- list(ids = "<unspecified>", strata = "<unspecified>")
  dc <- try(getCall(design), silent = TRUE)
  if (!inherits(dc, "try-error") && !is.null(dc)) {
    args <- as.list(dc)[-1]
    get_arg <- function(nm) if (!is.null(args[[nm]])) args[[nm]] else NULL
    ids_val <- get_arg("ids"); id_val <- get_arg("id")
    ids_obj <- if (!is.null(ids_val)) ids_val else id_val
    strata_obj <- get_arg("strata")
    deparse1 <- function(x) paste(deparse(x, width.cutoff = 200L), collapse = " ")
    if (!is.null(ids_obj)) ctx$ids <- deparse1(ids_obj)
    if (!is.null(strata_obj)) ctx$strata <- deparse1(strata_obj)
  }
  ctx
}

el_validate_survey_weights <- function(weights, n_obs) {
  if (!is.numeric(weights)) stop("Survey design weights must be numeric.", call. = FALSE)
  if (length(weights) != n_obs) stop("Survey design weights must align with the number of observations.", call. = FALSE)
  if (any(!is.finite(weights))) stop("Survey design weights must be finite (no NA/NaN/Inf).", call. = FALSE)
  min_w <- suppressWarnings(min(weights))
  if (is.finite(min_w) && min_w < 0) {
    stop(sprintf("Survey design weights must be nonnegative (min = %.6f).", min_w), call. = FALSE)
  }
  sum_w <- sum(weights)
  if (!is.finite(sum_w) || sum_w <= 0) {
    stop("Survey design weights must sum to a positive number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Extract strata factor
#'
#' Looks for strata already materialized in the \code{survey.design} object.
#' When unavailable, attempts to reconstruct strata from the original
#' \code{svydesign()} call. When multiple stratification variables are supplied,
#' their interaction is used.
#'
#' @keywords internal
el_extract_strata_factor <- function(design) {
  if (!inherits(design, "survey.design")) return(NULL)

  strata_df <- design$strata
  if (is.data.frame(strata_df) && nrow(strata_df) == nrow(design$variables) && ncol(strata_df) >= 1L) {
    norm_col <- function(x) {
      if (is.factor(x) || is.ordered(x)) return(droplevels(x))
      as.factor(x)
    }
    if (ncol(strata_df) == 1L) {
      return(norm_col(strata_df[[1L]]))
    }
    strata_df_norm <- as.data.frame(lapply(strata_df, norm_col))
    return(interaction(strata_df_norm, drop = TRUE))
  }

  dc <- try(getCall(design), silent = TRUE)
  if (inherits(dc, "try-error") || is.null(dc)) return(NULL)
  args <- as.list(dc)[-1L]
  get_arg <- function(nm) if (!is.null(args[[nm]])) args[[nm]] else NULL
  strata_expr <- get_arg("strata")
  if (is.null(strata_expr)) return(NULL)
  vars <- design$variables

  strata_formula <- if (inherits(strata_expr, "formula")) {
    strata_expr
  } else {
    tryCatch(stats::as.formula(strata_expr),
             error = function(e) return(NULL)
    )
  }
  if (is.null(strata_formula)) return(NULL)
  mf <- tryCatch(
    stats::model.frame(strata_formula, data = vars, na.action = stats::na.pass),
    error = function(e) NULL
  )
  if (is.null(mf) || nrow(mf) != nrow(vars)) return(NULL)
  if (ncol(mf) == 1L) {
    x <- mf[[1L]]
    if (is.factor(x) || is.ordered(x)) droplevels(x) else as.factor(x)
  } else {
    mf_norm <- as.data.frame(lapply(mf, function(x) {
      if (is.factor(x) || is.ordered(x)) droplevels(x) else as.factor(x)
    }))
    interaction(mf_norm, drop = TRUE)
  }
}
