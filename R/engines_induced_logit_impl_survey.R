#' Induced logistic regression estimator for survey designs
#'
#' @keywords internal
#' @noRd
il_validate_supported_survey_design <- function(design, survey_design_policy = c("strict", "warn"),
                                                context = "induced-logit survey path") {
  survey_design_policy <- match.arg(survey_design_policy)
  info <- nmar_detect_survey_design_assumptions(design)

  hard_fail <- info$has_poststrata || info$has_calibration
  if (isTRUE(hard_fail)) {
    stop(
      "Induced-logit survey path does not support calibrated/post-stratified designs.\n  ",
      "Detected: ", paste(info$risk_labels[info$risk_labels %in% c("post-stratification", "calibration")], collapse = ", "), ".\n  ",
      "Use an unadjusted survey.design for this estimator.",
      call. = FALSE
    )
  }

  if (isTRUE(info$has_probability_risk)) {
    msg <- paste0(
      context, ": induced-logit survey implementation assumes weight-based estimating equations.\n  ",
      "Detected design features: ", paste(info$risk_labels, collapse = ", "), ".\n  ",
      "These features can require additional probability/FPC handling beyond analysis weights."
    )
    if (identical(survey_design_policy, "strict")) {
      stop(msg, call. = FALSE)
    }
    warning(msg, call. = FALSE, immediate. = TRUE)
  }

  list(
    policy = survey_design_policy,
    has_probability_risk = isTRUE(info$has_probability_risk),
    risk_labels = info$risk_labels,
    flags = info$flags
  )
}

#' @keywords internal
#' @noRd
il_validate_survey_entry <- function(design) {
  il_require_survey()
  if (!inherits(design, "survey.design")) {
    stop("induced-logit survey path requires a survey.design object.", call. = FALSE)
  }
  if (!is.data.frame(design$variables)) {
    stop("Internal error: survey.design does not contain a data.frame in $variables.", call. = FALSE)
  }
  if (nrow(design$variables) == 0L) {
    stop("Input dataset is empty (0 rows).", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
il_survey_df <- function(design) {
  tryCatch(survey::degf(design), error = function(e) NA_real_)
}

#' @keywords internal
#' @noRd
il_survey_weight_total <- function(design) {
  tryCatch(sum(as.numeric(stats::weights(design))), error = function(e) NA_real_)
}
