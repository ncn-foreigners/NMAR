#' Input preprocessing
#'
#' Parses the two-part Formula, constructs EL design matrices, injects the
#' respondent delta indicator, attaches weights and survey metadata,
#' and returns the pieces needed by the EL core.
#'
#' Enforeces the following format required by the rest of el code:
#' \itemize{
#' \item LHS references exactly one outcome source variable in \code{data}; any
#' transforms are applied via the formula environment and must be defined
#' for all respondent rows.
#' \item The outcome is never allowed to appear on RHS1 (auxiliaries) or RHS2
#' (missingness predictors), either explicitly in the formula or implicitly
#' via dot (\code{.}) expansion. The missingness model uses the evaluated LHS
#' expression as a dedicated predictor column instead.
#' \item RHS1 always yields an intercept-free auxiliary design matrix with
#' k-1 coding for factor auxiliaries, regardless of user \code{+0}/\code{-1} syntax or
#' custom contrasts. Auxiliary columns are validated to be fully observed
#' and non-constant among respondents.
#' \item RHS2 always yields a missingness-design matrix for respondents that
#' includes an intercept column and zero-variance predictors emit
#' warnings. NA among respondents is rejected.
#' \item \code{respondent_mask} is defined from the raw outcome in \code{data}, not from
#' the transformed LHS. An injected \code{..nmar_delta..} indicator in
#' \code{analysis_data} must match this mask.
#' \item \code{N_pop} is the analysis-scale population size:
#' for IID it is \code{nrow(data)} unless overridden by \code{n_total}.
#' For survey designs it is \code{sum(weights)} or \code{n_total} when supplied.
#' }
#'
#' @keywords internal
el_prepare_inputs <- function(formula,
                              data,
                              weights = NULL,
                              n_total = NULL,
                              design_object = NULL) {
  if (missing(formula)) stop("`formula` must be supplied.", call. = FALSE)

  fml <- Formula::as.Formula(formula)
  parts <- length(fml)
  n_rhs_parts <- if (length(parts) >= 2L) parts[2L] else 1L
  if (n_rhs_parts > 2L) {
    stop("EL formulas support at most two RHS partitions (auxiliaries | missingness predictors).", call. = FALSE)
  }

  base_formula <- tryCatch(
    stats::as.formula(fml),
    error = function(e) stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  )
  lhs_expr <- tryCatch(base_formula[[2L]], error = function(e) NULL)
  if (is.null(lhs_expr)) {
    stop("`formula` must be a two-sided formula, e.g., y ~ x1 + x2.", call. = FALSE)
  }
  outcome_label <- paste(deparse(lhs_expr), collapse = " ")
  lhs_vars <- unique(all.vars(lhs_expr))
  if (length(lhs_vars) == 0) {
    stop("Left-hand side must reference a variable in `data`.", call. = FALSE)
  }
  if (length(lhs_vars) > 1) {
    stop(
      "Left-hand side may reference only one outcome variable; pre-compute other transforms before modeling.",
      call. = FALSE
    )
  }
  outcome_source <- lhs_vars[[1L]]
  if (!outcome_source %in% names(data)) {
    stop(sprintf("Variables not found in data: %s", outcome_source), call. = FALSE)
  }

  model_frame <- stats::model.frame(fml, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE)
  raw_outcome <- data[[outcome_source]]
  response_vector <- Formula::model.part(fml, data = model_frame, lhs = 1, drop = TRUE)
  if (length(response_vector) != nrow(model_frame)) {
    stop("Internal error: response extraction did not align with data.", call. = FALSE)
  }

  transformed_na <- which(!is.na(raw_outcome) & is.na(response_vector))
  if (length(transformed_na) > 0) {
    stop(
      sprintf(
        "LHS expression '%s' produced NA/NaN for observed outcome rows. Ensure the transform is defined for all respondents.",
        outcome_label
      ),
      call. = FALSE
    )
  }

  respondent_mask <- unname(!is.na(raw_outcome))
  if (!any(respondent_mask)) {
    stop("No respondents detected in data after preprocessing.", call. = FALSE)
  }

  if (is.logical(response_vector)) {
    warning("Coercing logical outcome to numeric 0/1.", call. = FALSE)
    response_vector <- as.numeric(response_vector)
  } else if (is.factor(response_vector) || is.ordered(response_vector)) {
    levs <- levels(response_vector)
    if (length(levs) == 2L) {
      warning(sprintf("Coercing two-level factor outcome (%s) to numeric 0/1.", paste(levs, collapse = "/")), call. = FALSE)
      response_vector <- as.numeric(response_vector) - 1
    } else {
      stop("Outcome variable must be numeric after evaluating the left-hand side.", call. = FALSE)
    }
  } else if (!is.numeric(response_vector)) {
    stop("Outcome variable must be numeric after evaluating the left-hand side.", call. = FALSE)
  }

  outcome_names <- c(outcome_source, outcome_label)

  rhs1_formula <- stats::formula(fml, lhs = 0, rhs = 1)
  rhs1_expr <- rhs1_formula[[2L]]
  rhs1_vars_explicit <- setdiff(all.vars(rhs1_expr), ".")
  rhs1_terms <- stats::terms(rhs1_formula, data = model_frame)
  el_assert_no_offset(rhs1_terms, "auxiliary predictors")

# Forces k-1 factor coding while keeping an intercept-free auxiliary block
  rhs1_tlabs <- attr(rhs1_terms, "term.labels")
  rhs1_forced <- if (length(rhs1_tlabs)) {
    stats::reformulate(rhs1_tlabs, response = NULL, intercept = TRUE)
  } else {
    ~1
  }
  rhs1_mat_full <- stats::model.matrix(rhs1_forced, data = model_frame)
  keep_rhs1 <- !(colnames(rhs1_mat_full) %in% c("(Intercept)", outcome_names))
  aux_design_full <- rhs1_mat_full[, keep_rhs1, drop = FALSE]
  if (outcome_source %in% rhs1_vars_explicit) {
    stop("The outcome cannot appear in auxiliary constraints.", call. = FALSE)
  }

  if (ncol(aux_design_full) > 0) {
    el_validate_matrix(
      aux_design_full[respondent_mask, , drop = FALSE],
      allow_na = FALSE,
      label = "Auxiliary covariate",
      severity = "error",
      row_map = which(respondent_mask)
    )
  }

  if (n_rhs_parts >= 2L) {
    rhs2_formula <- stats::formula(fml, lhs = 0, rhs = 2)
    rhs2_expr <- rhs2_formula[[2L]]
    rhs2_terms <- stats::terms(rhs2_formula, data = model_frame)
    el_assert_no_offset(rhs2_terms, "missingness predictors")

    rhs2_vars_explicit <- setdiff(all.vars(rhs2_expr), ".")
    if (outcome_source %in% rhs2_vars_explicit || outcome_label %in% rhs2_vars_explicit) {
      stop(
        "Outcome cannot appear in missingness predictors; encode its effect via the left-hand side.",
        call. = FALSE
      )
    }

    rhs2_mat_full <- stats::model.matrix(fml, data = model_frame, rhs = 2)

    if (any(colnames(rhs2_mat_full) %in% outcome_names)) {
      stop(
        "Outcome cannot appear in missingness predictors; encode its effect via the left-hand side.",
        call. = FALSE
      )
    }
  } else {
    rhs2_mat_full <- matrix(nrow = nrow(model_frame), ncol = 0)
  }

  if (!"(Intercept)" %in% colnames(rhs2_mat_full)) {
    rhs2_mat_full <- cbind("(Intercept)" = 1, rhs2_mat_full)
  }

  rhs2_mat_resp <- rhs2_mat_full[respondent_mask, , drop = FALSE]
  rhs2_validate <- rhs2_mat_resp
  if ("(Intercept)" %in% colnames(rhs2_validate)) {
    rhs2_validate <- rhs2_validate[, colnames(rhs2_validate) != "(Intercept)", drop = FALSE]
  }

  if (ncol(rhs2_validate) > 0) {
    el_validate_matrix(
      rhs2_validate,
      allow_na = FALSE,
      label = "Missingness-model predictor",
      severity = "warn",
      row_map = which(respondent_mask)
    )
  }

  intercept_col <- rhs2_mat_resp[, "(Intercept)", drop = FALSE]
  rhs2_predictors <- rhs2_mat_resp[, colnames(rhs2_mat_resp) != "(Intercept)", drop = FALSE]
  outcome_col <- matrix(response_vector[respondent_mask], ncol = 1)
  colnames(outcome_col) <- outcome_label
  missingness_design <- cbind(intercept_col, outcome_col, rhs2_predictors)

  if (length(respondent_mask) != nrow(data)) {
    stop("Internal error: respondent mask must align with data.", call. = FALSE)
  }
  delta_name <- "..nmar_delta.."
  if (delta_name %in% names(data)) {
    i <- 1L
    while (paste0(delta_name, i) %in% names(data)) i <- i + 1L
    delta_name <- paste0(delta_name, i)
  }
  data_aug <- data
  data_aug[[delta_name]] <- as.integer(respondent_mask)

  respondent_weights <- if (!is.null(weights)) {
    if (length(weights) != nrow(data)) {
      stop("`weights` must align with the number of rows in the data.", call. = FALSE)
    }
    weights[respondent_mask]
  } else {
    rep(1, sum(respondent_mask))
  }

  is_survey <- !is.null(design_object) && inherits(design_object, "survey.design")
  analysis_object <- if (isTRUE(is_survey)) {
    design_object$variables <- data_aug
    design_object
  } else {
    data_aug
  }

  default_pop <- if (isTRUE(is_survey) && !is.null(weights)) sum(weights) else nrow(data_aug)
  N_pop_val <- n_total %||% default_pop

  design <- list(
    missingness_design = missingness_design,
    aux_design_full = aux_design_full,
    respondent_mask = respondent_mask,
    outcome_expr = outcome_label,
    analysis_data = analysis_object,
    respondent_weights = respondent_weights,
    N_pop = N_pop_val
  )

  el_validate_design_spec(design, data_nrow = nrow(data))
  design
}
