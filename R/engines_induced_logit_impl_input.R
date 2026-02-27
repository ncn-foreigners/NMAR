#' Induced-logit input preprocessing (IID data.frame)
#'
#' Parses `formula` of the form `y_miss ~ mu_covariates | x1_covariates` and
#' performs strict validation consistent with the paper assumptions:
#' covariates are fully observed for all sampled units and only the outcome may
#' be missing (to indicate nonresponse).
#'
#' This function also precomputes and caches the outcome-model and missingness-
#' model design matrices (`model.matrix()`) in the returned spec, so downstream
#' fitting steps can reuse them without recomputing or risking inconsistencies.
#'
#' Version 1 restriction: LHS must be a single variable name (no transforms).
#'
#' @keywords internal
#' @noRd
induced_logit_prepare_inputs <- function(formula, data) {
  validator_assert_formula_two_sided(formula, name = "formula")

  if (!is.data.frame(data)) {
    stop("induced-logit IID path requires `data` to be a data.frame.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("Input dataset is empty (0 rows).", call. = FALSE)
  }

  fml <- Formula::as.Formula(formula)
  parts <- length(fml)
  n_rhs_parts <- if (length(parts) >= 2L) parts[2L] else 1L
  if (n_rhs_parts > 2L) {
    stop("induced-logit formulas support at most two RHS partitions (mu | x1).", call. = FALSE)
  }

  base_formula <- tryCatch(
    stats::as.formula(fml),
    error = function(e) stop("`formula` must be a two-sided formula, e.g., y_miss ~ x1 + x2.", call. = FALSE)
  )

  lhs_expr <- base_formula[[2L]]
  if (!is.symbol(lhs_expr)) {
    stop(
      "For induced-logit v1, the left-hand side must be a single variable name (no transforms).",
      call. = FALSE
    )
  }

  outcome_name <- as.character(lhs_expr)
  if (!outcome_name %in% names(data)) {
    stop(sprintf("Outcome variable `%s` not found in data.", outcome_name), call. = FALSE)
  }

  y_raw <- data[[outcome_name]]
  if (is.logical(y_raw)) {
    stop("Outcome must be numeric (logical outcomes are not supported in induced-logit v1).", call. = FALSE)
  }
  if (is.factor(y_raw) || is.ordered(y_raw)) {
    stop("Outcome must be numeric (factor outcomes are not supported in induced-logit v1).", call. = FALSE)
  }
  if (!is.numeric(y_raw)) {
    stop("Outcome must be numeric.", call. = FALSE)
  }

  respondent_mask <- !is.na(y_raw)
  if (!any(respondent_mask)) stop("No respondents detected (all outcomes are NA).", call. = FALSE)
  if (!any(!respondent_mask)) stop("No nonrespondents detected (no NA outcomes).", call. = FALSE)
  if (any(!is.finite(y_raw[respondent_mask]))) {
    stop("Observed (respondent) outcome values must be finite (no Inf/-Inf).", call. = FALSE)
  }

  rhs_mu <- stats::formula(fml, lhs = 0, rhs = 1)
  rhs_x1 <- if (n_rhs_parts >= 2L) stats::formula(fml, lhs = 0, rhs = 2) else ~1

  mu_vars <- all.vars(rhs_mu)
  x1_vars <- all.vars(rhs_x1)
  if ("." %in% mu_vars || "." %in% x1_vars) {
    stop("`.` expansion is not supported yet in induced-logit formulas. Specify covariates explicitly.", call. = FALSE)
  }
  if (outcome_name %in% mu_vars || outcome_name %in% x1_vars) {
    stop("Outcome cannot appear on the right-hand side for induced-logit.", call. = FALSE)
  }
  reserved <- c(IL_COL_R, IL_COL_MU_HAT)
  if (outcome_name %in% reserved) {
    stop(
      "Outcome name is reserved for internal induced-logit bookkeeping: ",
      outcome_name,
      call. = FALSE
    )
  }
  if (any(reserved %in% mu_vars) || any(reserved %in% x1_vars)) {
    stop(
      "The following names are reserved for internal induced-logit bookkeeping and cannot appear in the formula: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }

  model_frame <- stats::model.frame(fml, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE)

  mu_mat <- stats::model.matrix(rhs_mu, data = model_frame)
  x1_mat <- stats::model.matrix(rhs_x1, data = model_frame)

  if (!"(Intercept)" %in% colnames(x1_mat)) {
    stop(
      "Missingness-model (x1) must include an intercept.\n  ",
      "Remove `0` / `-1` from the missingness block, or omit the `|` block to use the default `| 1`.",
      call. = FALSE
    )
  }

  assert_finite_matrix <- function(mat, label) {
    if (!is.matrix(mat)) stop("Internal error: expected matrix for ", label, ".", call. = FALSE)
    if (anyNA(mat)) {
      bad_cols <- colnames(mat)[apply(is.na(mat), 2, any)]
      bad_cols <- bad_cols[!is.na(bad_cols)]
      if (length(bad_cols) == 0L) bad_cols <- colnames(mat)
      stop(sprintf(
        "%s covariates must be fully observed (no NA) for all sampled units.\n  Offending columns: %s",
        label, paste(bad_cols, collapse = ", ")
      ), call. = FALSE)
    }
    if (any(!is.finite(mat))) {
      bad_cols <- colnames(mat)[apply(!is.finite(mat), 2, any)]
      bad_cols <- bad_cols[!is.na(bad_cols)]
      if (length(bad_cols) == 0L) bad_cols <- colnames(mat)
      stop(sprintf(
        "%s covariates must be finite for all sampled units (no Inf/-Inf).\n  Offending columns: %s",
        label, paste(bad_cols, collapse = ", ")
      ), call. = FALSE)
    }
    invisible(TRUE)
  }

  assert_finite_matrix(mu_mat, label = "Outcome-model (mu)")
  assert_finite_matrix(x1_mat, label = "Missingness-model (x1)")

# Factor levels observed only among nonrespondents can make mu_hat undefined
# This fails early with a targeted error instead of propagating NA predictions
  mf_mu_full <- stats::model.frame(rhs_mu, data = model_frame, na.action = stats::na.pass, drop.unused.levels = TRUE)
  if (ncol(mf_mu_full) > 0) {
    mf_mu_resp <- mf_mu_full[respondent_mask, , drop = FALSE]
    factor_cols <- names(mf_mu_full)[vapply(mf_mu_full, function(x) is.factor(x) || is.ordered(x), logical(1))]
    if (length(factor_cols) > 0) {
      missing_levels <- lapply(factor_cols, function(nm) {
        full_levels <- levels(mf_mu_full[[nm]])
        present_resp <- unique(as.character(mf_mu_resp[[nm]]))
        present_resp <- present_resp[!is.na(present_resp)]
        setdiff(full_levels, present_resp)
      })
      bad <- which(vapply(missing_levels, length, integer(1)) > 0L)
      if (length(bad) > 0L) {
        pieces <- vapply(bad, function(i) {
          nm <- factor_cols[[i]]
          levs <- missing_levels[[i]]
          sprintf("%s: %s", nm, paste(levs, collapse = ", "))
        }, character(1))
        stop(
          "Outcome-model (mu) factor levels must be observed among respondents.\n  ",
          "The following levels appear only among nonrespondents:\n  ",
          paste(pieces, collapse = "\n  "),
          call. = FALSE
        )
      }
    }
  }

  mu_terms <- stats::terms(rhs_mu)
  x1_terms <- stats::terms(rhs_x1)

  list(
    outcome = outcome_name,
    mu_rhs = rhs_mu,
    x1_rhs = rhs_x1,
    mu_term_labels = attr(mu_terms, "term.labels"),
    x1_term_labels = attr(x1_terms, "term.labels"),
    mu_intercept = isTRUE(attr(mu_terms, "intercept") == 1L),
    x1_intercept = isTRUE(attr(x1_terms, "intercept") == 1L),
    respondent_mask = respondent_mask,
    r_vec = as.integer(respondent_mask),
    model_frame = model_frame,
    mu_mat_full = mu_mat,
    x1_mat_full = x1_mat
  )
}
