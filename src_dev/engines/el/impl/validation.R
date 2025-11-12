#' @keywords internal
el_validate_outcome <- function(data, outcome_var, require_na) {
  if (!outcome_var %in% names(data)) {
    stop(sprintf("Variables not found in data: %s", outcome_var), call. = FALSE)
  }
  if (!is.numeric(data[[outcome_var]])) {
    stop(sprintf("Outcome variable '%s' must be numeric.", outcome_var), call. = FALSE)
  }
  if (isTRUE(require_na) && !anyNA(data[[outcome_var]])) {
    stop(sprintf("Outcome variable '%s' must contain NA values to indicate nonresponse.", outcome_var), call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
el_validate_vars_present <- function(data, aux_vars, resp_vars) {
  missing_in_aux <- setdiff(aux_vars, names(data))
  if (length(missing_in_aux) > 0) {
    stop(sprintf("Variables not found in data: %s", paste(missing_in_aux, collapse = ", ")), call. = FALSE)
  }
  missing_in_resp <- setdiff(resp_vars, names(data))
  if (length(missing_in_resp) > 0) {
    stop(sprintf("Variables not found in data: %s", paste(missing_in_resp, collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
el_validate_aux_no_na <- function(data, aux_expr, respondent_mask, auxiliary_means, env = parent.frame()) {
  if (is.null(aux_expr) || length(all.vars(aux_expr)) == 0) return(invisible(NULL))
# Build one-sided, no-intercept formula (~ 0 + aux_expr) and set environment
  aux_check_fml <- as.formula(call("~", call("+", 0, aux_expr)))
  environment(aux_check_fml) <- env
  aux_mf <- tryCatch(
    model.frame(aux_check_fml, data = data, na.action = stats::na.pass),
    error = function(e) NULL
  )
  if (!is.null(aux_mf) && ncol(aux_mf) > 0) {
# Always forbid NA among respondents (aux used in constraints for respondents)
    if (!is.null(respondent_mask) && anyNA(respondent_mask)) {
      stop("Internal error: respondent mask contains NA.", call. = FALSE)
    }
    if (!is.null(respondent_mask) && length(respondent_mask) == nrow(aux_mf)) {
      aux_resp <- aux_mf[respondent_mask, , drop = FALSE]
      if (anyNA(aux_resp)) {
        na_by_col <- vapply(seq_len(ncol(aux_resp)), function(j) anyNA(aux_resp[[j]]), logical(1))
        bad_col <- colnames(aux_resp)[which(na_by_col)[1]]
        bad_row_local <- which(is.na(aux_resp[[bad_col]]))[1]
# Map back to global row index for clarity
        bad_row <- which(respondent_mask)[bad_row_local]
        msg <- sprintf("Covariate '%s' contains NA values among respondents.%s",
                       bad_col,
                       if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else "")
        stop(msg, call. = FALSE)
      }
    }
# If auxiliary means are not provided, we must compute μ_x from full data => forbid NA anywhere
    if (is.null(auxiliary_means) && anyNA(aux_mf)) {
      na_by_col <- vapply(seq_len(ncol(aux_mf)), function(j) anyNA(aux_mf[[j]]), logical(1))
      bad_col <- colnames(aux_mf)[which(na_by_col)[1]]
      bad_row <- which(is.na(aux_mf[[bad_col]]))[1]
      msg <- sprintf("Auxiliary variables contain NA values in full data (needed to compute population means).%s\n",
                     if (is.finite(bad_row)) sprintf("\nFirst NA in '%s' at row %d", bad_col, bad_row) else "")
      msg <- paste0(msg, "Either provide auxiliary_means=... or remove NA values from auxiliary variables.")
      stop(msg, call. = FALSE)
    }
  }
  invisible(NULL)
}
