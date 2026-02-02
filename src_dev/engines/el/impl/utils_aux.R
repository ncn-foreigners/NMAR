#' Check auxiliary means consistency against respondents sample support.
#'
#' Computes a simple z-score diagnostic comparing user-supplied auxiliary means
#' to the respondents' sample means.
#'
#' @param auxiliary_matrix_resp Respondent-side auxiliary design matrix.
#' @param provided_means Optional named numeric vector of auxiliary means aligned to the matrix columns.
#' @return list(max_z = numeric(1) or NA, cols = character())
#'
#' @keywords internal
el_check_auxiliary_inconsistency_matrix <- function(auxiliary_matrix_resp, provided_means = NULL) {
  out <- list(max_z = NA_real_, cols = character(0))
  if (is.null(auxiliary_matrix_resp) || !is.matrix(auxiliary_matrix_resp) || ncol(auxiliary_matrix_resp) == 0) return(out)

  mm <- auxiliary_matrix_resp
# Remove nearly constant columns to avoid undefined z-scores
  sds <- apply(mm, 2, stats::sd)
  keep <- which(is.finite(sds) & sds >= 1e-8)
  if (length(keep) == 0) return(out)
  mm <- mm[, keep, drop = FALSE]
  sample_means <- colMeans(mm)
  sample_sds <- apply(mm, 2, stats::sd)
  out$cols <- colnames(mm)
  if (!is.null(provided_means)) {
    pm <- provided_means[out$cols]
    z <- abs((pm - sample_means) / pmax(sample_sds, 1e-8))
    out$max_z <- suppressWarnings(max(z, na.rm = TRUE))
  } else {
    out$max_z <- NA_real_
  }
  out
}

#' Assert that terms object lacks offsets
#' @keywords internal
el_assert_no_offset <- function(terms_obj, label) {
  if (is.null(terms_obj)) return(invisible(NULL))
  offsets <- attr(terms_obj, "offset")
  if (!is.null(offsets) && length(offsets) > 0) {
    stop(
      sprintf(
        "Offsets (offset()) are not supported for %s. Remove offset() from the formula.",
        label
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Validate matrix columns for NA and zero variance
#' @keywords internal
el_validate_matrix <- function(mat,
                               allow_na,
                               label,
                               severity,
                               row_map = NULL,
                               scope_note = NULL,
                               plural_label = FALSE) {
  severity <- match.arg(severity, c("error", "warn"))
  if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0 || nrow(mat) == 0) return(invisible(NULL))

  if (!allow_na && anyNA(mat)) {
    na_loc <- which(is.na(mat), arr.ind = TRUE)[1, , drop = TRUE]
    bad_col <- colnames(mat)[na_loc[2]]
    row_label <- if (!is.null(row_map) && length(row_map) >= na_loc[1]) row_map[na_loc[1]] else na_loc[1]
    note <- scope_note %||% if (!is.null(row_map)) " among respondents" else ""
    msg <- if (isTRUE(plural_label)) {
      sprintf(
        "%s contain NA values%s.%s",
        label,
        note,
        if (isTRUE(is.finite(row_label))) sprintf("\nFirst NA in '%s' at row %d", bad_col, row_label) else ""
      )
    } else {
      sprintf(
        "%s '%s' contains NA values%s.%s",
        label,
        bad_col,
        note,
        if (isTRUE(is.finite(row_label))) sprintf("\nFirst NA at row %d", row_label) else ""
      )
    }
    stop(msg, call. = FALSE)
  }

  is_constant <- vapply(seq_len(ncol(mat)), function(j) {
    col <- mat[, j]
    if (all(is.na(col))) return(FALSE)
    all(col == col[which(!is.na(col))[1]], na.rm = TRUE)
  }, logical(1))

  if (!any(is_constant)) return(invisible(NULL))
  cols <- colnames(mat)[is_constant]
  msg <- sprintf(
    "%s %s zero variance among respondents: %s",
    label,
    if (length(cols) == 1) "has" else "have",
    paste(cols, collapse = ", ")
  )
  if (identical(severity, "warn")) {
    warning(msg, call. = FALSE)
    return(invisible(NULL))
  }
  stop(msg, call. = FALSE)
}

#' Validate design dimensions
#' @keywords internal
el_validate_design_spec <- function(design, data_nrow) {
  mask <- design$respondent_mask
  if (!is.logical(mask)) {
    stop("Internal error: respondent mask must be logical.", call. = FALSE)
  }
  if (length(mask) != data_nrow) {
    stop(sprintf("Internal error: respondent mask length (%d) must equal data rows (%d).", length(mask), data_nrow), call. = FALSE)
  }
  missingness_design <- design$missingness_design
  if (!is.null(missingness_design) && nrow(missingness_design) != sum(mask)) {
    stop("Internal error: missingness design rows must equal respondent count.", call. = FALSE)
  }
  if (!is.null(design$aux_design_full) && nrow(design$aux_design_full) != data_nrow) {
    stop("Internal error: aux_design_full must align with original data rows.", call. = FALSE)
  }
  invisible(design)
}
