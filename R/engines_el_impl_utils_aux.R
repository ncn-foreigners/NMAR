#' EL auxiliary inconsistency checks (engine-local)
#' @name el_utils_aux
#' @keywords internal
#' @noRd
NULL

#' Check auxiliary means consistency against respondents' sample support
#'
#' @param aux_matrix_resp Respondent-side auxiliary design matrix.
#' @param provided_means Optional named numeric vector of auxiliary means aligned to the matrix columns.
#' @param threshold numeric, z-score threshold for flagging
#' @return list(max_z = numeric(1) or NA, cols = character())
#' @keywords internal
el_check_aux_inconsistency_matrix <- function(aux_matrix_resp, provided_means = NULL, threshold = 8) {
  out <- list(max_z = NA_real_, cols = character(0))
  if (is.null(aux_matrix_resp) || !is.matrix(aux_matrix_resp) || ncol(aux_matrix_resp) == 0) return(out)

  mm <- aux_matrix_resp
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

## el_resolve_auxiliaries moved to impl/auxiliaries.R to co-locate with validation
