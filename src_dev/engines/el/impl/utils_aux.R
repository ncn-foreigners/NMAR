#' EL auxiliary inconsistency checks (engine-local)
#' @name el_utils_aux
#' @keywords internal
#' @noRd
NULL

#' Check auxiliary means consistency against respondents' sample support
#'
#' @param respondent_df data.frame of respondents (no NAs in outcome indicator)
#' @param aux_formula RHS-only formula for auxiliaries (no intercept)
#' @param provided_means optional named numeric vector of auxiliary means on the same columns as model matrix
#' @param threshold numeric, z-score threshold for flagging
#' @return list(max_z = numeric(1) or NA, cols = character())
#' @keywords internal
el_check_aux_inconsistency <- function(respondent_df, aux_formula, provided_means = NULL, threshold = 8) {
  out <- list(max_z = NA_real_, cols = character(0))
  if (is.null(aux_formula)) return(out)
  mm <- tryCatch(model.matrix(aux_formula, data = respondent_df), error = function(e) NULL)
  if (is.null(mm) || ncol(mm) == 0) return(out)
# Drop near-constant columns
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
