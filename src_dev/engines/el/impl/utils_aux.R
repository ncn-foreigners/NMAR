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

#' Resolve auxiliary design and population means
#'
#' Centralizes the construction of the auxiliary model matrix on respondents
#' and the corresponding population means used for centering (X - mu_x).
#'
#' Rules:
#' - If \code{aux_formula} is NULL or results in zero columns, auxiliaries are
#'   disabled (returns an empty matrix and NULL means).
#' - If \code{auxiliary_means} is provided, it is matched/reordered to the
#'   respondent design columns; unmatched names are dropped. If nothing matches,
#'   auxiliaries are disabled.
#' - Else, compute means from the full data:
#'     * survey.design: design-weighted column means using \code{weights(full_data)}
#'     * data.frame: simple column means
#'   Factor levels absent among respondents are dropped via intersection of
#'   model-matrix column names between full and respondent data.
#'
#' @param full_data data.frame or survey.design with all sampled units.
#' @param respondent_data data.frame of respondents only.
#' @param aux_formula RHS-only formula for auxiliaries (no intercept) or NULL.
#' @param auxiliary_means optional named numeric vector of population means in
#'   the auxiliary model-matrix columns.
#' @return list(matrix, means, has_aux) where `matrix` is the respondent-side
#'   auxiliary design on the unscaled space, `means` is a named numeric vector
#'   aligned to its columns (or NULL), and `has_aux` is a logical flag.
#' @keywords internal
el_resolve_auxiliaries <- function(full_data,
                                   respondent_data,
                                   aux_formula,
                                   auxiliary_means) {
# No auxiliaries requested
  if (is.null(aux_formula)) {
    return(list(matrix = matrix(nrow = nrow(respondent_data), ncol = 0),
                means = NULL,
                has_aux = FALSE))
  }

  aux_resp <- tryCatch(model.matrix(aux_formula, data = respondent_data),
                       error = function(e) NULL)
  if (is.null(aux_resp) || ncol(aux_resp) == 0) {
    return(list(matrix = matrix(nrow = nrow(respondent_data), ncol = 0),
                means = NULL,
                has_aux = FALSE))
  }

# User-supplied population means take precedence
  if (!is.null(auxiliary_means)) {
    keep <- intersect(colnames(aux_resp), names(auxiliary_means))
    if (length(keep) == 0L) {
      return(list(matrix = matrix(nrow = nrow(respondent_data), ncol = 0),
                  means = NULL,
                  has_aux = FALSE))
    }
    aux_resp <- aux_resp[, keep, drop = FALSE]
    mu <- as.numeric(auxiliary_means[keep])
    names(mu) <- keep
    return(list(matrix = aux_resp, means = mu, has_aux = TRUE))
  }

# Otherwise compute from the full data (design-weighted if survey)
  if (inherits(full_data, "survey.design")) {
    mm_full <- tryCatch(model.matrix(aux_formula, data = full_data$variables),
                        error = function(e) NULL)
    if (is.null(mm_full) || ncol(mm_full) == 0) {
      return(list(matrix = matrix(nrow = nrow(respondent_data), ncol = 0),
                  means = NULL,
                  has_aux = FALSE))
    }
    common_cols <- intersect(colnames(aux_resp), colnames(mm_full))
    if (length(common_cols) == 0L) {
      return(list(matrix = matrix(nrow = nrow(respondent_data), ncol = 0),
                  means = NULL,
                  has_aux = FALSE))
    }
    mm_full <- mm_full[, common_cols, drop = FALSE]
    w_full <- as.numeric(weights(full_data))
    mu <- as.numeric(colSums(mm_full * w_full) / sum(w_full))
    names(mu) <- colnames(mm_full)
    aux_resp <- aux_resp[, common_cols, drop = FALSE]
    mu <- mu[colnames(aux_resp)]
    return(list(matrix = aux_resp, means = mu, has_aux = TRUE))
  } else {
    mm_full <- tryCatch(model.matrix(aux_formula, data = full_data),
                        error = function(e) NULL)
    if (is.null(mm_full) || ncol(mm_full) == 0) {
      return(list(matrix = matrix(nrow = nrow(respondent_data), ncol = 0),
                  means = NULL,
                  has_aux = FALSE))
    }
    common_cols <- intersect(colnames(aux_resp), colnames(mm_full))
    if (length(common_cols) == 0L) {
      return(list(matrix = matrix(nrow = nrow(respondent_data), ncol = 0),
                  means = NULL,
                  has_aux = FALSE))
    }
    mm_full <- mm_full[, common_cols, drop = FALSE]
    mu <- colMeans(mm_full)
    mu <- as.numeric(mu[common_cols])
    names(mu) <- common_cols
    aux_resp <- aux_resp[, common_cols, drop = FALSE]
    mu <- mu[colnames(aux_resp)]
    return(list(matrix = aux_resp, means = mu, has_aux = TRUE))
  }
}
