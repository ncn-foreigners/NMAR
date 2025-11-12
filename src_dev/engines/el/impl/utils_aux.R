#' EL auxiliary inconsistency checks (engine-local)
#' @name el_utils_aux
#' @keywords internal
#' @noRd
NULL

#' Check auxiliary means consistency against respondents' sample support
#'
#' @param respondent_df data.frame of respondents (no NAs in outcome indicator)
#' @param aux_formula RHS-only formula for auxiliaries (no intercept)
#' @param provided_means optional named numeric vector of auxiliary means on the same columns as the matrix
#' @param threshold numeric, z-score threshold for flagging
#' @return list(max_z = numeric(1) or NA, cols = character())
#' @keywords internal
el_check_aux_inconsistency_matrix <- function(aux_matrix_resp, provided_means = NULL, threshold = 8) {
  out <- list(max_z = NA_real_, cols = character(0))
  if (is.null(aux_matrix_resp) || !is.matrix(aux_matrix_resp) || ncol(aux_matrix_resp) == 0) return(out)

  mm <- aux_matrix_resp
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
#' Centralizes auxiliary matrix handling once the canonical Formula pipeline has
#' produced respondent/full model matrices. Computes the mean vector used for
#' centering (X - mu_x) and validates any user-provided population means.
#'
#' Rules:
#' - If \code{auxiliary_matrix_resp} has zero columns, auxiliaries are disabled.
#' - If \code{auxiliary_means} is provided, it must supply entries for every
#'   auxiliary design column (extras are ignored with a warning).
#' - Else, compute means from \code{auxiliary_matrix_full}; if \code{weights_full}
#'   is supplied (survey designs), use the weighted column means, otherwise use
#'   simple averages. Factor levels absent among respondents are automatically
#'   dropped by the shared Formula pipeline, so callers must pass aligned
#'   matrices.
#'
#' @param auxiliary_matrix_resp Respondent-side auxiliary design matrix.
#' @param auxiliary_matrix_full Full-sample auxiliary design matrix used for
#'   computing population means.
#' @param auxiliary_means optional named numeric vector of population means in
#'   the auxiliary model-matrix columns.
#' @param weights_full Optional vector of design weights aligned with rows of
#'   `auxiliary_matrix_full`.
#' @return list(matrix, means) where `matrix` is the respondent-side
#'   auxiliary design on the unscaled space and `means` is a named numeric vector
#'   aligned to its columns (or NULL).
#' @keywords internal
el_resolve_auxiliaries <- function(auxiliary_matrix_resp,
                                   auxiliary_matrix_full,
                                   auxiliary_means,
                                   weights_full = NULL) {
  n_resp <- if (!is.null(auxiliary_matrix_resp) && is.matrix(auxiliary_matrix_resp)) nrow(auxiliary_matrix_resp) else 0

# No auxiliaries requested
  if (is.null(auxiliary_matrix_resp) || !is.matrix(auxiliary_matrix_resp) || ncol(auxiliary_matrix_resp) == 0) {
    return(list(
      matrix = matrix(nrow = n_resp, ncol = 0),
      means = NULL
    ))
  }

  aux_resp <- auxiliary_matrix_resp
  aux_full <- auxiliary_matrix_full

  if ("(Intercept)" %in% colnames(aux_resp)) {
    aux_resp <- aux_resp[, colnames(aux_resp) != "(Intercept)", drop = FALSE]
  }
  if (!is.null(aux_full) && "(Intercept)" %in% colnames(aux_full)) {
    aux_full <- aux_full[, colnames(aux_full) != "(Intercept)", drop = FALSE]
  }

  if (ncol(aux_resp) == 0) {
    return(list(
      matrix = matrix(nrow = n_resp, ncol = 0),
      means = NULL
    ))
  }

  if (!is.null(auxiliary_means)) {
    provided_names <- names(auxiliary_means)
    missing_names <- setdiff(colnames(aux_resp), provided_names)
    if (length(missing_names) > 0) {
      stop(
        sprintf(
          "auxiliary_means must supply entries for: %s",
          paste(missing_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    extra <- setdiff(provided_names, colnames(aux_resp))
    if (length(extra) > 0) {
      warning(
        sprintf(
          "Ignoring unused names in 'auxiliary_means': %s",
          paste(extra, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    ordered_means <- auxiliary_means[colnames(aux_resp)]
    mu <- as.numeric(ordered_means)
    names(mu) <- colnames(aux_resp)
    return(list(matrix = aux_resp, means = mu))
  }

  if (is.null(aux_full) || ncol(aux_full) == 0) {
    return(list(
      matrix = matrix(nrow = nrow(aux_resp), ncol = 0),
      means = NULL
    ))
  }

  if (!identical(colnames(aux_resp), colnames(aux_full))) {
    common_cols <- intersect(colnames(aux_resp), colnames(aux_full))
    if (length(common_cols) == 0L) {
      return(list(
        matrix = matrix(nrow = nrow(aux_resp), ncol = 0),
        means = NULL
      ))
    }
    aux_resp <- aux_resp[, common_cols, drop = FALSE]
    aux_full <- aux_full[, common_cols, drop = FALSE]
  }

  if (!is.null(weights_full)) {
    mu <- as.numeric(colSums(aux_full * weights_full) / sum(weights_full))
  } else {
    mu <- as.numeric(colMeans(aux_full))
  }
  names(mu) <- colnames(aux_resp)
  list(matrix = aux_resp, means = mu)
}
