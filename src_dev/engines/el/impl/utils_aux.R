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
el_check_aux_inconsistency <- function(respondent_df,
                                       aux_formula,
                                       provided_means = NULL,
                                       threshold = 8,
                                       precomputed_resp = NULL) {
  out <- list(max_z = NA_real_, cols = character(0))
  if (is.null(aux_formula)) return(out)
  mm <- precomputed_resp %||% tryCatch(model.matrix(aux_formula, data = respondent_df), error = function(e) NULL)
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
                                   auxiliary_means,
                                   precomputed_resp = NULL,
                                   precomputed_full = NULL) {
# No auxiliaries requested
  if (is.null(aux_formula)) {
    return(list(matrix = matrix(nrow = nrow(respondent_data), ncol = 0),
                means = NULL,
                has_aux = FALSE))
  }

  aux_resp <- precomputed_resp %||%
    tryCatch(model.matrix(aux_formula, data = respondent_data),
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
    mm_full <- precomputed_full %||%
      tryCatch(model.matrix(aux_formula, data = full_data$variables),
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
    mm_full <- precomputed_full %||%
      tryCatch(model.matrix(aux_formula, data = full_data),
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

el_build_precomputed_design <- function(design_matrices,
                                        estimation_data,
                                        outcome_var,
                                        respondent_indices) {
  if (is.null(design_matrices)) {
    return(NULL)
  }
  n_resp <- length(respondent_indices)
  response_pred_full <- design_matrices$response
  aux_full <- design_matrices$auxiliary

  if (!n_resp) {
    response_pred <- matrix(nrow = 0, ncol = if (is.null(response_pred_full)) 0 else ncol(response_pred_full))
    aux_resp <- matrix(nrow = 0, ncol = if (is.null(aux_full)) 0 else ncol(aux_full))
    response_matrix <- cbind(`(Intercept)` = numeric(0), matrix(nrow = 0, ncol = 0))
  } else {
    response_pred <- if (is.null(response_pred_full)) {
      matrix(nrow = n_resp, ncol = 0)
    } else {
      response_pred_full[respondent_indices, , drop = FALSE]
    }
    outcome_vals <- estimation_data[respondent_indices, outcome_var, drop = TRUE]
    outcome_col <- matrix(outcome_vals, ncol = 1)
    colnames(outcome_col) <- outcome_var
    if (!outcome_var %in% colnames(response_pred)) {
      response_pred <- cbind(outcome_col, response_pred)
    } else {
      response_pred[, outcome_var] <- outcome_vals
    }
    response_matrix <- cbind(`(Intercept)` = rep(1, n_resp), response_pred)
    aux_resp <- if (is.null(aux_full)) {
      matrix(nrow = n_resp, ncol = 0)
    } else {
      aux_full[respondent_indices, , drop = FALSE]
    }
  }

# Invariants for response design (defensive, method-agnostic assumption):
# - Intercept column exists
# - Outcome column exists and equals the provided outcome values (for respondent rows)
  if (n_resp > 0) {
    if (!"(Intercept)" %in% colnames(response_matrix)) {
      stop("Response design matrix must include an '(Intercept)' column.", call. = FALSE)
    }
    if (!outcome_var %in% colnames(response_matrix)) {
      stop("Response design matrix must include the outcome column '", outcome_var, "'.", call. = FALSE)
    }
# Check equality of outcome values
    yy <- estimation_data[respondent_indices, outcome_var, drop = TRUE]
    if (length(yy) != nrow(response_matrix) || any(!isTRUE(all.equal(as.numeric(response_matrix[, outcome_var]), as.numeric(yy))))) {
      stop("Response design matrix outcome column is inconsistent with data.", call. = FALSE)
    }
  }

  list(
    response_matrix = response_matrix,
    auxiliary_resp = aux_resp,
    auxiliary_full = aux_full
  )
}
