#' EL auxiliary design resolution and population means
#'
#' Computes the respondent-side auxiliary matrix and the population means vector
#' used for centering (X - mu_x). When `auxiliary_means` is supplied, only
#' respondent rows are required to be fully observed; NA values are permitted on
#' nonrespondent rows. When `auxiliary_means` is NULL, auxiliaries must be fully
#' observed in the full data used to estimate population means.
#'
#' @keywords internal
el_resolve_auxiliaries <- function(aux_design_full,
                                   respondent_mask,
                                   auxiliary_means,
                                   weights_full = NULL) {
  n_resp <- if (!is.null(respondent_mask)) sum(respondent_mask) else 0
  if (is.null(aux_design_full) || !is.matrix(aux_design_full) || ncol(aux_design_full) == 0) {
    return(list(
      auxiliary_design = matrix(nrow = n_resp, ncol = 0),
      means = NULL
    ))
  }

  aux_resp <- aux_design_full[respondent_mask, , drop = FALSE]

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
    return(list(auxiliary_design = aux_resp, means = mu))
  }

  el_validate_matrix(
    aux_design_full,
    allow_na = FALSE,
    label = "Auxiliary variables",
    severity = "error",
    row_map = NULL,
    scope_note = " in full data (needed to compute population means)",
    plural_label = TRUE
  )

  if (!is.null(weights_full)) {
    mu <- as.numeric(colSums(aux_design_full * weights_full) / sum(weights_full))
  } else {
    mu <- as.numeric(colMeans(aux_design_full))
  }
  names(mu) <- colnames(aux_resp)
  list(auxiliary_design = aux_resp, means = mu)
}

#' Strata augmentation for survey designs
#'
#' Augments the auxiliary design with strata dummies (dropping one level) and
#' appends stratum-share means when weights and N_pop are available. Intended for
#' survey workflows only.
#'
#' @keywords internal
el_augment_strata_aux <- function(aux_design_full,
                                  strata_factor,
                                  weights_full,
                                  N_pop,
                                  auxiliary_means) {
  if (is.null(aux_design_full) || !is.matrix(aux_design_full)) {
    return(list(mat = aux_design_full, means = auxiliary_means))
  }
  if (is.null(strata_factor) || length(unique(strata_factor)) <= 1L) {
    return(list(mat = aux_design_full, means = auxiliary_means))
  }
  if (nrow(aux_design_full) != length(strata_factor)) {
    stop("Strata factor must align with auxiliary design rows.", call. = FALSE)
  }
  strata_levels <- levels(strata_factor)
  if (length(strata_levels) <= 1L) {
    return(list(mat = aux_design_full, means = auxiliary_means))
  }
  dummy_levels <- strata_levels[-1L]
  strata_mat <- stats::model.matrix(~strata_factor)[, -1, drop = FALSE]
  colnames(strata_mat) <- paste0("strata_", dummy_levels)
  mat_aug <- cbind(aux_design_full, strata_mat)

  means_aug <- auxiliary_means
  if (!is.null(weights_full) && !is.null(N_pop) && is.finite(N_pop)) {
    W_h <- vapply(strata_levels, function(lev) {
      sum(weights_full[strata_factor == lev])
    }, numeric(1))
    W_h <- W_h / N_pop
    strata_means <- W_h[dummy_levels]
    names(strata_means) <- paste0("strata_", dummy_levels)
    means_aug <- c(means_aug, strata_means)
  }

  list(mat = mat_aug, means = means_aug)
}
