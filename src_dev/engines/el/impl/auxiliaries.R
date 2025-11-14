#' EL auxiliary design resolution and population means
#'
#' Computes the respondent-side auxiliary matrix and the population means vector
#' used for centering (X - mu_x). When `auxiliary_means` is supplied, only
#' respondent rows are required to be fully observed; NA values are permitted on
#' nonrespondent rows. When `auxiliary_means` is NULL, auxiliaries must be fully
#' observed in the full data used to estimate population means.
#'
#' @keywords internal
el_resolve_auxiliaries <- function(auxiliary_design_full,
                                   respondent_mask,
                                   auxiliary_means,
                                   weights_full = NULL) {
  n_resp <- if (!is.null(respondent_mask)) sum(respondent_mask) else 0
  if (is.null(auxiliary_design_full) || !is.matrix(auxiliary_design_full) || ncol(auxiliary_design_full) == 0) {
    return(list(
      auxiliary_design = matrix(nrow = n_resp, ncol = 0),
      means = NULL
    ))
  }

  aux_resp <- auxiliary_design_full[respondent_mask, , drop = FALSE]

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

  el_assert_no_na(
    auxiliary_design_full,
    row_map = seq_len(nrow(auxiliary_design_full)),
    label = "Auxiliary variables",
    scope_note = " in full data (needed to compute population means)",
    plural_label = TRUE
  )

  if (!is.null(weights_full)) {
    mu <- as.numeric(colSums(auxiliary_design_full * weights_full) / sum(weights_full))
  } else {
    mu <- as.numeric(colMeans(auxiliary_design_full))
  }
  names(mu) <- colnames(aux_resp)
  list(auxiliary_design = aux_resp, means = mu)
}
