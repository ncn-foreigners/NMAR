#' EL auxiliary design resolution and population means
#'
#' Computes the respondent-side auxiliary matrix and the population means vector
#' used for centering (X - mu_x). Validates user-provided `auxiliary_means`
#' when supplied.
#'
#' @keywords internal
el_resolve_auxiliaries <- function(auxiliary_matrix_resp,
                                   auxiliary_matrix_full,
                                   auxiliary_means,
                                   weights_full = NULL) {
  n_resp <- if (!is.null(auxiliary_matrix_resp) && is.matrix(auxiliary_matrix_resp)) nrow(auxiliary_matrix_resp) else 0

  if (is.null(auxiliary_matrix_resp) || !is.matrix(auxiliary_matrix_resp) || ncol(auxiliary_matrix_resp) == 0) {
    return(list(
      matrix = matrix(nrow = n_resp, ncol = 0),
      means = NULL
    ))
  }

  aux_resp <- auxiliary_matrix_resp
  aux_full <- auxiliary_matrix_full

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
    stop(
      "Internal error: auxiliary full matrix is missing while auxiliaries are requested.",
      call. = FALSE
    )
  }
  el_validate_auxiliary_full(aux_full)

  if (!is.null(weights_full)) {
    mu <- as.numeric(colSums(aux_full * weights_full) / sum(weights_full))
  } else {
    mu <- as.numeric(colMeans(aux_full))
  }
  names(mu) <- colnames(aux_resp)
  list(matrix = aux_resp, means = mu)
}

#' Validate that full auxiliary data has no NA when means must be computed
#' @keywords internal
el_validate_auxiliary_full <- function(aux_full) {
  if (is.null(aux_full) || !is.matrix(aux_full) || ncol(aux_full) == 0) return(invisible(NULL))
  if (!anyNA(aux_full)) return(invisible(NULL))
  na_loc <- which(is.na(aux_full), arr.ind = TRUE)[1, , drop = TRUE]
  bad_col <- colnames(aux_full)[na_loc[2]]
  bad_row <- na_loc[1]
  msg <- sprintf(
    "Auxiliary variables contain NA values in full data (needed to compute population means).%s\nEither provide auxiliary_means=... or remove NA values from auxiliary variables.",
    if (is.finite(bad_row)) sprintf("\nFirst NA in '%s' at row %d", bad_col, bad_row) else ""
  )
  stop(msg, call. = FALSE)
}
