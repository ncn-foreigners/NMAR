#' EL preparation utilities shared across IID and survey entry points
#'
#' Helper functions used while building design matrices: dropping intercept
#' columns, detecting explicit `+ 1`, validating respondent/full matrices, adding
#' the NMAR delta indicator, and translating `model.frame()` errors into user
#' messages.
#'
#' @keywords internal

el_drop_intercept_col <- function(mm) {
  if (is.null(mm) || !is.matrix(mm) || ncol(mm) == 0) return(mm)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  mm
}

el_formula_has_explicit_intercept <- function(node) {
  if (is.null(node)) return(FALSE)
  if (is.numeric(node) && length(node) == 1L && isTRUE(all.equal(node, 1))) return(TRUE)
  if (is.call(node)) {
    op <- as.character(node[[1L]])
    if (op %in% c("+", "-", "~")) {
      for (i in 2L:length(node)) {
        if (Recall(node[[i]])) return(TRUE)
      }
    }
  }
  FALSE
}

el_validate_response_matrix <- function(response_matrix, respondent_indices) {
  if (!is.matrix(response_matrix) || ncol(response_matrix) == 0 || nrow(response_matrix) == 0) {
    return(invisible(NULL))
  }
  if (!anyNA(response_matrix)) return(invisible(NULL))

  na_loc <- which(is.na(response_matrix), arr.ind = TRUE)[1, , drop = TRUE]
  bad_col <- colnames(response_matrix)[na_loc[2]]
  bad_row <- if (length(respondent_indices) >= na_loc[1]) respondent_indices[na_loc[1]] else NA_integer_

  msg <- sprintf(
    "Response-model predictor '%s' contains NA values among respondents.%s",
    bad_col,
    if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else ""
  )
  stop(msg, call. = FALSE)
}

el_validate_auxiliary_resp <- function(aux_resp, respondent_indices) {
  if (is.null(aux_resp) || !is.matrix(aux_resp) || ncol(aux_resp) == 0) return(invisible(NULL))

  if (anyNA(aux_resp)) {
    na_loc <- which(is.na(aux_resp), arr.ind = TRUE)[1, , drop = TRUE]
    bad_col <- colnames(aux_resp)[na_loc[2]]
    bad_row <- if (length(respondent_indices) >= na_loc[1]) respondent_indices[na_loc[1]] else NA_integer_
    msg <- sprintf(
      "Covariate '%s' contains NA values among respondents.%s",
      bad_col,
      if (is.finite(bad_row)) sprintf("\nFirst NA at row %d", bad_row) else ""
    )
    stop(msg, call. = FALSE)
  }
  invisible(NULL)
}

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

el_make_delta_column <- function(data, outcome_var, respondent_mask = NULL) {
  if (is.null(respondent_mask)) {
    respondent_mask <- !is.na(data[[outcome_var]])
  }
  if (length(respondent_mask) != nrow(data)) {
    stop("Internal error: respondent mask must align with data.", call. = FALSE)
  }

  delta_name <- "..nmar_delta.."
  if (delta_name %in% names(data)) {
    i <- 1L
    while (paste0(delta_name, i) %in% names(data)) i <- i + 1L
    delta_name <- paste0(delta_name, i)
  }

  data[[delta_name]] <- as.integer(respondent_mask)
  list(data = data, delta_name = delta_name)
}

el_rethrow_data_error <- function(err) {
  msg <- conditionMessage(err)
  missing_pattern <- "object '([^']+)' not found"
  if (grepl(missing_pattern, msg, perl = TRUE)) {
    missing_var <- sub(missing_pattern, "\\1", msg, perl = TRUE)
    stop(sprintf("Variables not found in data: %s", missing_var), call. = FALSE)
  }
  stop(msg, call. = FALSE)
}
