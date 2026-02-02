et_np_validate_df <- function(X, Y, Z, refuse_col) {



  is_df_empty <- function(df) {
    is.null(df) || nrow(df) == 0 || ncol(df) == 0
  }

  is_X_empty <- is_df_empty(X)
  is_Y_empty <- is_df_empty(Y)
  is_Z_empty <- is_df_empty(Z)

# --- validate 'refuse_col' ---
# Must be a single, non-empty, non-NA string
  if (is.null(refuse_col) ||
      !is.character(refuse_col) ||
      length(refuse_col) != 1 ||
      is.na(refuse_col) ||
      refuse_col == "") {
    stop(paste("Validation Error (1): 'refuse_col' must be a single,",
               "non-empty, non-NA string."))
  }

# --- validate column names (disjoint sets) ---

  all_sets <- list()
  if (!is_X_empty) all_sets$X <- colnames(X)
  if (!is_Y_empty) all_sets$Y <- colnames(Y)
  if (!is_Z_empty) all_sets$Z <- colnames(Z)

  all_sets$refuse <- refuse_col

  all_names_vector <- unlist(all_sets, use.names = FALSE)

  if (any(duplicated(all_names_vector))) {
    duplicates <- unique(all_names_vector[duplicated(all_names_vector)])
    stop(paste("Validation Error (2): Column names and 'refuse_col' must be disjoint.",
               "Found shared names:", paste(duplicates, collapse = ", ")))
  }

# --- no NA values in any data frame ---

  if (!is_X_empty && any(is.na(X))) {
    stop("Validation Error (3): X contains NA values.")
  }

  if (!is_Y_empty && any(is.na(Y))) {
    stop("Validation Error (3): Y contains NA values.")
  }

  if (!is_Z_empty && any(is.na(Z))) {
    stop("Validation Error (3): Z contains NA values.")
  }

# --- no non-finite values in numeric inputs ---

  if (!is_X_empty && any(!is.finite(as.matrix(X)))) {
    stop("Validation Error (3b): X contains non-finite values.")
  }

  if (!is_Y_empty && any(!is.finite(as.matrix(Y)))) {
    stop("Validation Error (3b): Y contains non-finite values.")
  }

# --- all columns (expect Z) must be numeric ---

  if (!is_X_empty && !all(sapply(X, is.numeric))) {
    stop("Validation Error (4): Not all columns in X are numeric.")
  }

  if (!is_Y_empty && !all(sapply(Y, is.numeric))) {
    stop("Validation Error (4): Not all columns in Y are numeric.")
  }

# if (!is_Z_empty && !all(sapply(Z, is.numeric))) {
#   stop("Validation Error (4): Not all columns in Z are numeric.")
# }

# --- Success ---
# message("Validation successful.")
  return(invisible(TRUE))
}
