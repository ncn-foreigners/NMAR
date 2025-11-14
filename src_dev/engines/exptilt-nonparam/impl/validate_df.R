et_np_validate_df <- function(X, Y, Z, refuse_col) {

# --- Helper function and checks for empty data frames ---
# We define "empty" as NULL, or having 0 rows, or 0 columns
  is_df_empty <- function(df) {
    is.null(df) || nrow(df) == 0 || ncol(df) == 0
  }

  is_X_empty <- is_df_empty(X)
  is_Y_empty <- is_df_empty(Y)
  is_Z_empty <- is_df_empty(Z)

# --- 1. Validate 'refuse_col' itself ---
# Must be a single, non-empty, non-NA string
  if (is.null(refuse_col) ||
      !is.character(refuse_col) ||
      length(refuse_col) != 1 ||
      is.na(refuse_col) ||
      refuse_col == "") {
    stop(paste("Validation Error (1): 'refuse_col' must be a single,",
               "non-empty, non-NA string."))
  }

# --- 2. Validate Column Names (Disjoint Sets) ---
# Colnames across all *non-empty* dfs + refuse_col must be unique.

  all_sets <- list()
  if (!is_X_empty) all_sets$X <- colnames(X)
  if (!is_Y_empty) all_sets$Y <- colnames(Y)
  if (!is_Z_empty) all_sets$Z <- colnames(Z)

# Add the refuse_col as its own "set"
  all_sets$refuse <- refuse_col

# Combine all names into a single vector
  all_names_vector <- unlist(all_sets, use.names = FALSE)

# Check if the length of unique names is less than the total length
# (or, more simply, if any are duplicated)
  if (any(duplicated(all_names_vector))) {
    duplicates <- unique(all_names_vector[duplicated(all_names_vector)])
    stop(paste("Validation Error (2): Column names and 'refuse_col' must be disjoint.",
               "Found shared names:", paste(duplicates, collapse = ", ")))
  }

# --- 3. No NA values in any data frame ---
# This check only applies to non-empty data frames.

  if (!is_X_empty && any(is.na(X))) {
    stop("Validation Error (3): X contains NA values.")
  }

  if (!is_Y_empty && any(is.na(Y))) {
    stop("Validation Error (3): Y contains NA values.")
  }

  if (!is_Z_empty && any(is.na(Z))) {
    stop("Validation Error (3): Z contains NA values.")
  }

# --- 4. All columns (Expect Z) must be numeric ---
# This check only applies to non-empty data frames.

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
