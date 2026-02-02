#' @exportS3Method NULL

et_validate_df <- function(X, Y, Z) {

# --- Helper checks for empty data frames ---
# We define "empty" as having 0 rows or 0 columns
  is_X_empty <- is.null(X) || nrow(X) == 0 || ncol(X) == 0
  is_Y_empty <- is.null(Y) || nrow(Y) == 0 || ncol(Y) == 0
  is_Z_empty <- is.null(Z) || nrow(Z) == 0 || ncol(Z) == 0

# --- 1. Validate Column Names (Disjoint Sets) ---
# Colnames across all *non-empty* data frames must be unique.
# No column name can appear in more than one data frame.

  all_sets <- list()
  if (!is_X_empty) all_sets$X <- colnames(X)
  if (!is_Y_empty) all_sets$Y <- colnames(Y)
  if (!is_Z_empty) all_sets$Z <- colnames(Z)

# Only check if there is more than one non-empty df
  if (length(all_sets) > 1) {
# Combine all column names into a single vector
    all_names_vector <- unlist(all_sets, use.names = FALSE)

# Check if there are any duplicates in the combined vector
    if (any(duplicated(all_names_vector))) {

# Find the specific duplicates for a better error message
      duplicates <- unique(all_names_vector[duplicated(all_names_vector)])

      stop(paste("Validation Error (1): Column names must be disjoint (unique) across X, Y, and Z.",
                 "Found shared columns:", paste(duplicates, collapse = ", ")))
    }
  }

# --- 2. At least 1 NA in Y ---
# This check only applies if Y is *not* empty.

  if (!is_Y_empty) {
    if (!any(is.na(Y))) {
      stop("Validation Error (2): Y is not empty but contains no NA values.")
    }
  }

# --- 3. At least one of X or Z is not empty ---
# Fails if *both* X and Z are empty.

  if (is_X_empty && is_Z_empty) {
    stop("Validation Error (3): Both X and Z are empty. At least one must be non-empty.")
  }

# --- 4. Any NA in X, Z ---
# This is interpreted as a *failure* condition.
# Fail if any NA is found in non-empty X or non-empty Z.

  if (!is_X_empty && any(is.na(X))) {
    stop("Validation Error (4): X contains NA values.")
  }

  if (!is_Z_empty && any(is.na(Z))) {
    stop("Validation Error (4): Z contains NA values.")
  }

# --- 5. Finite numeric values in X, Y (observed), Z ---
# The exptilt estimator assumes numeric, finite design matrices. NA handling in
# Y is special: NA values are allowed and required (at least one NA) to encode
# nonresponse; observed Y values must be finite.

  if (!is_X_empty) {
    if (!is.numeric(X)) {
      stop("Validation Error (5): X must be numeric.")
    }
    if (any(!is.finite(X))) {
      stop("Validation Error (5): X contains non-finite values.")
    }
  }

  if (!is_Z_empty) {
    if (!is.numeric(Z)) {
      stop("Validation Error (5): Z must be numeric.")
    }
    if (any(!is.finite(Z))) {
      stop("Validation Error (5): Z contains non-finite values.")
    }
  }

  if (!is_Y_empty) {
    if (!is.numeric(Y)) {
      stop("Validation Error (5): Y must be numeric.")
    }
    y_obs <- Y[!is.na(Y)]
    if (length(y_obs) > 0 && any(!is.finite(y_obs))) {
      stop("Validation Error (5): Observed Y contains non-finite values.")
    }
  }

# --- Success ---
#   message("Validation successful.")
  return(invisible(TRUE))
}
