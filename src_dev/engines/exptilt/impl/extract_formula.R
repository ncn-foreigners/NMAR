#' @importFrom stats formula model.frame na.pass terms
et_extract_formula <- function(formula_str, data, allow_z_categorical = FALSE) {
# browser()


  f <- Formula::Formula(as.formula(formula_str))
  num_rhs_parts <- length(f)[2]

# --- Y:  (LHS) ---
  mf_Y <- model.frame(f, data = data, lhs = 1, rhs = 0, na.action = na.pass)
  Y_mat <- data.matrix(mf_Y)

# --- X:  (RHS), Part 1 ---
  f_X <- formula(f, lhs = 0, rhs = 1)
  X_terms <- attr(terms(f_X), "term.labels")

  X_mat <- model.matrix(f_X, data = data)
  X_mat <- X_mat[, colnames(X_mat) != "(Intercept)", drop = FALSE]

# Check for name mismatch
  if (ncol(X_mat) == length(X_terms)) {
    colnames(X_mat) <- X_terms
  }

# --- Z: (RHS), Part 2 ---
  if (num_rhs_parts > 1) {
    f_Z <- formula(f, lhs = 0, rhs = 2)
    Z_terms <- attr(terms(f_Z), "term.labels")
    if (allow_z_categorical) {
     Z_mat = data[, Z_terms, drop = FALSE]
    }
    else {

# Z_mat <- model.matrix(f_Z, data = data, na.action = na.pass)
    mf_Z <- model.frame(f_Z, data = data, na.action = na.pass)
    Z_mat <- model.matrix(f_Z, data = mf_Z)

# This removes the intercept, leaving the two broken columns
    Z_mat <- Z_mat[, colnames(Z_mat) != "(Intercept)", drop = FALSE]
}

# It checks if 2 columns exist (they do: "" and "x3")
# It then renames them to c("Y", "x3")
    if (ncol(Z_mat) == length(Z_terms)) {
      colnames(Z_mat) <- Z_terms
    } else {
      warning("Column/Term mismatch in Z matrix. Names may be incorrect.")
    }

  } else {
    Z_mat <- NULL
  }

  return(list(Y = Y_mat, X = X_mat, Z = Z_mat))
}
