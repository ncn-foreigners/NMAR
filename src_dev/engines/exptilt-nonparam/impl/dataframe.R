#' @importFrom stats aggregate
#' @exportS3Method exptilt_nonparam data.frame
exptilt_nonparam.data.frame <- function(
    data,
    formula,
    trace_level = 0,
    refusal_col, # This must be the name of a 0/1 *indicator* column
    max_iter = 100,
    tol = 1e-6,
    design_weights = NULL # <-- NEW: Argument for sampling weights
) {

# --- 1. EXTRACT and VALIDATE ---

  res = et_extract_formula(format(formula), data, allow_z_categorical = TRUE)
  Y = res$Y
  Y = Y[, !colnames(Y) %in% c("(Intercept)"), drop = FALSE]
  X = res$X
  Z = res$Z

  if (refusal_col %in% colnames(data) == FALSE) {
    stop(paste("Refusal column", refusal_col, "not found in data frame."))
  }
  if (is.null(Z)) {
    Z = matrix(1, nrow = nrow(X), ncol = 1)
    colnames(Z) = "(Intercept)"
  }
  outcome_cols <- colnames(Y)

  X = X[, !colnames(X) %in% c("(Intercept)", outcome_cols), drop = FALSE]
  Z = Z[, !colnames(Z) %in% c("(Intercept)", outcome_cols), drop = FALSE]

# et_np_validate_df(X, Y, Z, refusal_col) # Assumes you have this
  x1_cols = colnames(X)
  x2_cols = colnames(Z)

# --- 2. PREPARE KEYS and WEIGHTS (NEW LOGIC) ---

# Initialize weights if not provided
  if (is.null(design_weights)) {
    design_weights <- rep(1, nrow(data))
  }
  if (length(design_weights) != nrow(data)) {
    stop("Length of design_weights must match nrow(data).")
  }

# Add weights to the data for aggregation
  data$`_temp_weights_` <- design_weights

# Create a unique key for x1 columns (the response model)
  if (length(x1_cols) > 0) {
    data$x1_key <- apply(data[, x1_cols, drop = FALSE], 1, paste, collapse = "_")
  } else {
    data$x1_key <- "all"
  }

# Create a unique key for ALL x columns (x1 + x2)
  all_x_cols <- c(x1_cols, x2_cols)
  if (length(all_x_cols) > 0) {
    data$x_key <- apply(data[, all_x_cols, drop = FALSE], 1, paste, collapse = "_")
  } else {
    data$x_key <- "all"
  }

# Map for relating x_key back to x1_key
  x_to_x1_map <- unique(data.frame(x_key = data$x_key, x1_key = data$x1_key))
  rownames(x_to_x1_map) <- x_to_x1_map$x_key

# Get unique x1 keys in a stable order
  unique_x1_keys <- unique(data$x1_key)

# --- 3. CALCULATE FIXED MATRICES (WEIGHTED) ---

# "n matrix" (n_{y*x*}) - WEIGHTED respondent counts
# Assumes Y columns are 0/1 indicators
  n_y_x_data <- data.frame(
    x_key = data$x_key,
# Multiply each outcome by the weight
    as.matrix(data[, outcome_cols]) * data$`_temp_weights_`
  )
  n_y_x_matrix <- stats::aggregate(. ~ x_key, data = n_y_x_data, FUN = sum)
  rownames(n_y_x_matrix) <- n_y_x_matrix$x_key
# Ensure matrix rows are in a consistent order
  n_y_x_matrix <- as.matrix(n_y_x_matrix[x_to_x1_map$x_key, -1])

# "m vector" (m_{x*}) - WEIGHTED nonrespondent counts
# Assumes refusal_col is a 0/1 indicator
  m_x_data <- data.frame(
    x_key = data$x_key,
    m_x = data[, refusal_col] * data$`_temp_weights_` # Apply weights
  )
  m_x_vec_df <- stats::aggregate(m_x ~ x_key, data = m_x_data, FUN = sum)
  rownames(m_x_vec_df) <- m_x_vec_df$x_key
# Ensure vector is in the same consistent order
  m_x_vec <- m_x_vec_df[x_to_x1_map$x_key, "m_x", drop = TRUE]

# "p matrix" (p_hat_{y*|x*}) - Respondent proportions (from weighted counts)
  total_resp_x <- rowSums(n_y_x_matrix)
  total_resp_x[total_resp_x == 0] <- 1 # Avoid division by zero
  p_hat_matrix <- n_y_x_matrix / total_resp_x

# "n_y_x1 matrix" (n_{y*x1*}) - WEIGHTED respondent counts, summed by x1
  n_y_x1_data <- data.frame(
    x1_key = data$x1_key,
    as.matrix(data[, outcome_cols]) * data$`_temp_weights_` # Apply weights
  )
  n_y_x1_matrix <- stats::aggregate(. ~ x1_key, data = n_y_x1_data, FUN = sum)
  rownames(n_y_x1_matrix) <- n_y_x1_matrix[, 1]
# Ensure matrix rows are in a consistent order
  n_y_x1_matrix <- as.matrix(n_y_x1_matrix[unique_x1_keys, -1])

# --- 4. THE EM ALGORITHM (The "Loop") ---
# This part is identical to your notebook code, as it now runs on
# the correctly weighted, aggregated matrices.

  odds_matrix <- matrix(1.0,
                        nrow = nrow(n_y_x1_matrix),
                        ncol = ncol(n_y_x1_matrix),
                        dimnames = dimnames(n_y_x1_matrix))

  m_y_x1_matrix <- matrix(NA,
                          nrow = nrow(n_y_x1_matrix),
                          ncol = ncol(n_y_x1_matrix),
                          dimnames = dimnames(n_y_x1_matrix))

# Get the x1_key for each row of the p_hat_matrix
  p_hat_x1_keys <- x_to_x1_map[rownames(p_hat_matrix), "x1_key"]

  for (iter in 1:max_iter) {
    odds_old <- odds_matrix

# E-STEP (Calculate myx):
    odds_joined_matrix <- odds_matrix[p_hat_x1_keys, ]
    denominator <- rowSums(p_hat_matrix * odds_joined_matrix)
    denominator[denominator == 0] <- 1
    m_y_x_matrix <- m_x_vec * (p_hat_matrix * odds_joined_matrix) / denominator

# M-STEP (Calculate O):
    m_y_x_agg_data <- data.frame(x1_key = p_hat_x1_keys, m_y_x_matrix)
    m_y_x1_matrix <- stats::aggregate(. ~ x1_key, data = m_y_x_agg_data, FUN = sum)
    rownames(m_y_x1_matrix) <- m_y_x1_matrix[, 1]
    m_y_x1_matrix <- as.matrix(m_y_x1_matrix[unique_x1_keys, -1]) # Order rows

    n_safe <- n_y_x1_matrix
    n_safe[n_safe == 0] <- 1
    odds_matrix <- m_y_x1_matrix / n_safe

    diff <- sum(abs(odds_matrix - odds_old), na.rm = TRUE)
    if (diff < tol) break
  }

  if (iter == max_iter) {
    warning(paste("Algorithm did not converge after", max_iter, "iterations."))
  }

# --- 5. CALCULATE FINAL m_y_x & RETURN ---
  odds_joined_matrix <- odds_matrix[p_hat_x1_keys, ]
  denominator <- rowSums(p_hat_matrix * odds_joined_matrix)
  denominator[denominator == 0] <- 1
  m_y_x_matrix <- m_x_vec * (p_hat_matrix * odds_joined_matrix) / denominator

# Add back the original x_key identifiers
  rownames(m_y_x_matrix) <- rownames(p_hat_matrix)

# Create a final adjusted data.frame (aggregated)
  data_to_return <- as.data.frame(n_y_x_matrix + m_y_x_matrix)
# Add back the x1 and x2 columns
  data_to_return <- cbind(x_to_x1_map[rownames(n_y_x_matrix), ], data_to_return)

# return(list(
#   ... (your old list return) ...
# ))

  return(new_nmar_result_exptilt_nonparam(
    n_y_x_matrix = n_y_x_matrix,
    m_x_vec = m_x_vec,
    p_hat_y_given_x_matrix = p_hat_matrix,
    data_to_return = data_to_return,
    x1_cols = x1_cols,
    x2_cols = x2_cols,
    outcome_cols = outcome_cols,
    engine_params = list(
      refusal_col = refusal_col,
      max_iter = max_iter,
      tol = tol
    ),
    fit_stats = list(
      iterations = iter,
      converged = (iter < max_iter)
    ),
    loss_value = diff
  ))
}
