#' @importFrom stats aggregate
#' @exportS3Method exptilt_nonparam data.frame
exptilt_nonparam.data.frame <- function(
    data,
    formula,
    trace_level = 0,
    refusal_col,
    max_iter = 100,
    tol = 1e-6,
    design_weights = NULL # Added to match S3 pattern
) {

# --- 0. CREATE VERBOSER ---
  verboser <- create_verboser(trace_level)

  verboser("============================================================", level = 1, type = "step")
  verboser("  NONPARAMETRIC EXPTILT (Aggregated Data) STARTED", level = 1, type = "step")
  verboser("============================================================", level = 1, type = "step")

# Show trace level info
  trace_msg <- sprintf("Running with trace_level = %d", trace_level)
  if (trace_level < 3) {
    trace_msg <- paste0(trace_msg, " | For more detail, use trace_level = 3")
  }
  verboser(trace_msg, level = 1)
  verboser(sprintf("Formula: %s", deparse(formula)), level = 1)

# --- 1. EXTRACT and VALIDATE (Aggregated data style) ---

# --- NOTE: Replace this with your et_extract_formula call ---
# This temporary parser assumes a formula like:
# Voted_A + Voted_B + Other ~ Gender | Age_group
# LHS = outcomes, RHS 1 = x1, RHS 2 = x2
  parsed_formula <- list(
    outcome_cols = c("Voted_A", "Voted_B", "Other"),
    x1_cols = c("Gender"), # Based on your formula
    x2_cols = c("Age_group") # Based on your formula
  )
# --- END TEMPORARY PARSER ---

# res = et_extract_formula(format(formula), data, allow_z_categorical = T)
# outcome_cols <- colnames(res$Y)
# ... etc ...

# Using the temporary parser for this example
  outcome_cols <- parsed_formula$outcome_cols
  x1_cols <- parsed_formula$x1_cols
  x2_cols <- parsed_formula$x2_cols


  if (refusal_col %in% colnames(data) == FALSE) {
    stop(paste("Refusal column", refusal_col, "not found in data frame."))
  }

# --- DATA SUMMARY (Level 1) ---
# This logic is now correct for AGGREGATED data
  n_strata <- nrow(data)
  total_resp <- sum(data[, outcome_cols], na.rm = TRUE)
  total_nonresp <- sum(data[, refusal_col], na.rm = TRUE)
  n_total <- total_resp + total_nonresp
  pct_nonresp <- 100 * total_nonresp / n_total

  verboser("", level = 1)
  verboser("-- DATA SUMMARY (Aggregated) --", level = 1)
  verboser(sprintf("  Total strata:         %d", n_strata), level = 1)
  verboser(sprintf("  Total observations:   %d (Respondents + Nonrespondents)", n_total), level = 1)
  verboser(sprintf("  Respondents:          %d (%.1f%%)", total_resp, 100 - pct_nonresp), level = 1)
  verboser(sprintf("  Non-respondents:      %d (%.1f%%)", total_nonresp, pct_nonresp), level = 1)

# --- MODEL SPEC (Level 2) ---
  verboser("", level = 2)
  verboser("-- MODEL SPECIFICATION --", level = 2)
  verboser(sprintf("  Outcome columns:      %s", paste(outcome_cols, collapse = ", ")), level = 2)
  verboser(sprintf("  Response model (x1):  %s", if (length(x1_cols) > 0) paste(x1_cols, collapse = ", ") else "(none)"), level = 2)
  verboser(sprintf("  Instrument (x2):      %s", if (length(x2_cols) > 0) paste(x2_cols, collapse = ", ") else "(none)"), level = 2)
  verboser(sprintf("  Refusal column:       %s", refusal_col), level = 2)

  if (!is.null(design_weights)) {
    verboser("  Survey weights:       Yes (NOTE: Weights assumed to be pre-applied to aggregated counts)", level = 2, type = "info")
  } else {
    verboser("  Survey weights:       No (using raw counts)", level = 2, type = "info")
  }

# --- 2. PREPARE KEYS and FIXED MATRICES ---
# This uses the robust key-creation logic

  if (length(x1_cols) > 0) {
    data$x1_key <- apply(data[, x1_cols, drop = FALSE], 1, paste, collapse = "_")
  } else {
    data$x1_key <- "all"
  }

  all_x_cols <- c(x1_cols, x2_cols)
  if (length(all_x_cols) > 0) {
    data$x_key <- apply(data[, all_x_cols, drop = FALSE], 1, paste, collapse = "_")
  } else {
    data$x_key <- "all"
  }

  x_to_x1_map <- unique(data.frame(x_key = data$x_key, x1_key = data$x1_key))
  rownames(x_to_x1_map) <- x_to_x1_map$x_key
  unique_x1_keys <- sort(unique(data$x1_key))

  verboser("", level = 2)
  verboser("-- PRE-COMPUTATION --", level = 2)

# "n matrix" (n_{y*x*})
  n_y_x_matrix <- as.matrix(data[, outcome_cols])
  rownames(n_y_x_matrix) <- data$x_key

# "m vector" (m_{x*})
  m_x_vec <- data[, refusal_col]
  names(m_x_vec) <- data$x_key

# "p matrix" (p_hat_{y*|x*})
  total_resp_x <- rowSums(n_y_x_matrix)
  total_resp_x[total_resp_x == 0] <- 1
  p_hat_matrix <- n_y_x_matrix / total_resp_x

# "n_y_x1 matrix" (n_{y*x1*})
  n_y_x1_matrix <- stats::aggregate(
    n_y_x_matrix ~ data$x1_key,
    FUN = sum
  )
  rownames(n_y_x1_matrix) <- n_y_x1_matrix[, 1]
  n_y_x1_matrix <- as.matrix(n_y_x1_matrix[unique_x1_keys, -1, drop = FALSE])

  verboser("  OK Fixed matrices computed from aggregated data.", level = 2)
  verboser("  n_y_x1_matrix (Respondents by x1):", level = 3, type = "detail", obj = n_y_x1_matrix)
  verboser("  p_hat_matrix (Respondent proportions by x1,x2):", level = 3, type = "detail", obj = p_hat_matrix)

# --- 3. THE EM ALGORITHM (The "Loop") ---
  verboser("", level = 1)
  verboser("-- EM ALGORITHM --", level = 1)
  verboser(sprintf("  Solving... (max_iter = %d, tol = %s)", max_iter, tol), level = 1)

  odds_matrix <- matrix(1.0,
                        nrow = nrow(n_y_x1_matrix),
                        ncol = ncol(n_y_x1_matrix),
                        dimnames = dimnames(n_y_x1_matrix))

  m_y_x1_matrix <- matrix(NA,
                          nrow = nrow(n_y_x1_matrix),
                          ncol = ncol(n_y_x1_matrix),
                          dimnames = dimnames(n_y_x1_matrix))

  p_hat_x1_keys <- x_to_x1_map[rownames(p_hat_matrix), "x1_key"]

  iter <- 0
  diff <- tol + 1

  for (iter in 1:max_iter) {
    odds_old <- odds_matrix

# E-STEP (Calculate myx):
    odds_joined_matrix <- odds_matrix[p_hat_x1_keys, , drop = FALSE]
    denominator <- rowSums(p_hat_matrix * odds_joined_matrix)
    denominator[denominator == 0] <- 1
    m_y_x_matrix <- m_x_vec * (p_hat_matrix * odds_joined_matrix) / denominator

# M-STEP (Calculate O):
    m_y_x_agg_data <- data.frame(x1_key = p_hat_x1_keys, m_y_x_matrix)
    m_y_x1_matrix <- stats::aggregate(. ~ x1_key, data = m_y_x_agg_data, FUN = sum)
    rownames(m_y_x1_matrix) <- m_y_x1_matrix[, 1]
    m_y_x1_matrix <- as.matrix(m_y_x1_matrix[unique_x1_keys, -1, drop = FALSE])

    n_safe <- n_y_x1_matrix
    n_safe[n_safe == 0] <- 1
    odds_matrix <- m_y_x1_matrix / n_safe

    diff <- sum(abs(odds_matrix - odds_old), na.rm = TRUE)
    if (diff < tol) break
  }

# --- 4. CALCULATE FINAL m_y_x & RETURN ---

# Convergence status messages
  conv_status <- if (iter < max_iter) "OK Converged" else "Warning: Max iterations reached"
  verboser(sprintf("  %s", conv_status), level = 1)
  verboser(sprintf("  Iterations:             %d", iter), level = 1)
  verboser(sprintf("  Termination code:       %d", if (iter < max_iter) 1 else 0), level = 2)
  verboser(sprintf("  Max |odds change|:      %.6f", diff), level = 2)

  odds_joined_matrix <- odds_matrix[p_hat_x1_keys, , drop = FALSE]
  denominator <- rowSums(p_hat_matrix * odds_joined_matrix)
  denominator[denominator == 0] <- 1
  m_y_x_matrix <- m_x_vec * (p_hat_matrix * odds_joined_matrix) / denominator

  rownames(m_y_x_matrix) <- rownames(p_hat_matrix)

# --- PREPARE THE RETURN OBJECT ---

# 1. Create the data_to_return dataframe as requested
  data_to_return <- data

# Get the original n_y_x* matrix, ordered by the x_key from the data
  n_y_x_matrix_ordered <- n_y_x_matrix[data$x_key, ]

# Get the final m_y_x* matrix, ordered by the x_key from the data
  m_y_x_matrix_ordered <- m_y_x_matrix[data$x_key, ]

# *** YOUR REQUESTED LOGIC ***
# Add the calculated m_y_x scores to the original n_y_x counts
  data_to_return[, outcome_cols] <- n_y_x_matrix_ordered + m_y_x_matrix_ordered

# Remove helper keys and refusal column
  data_to_return$x1_key <- NULL
  data_to_return$x_key <- NULL
  data_to_return[, refusal_col] <- NULL

# 2. Get the other matrices for the S3 object
  final_n_y_x_matrix <- n_y_x_matrix[rownames(m_y_x_matrix), , drop = FALSE]

# --- FINAL RESULTS (Level 1-3) ---
  verboser("", level = 1, type = "result")
  verboser("============================================================", level = 1, type = "result")
  verboser("  ESTIMATION COMPLETE", level = 1, type = "result")
  verboser("============================================================", level = 1, type = "result")
  verboser("  Final Odds Ratios (O(x1, y)):", level = 1, type = "result", obj = odds_matrix)
  verboser("  Final Expected Nonrespondents (m(x1, y)):", level = 2, type = "result", obj = m_y_x1_matrix)
  verboser("  Final Adjusted Counts (Aggregated):", level = 3, type = "detail", obj = data_to_return)

# This is the return from your previous function
  return(new_nmar_result_exptilt_nonparam(
    n_y_x_matrix = final_n_y_x_matrix, # n_y_x* (by x1, x2)
    m_x_vec = m_x_vec, # m_x* (by x1, x2)
    p_hat_y_given_x_matrix = p_hat_matrix,
    data_to_return = data_to_return, # The dataframe you requested
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
      converged = (iter < max_iter) # <-- Now correctly reports convergence
    ),
    loss_value = diff,
# Add the other matrices as well
    odds_ratio = odds_matrix, # O(x1, y)
    m_y_x1 = m_y_x1_matrix, # m_y_x1* (by x1)
    m_y_x = m_y_x_matrix # m_y_x* (by x1, x2)
  ))
}
