#' @importFrom stats aggregate
#' @exportS3Method exptilt_nonparam data.frame
exptilt_nonparam.data.frame <- function(
    data,
    formula,
    trace_level = 0,
    refusal_col,
    max_iter = 100,
    tol_value = 1e-6,
    design_weights = NULL,
    ...
) {


  verboser <- create_verboser(trace_level)
  verboser("============================================================", level = 1, type = "step")
  verboser("  NONPARAMETRIC EXPTILT (Aggregated Data) STARTED", level = 1, type = "step")
  verboser("============================================================", level = 1, type = "step")


  if (refusal_col %in% colnames(data) == FALSE) {
    stop(paste("Refusal column", refusal_col, "not found in data frame."))
  }


  res <- et_extract_formula(format(formula), data, allow_z_categorical = TRUE)


  Y <- as.data.frame(res$Y)
  X <- as.data.frame(res$X)

  if (!is.null(res$Z)) {
    Z <- as.data.frame(res$Z)
  } else {
# if no instrument provided, create dummy intercept
    Z <- data.frame(`(Intercept)` = rep(1, nrow(data)))
    colnames(Z) <- "(Intercept)"
  }

  et_np_validate_df(X, Y, Z, refuse_col = refusal_col)

  outcome_cols <- colnames(Y)
  x1_cols <- colnames(X)
  x2_cols <- colnames(Z)


  data_orig <- data[, c(outcome_cols, x1_cols, x2_cols, refusal_col)]

# build working dataset (numeric X/Y, potentially categorical Z)
  refusal_vec <- data[[refusal_col]]
  data_working <- cbind(Y, X, Z)
  data_working[[refusal_col]] <- refusal_vec

# data summary for stats
  n_strata <- nrow(data_working)
  total_resp <- sum(data_working[, outcome_cols], na.rm = TRUE)
  total_nonresp <- sum(data_working[, refusal_col], na.rm = TRUE)
  n_total <- total_resp + total_nonresp
  pct_nonresp <- 100 * total_nonresp / n_total

  verboser("", level = 1)
  verboser("-- DATA SUMMARY (Aggregated) --", level = 1)
  verboser(sprintf("  Total strata:         %d", n_strata), level = 1)
  verboser(sprintf("  Total observations:   %d", n_total), level = 1)

  verboser("", level = 2)
  verboser("-- MODEL SPECIFICATION --", level = 2)
  verboser(sprintf("  Outcome columns:      %s", paste(outcome_cols, collapse = ", ")), level = 2)

  if (!is.null(design_weights)) {
    verboser("  Survey weights:       Yes (scaling n and m counts by weights)", level = 2, type = "info")
  } else {
    verboser("  Survey weights:       No (using raw counts)", level = 2, type = "info")
  }


  if (is.null(design_weights)) {
    design_weights <- rep(1, nrow(data_working))
  }
  if (length(design_weights) != nrow(data_working)) {
    stop("Length of design_weights does not match number of rows in data.")
  }

  if (length(x1_cols) > 0) {
    data_working$x1_key <- apply(data_working[, x1_cols, drop = FALSE], 1, paste, collapse = "_")
  } else {
    data_working$x1_key <- "all"
  }

  all_x_cols <- c(x1_cols, x2_cols)
  if (length(all_x_cols) > 0) {
    data_working$x_key <- apply(data_working[, all_x_cols, drop = FALSE], 1, paste, collapse = "_")
  } else {
    data_working$x_key <- "all"
  }

  x_to_x1_map <- unique(data.frame(x_key = data_working$x_key, x1_key = data_working$x1_key))
  rownames(x_to_x1_map) <- x_to_x1_map$x_key
  unique_x1_keys <- sort(unique(data_working$x1_key))
  unique_x_keys <- sort(unique(data_working$x_key))

  verboser("", level = 2)
  verboser("-- PRE-COMPUTATION (Weighted Aggregation) --", level = 2)

# apply weights
  weighted_n_data <- as.matrix(data_working[, outcome_cols]) * design_weights
  weighted_m_data <- data_working[[refusal_col]] * design_weights

  agg_df <- data.frame(
    x_key = data_working$x_key,
    x1_key = data_working$x1_key,
    m_weighted = weighted_m_data
  )
  agg_df <- cbind(agg_df, weighted_n_data)

# n matrix
  n_y_x_df <- stats::aggregate(. ~ x_key, data = agg_df[, c("x_key", outcome_cols)], FUN = sum)
  rownames(n_y_x_df) <- n_y_x_df$x_key
  n_y_x_matrix <- as.matrix(n_y_x_df[unique_x_keys, outcome_cols])

# m vector
  m_x_df <- stats::aggregate(m_weighted ~ x_key, data = agg_df[, c("x_key", "m_weighted")], FUN = sum)
  rownames(m_x_df) <- m_x_df$x_key
  m_x_vec <- m_x_df[unique_x_keys, "m_weighted"]
  names(m_x_vec) <- unique_x_keys

# p matrix
  total_resp_x <- rowSums(n_y_x_matrix)
  total_resp_x[total_resp_x == 0] <- 1
  p_hat_matrix <- n_y_x_matrix / total_resp_x

# n_y_x1 matrix
  n_y_x1_df <- stats::aggregate(. ~ x1_key, data = agg_df[, c("x1_key", outcome_cols)], FUN = sum)
  rownames(n_y_x1_df) <- n_y_x1_df$x1_key
  n_y_x1_matrix <- as.matrix(n_y_x1_df[unique_x1_keys, outcome_cols])

  verboser("  OK Fixed matrices computed.", level = 2)

# EM
  verboser("", level = 1)
  verboser("-- EM ALGORITHM --", level = 1)
  verboser(sprintf("  Solving... (max_iter = %d, tol_value = %s)", max_iter, tol_value), level = 1)

  odds_matrix <- matrix(1.0, nrow = nrow(n_y_x1_matrix), ncol = ncol(n_y_x1_matrix), dimnames = dimnames(n_y_x1_matrix))
  p_hat_x1_keys <- x_to_x1_map[rownames(p_hat_matrix), "x1_key"]

  iter <- 0
  diff <- tol_value + 1

  for (iter in 1:max_iter) {
    odds_old <- odds_matrix

# E
    odds_joined_matrix <- odds_matrix[p_hat_x1_keys, , drop = FALSE]
    denominator <- rowSums(p_hat_matrix * odds_joined_matrix)
    denominator[denominator == 0] <- 1
    m_y_x_matrix <- m_x_vec * (p_hat_matrix * odds_joined_matrix) / denominator

# M
    m_agg_df <- as.data.frame(m_y_x_matrix)
    m_agg_df$x1_key <- p_hat_x1_keys
    m_y_x1_df <- stats::aggregate(. ~ x1_key, data = m_agg_df, FUN = sum)
    rownames(m_y_x1_df) <- m_y_x1_df$x1_key
    m_y_x1_matrix <- as.matrix(m_y_x1_df[unique_x1_keys, colnames(odds_matrix)])

    n_safe <- n_y_x1_matrix
    n_safe[n_safe == 0] <- 1
    odds_matrix <- m_y_x1_matrix / n_safe

    diff <- sum(abs(odds_matrix - odds_old), na.rm = TRUE)
    if (diff < tol_value) break
  }


  odds_joined_matrix <- odds_matrix[p_hat_x1_keys, , drop = FALSE]
  denominator <- rowSums(p_hat_matrix * odds_joined_matrix)
  denominator[denominator == 0] <- 1
  m_y_x_matrix <- m_x_vec * (p_hat_matrix * odds_joined_matrix) / denominator
  rownames(m_y_x_matrix) <- rownames(p_hat_matrix)

  data_to_return <- data_orig

  m_y_x_mapped <- m_y_x_matrix[data_working$x_key, ]

  original_counts <- as.matrix(data_to_return[, outcome_cols])

  if (!is.null(design_weights) && !all(design_weights == 1)) {
    original_counts <- original_counts * design_weights
  }

  data_to_return[, outcome_cols] <- original_counts + m_y_x_mapped

  data_to_return[, refusal_col] <- NULL

  final_n_y_x_matrix <- n_y_x_matrix[rownames(m_y_x_matrix), , drop = FALSE]

  verboser("", level = 1, type = "result")
  verboser("============================================================", level = 1, type = "result")
  verboser("  ESTIMATION COMPLETE", level = 1, type = "result")
  verboser("============================================================", level = 1, type = "result")
  verboser("  Final Adjusted Counts (Aggregated):", level = 3, type = "detail", obj = data_to_return)

  return(new_nmar_result_exptilt_nonparam(
    n_y_x_matrix = final_n_y_x_matrix,
    m_x_vec = m_x_vec,
    p_hat_y_given_x_matrix = p_hat_matrix,
    data_to_return = data_to_return,
    x1_cols = x1_cols,
    x2_cols = x2_cols,
    outcome_cols = outcome_cols,
    engine_params = list(
      refusal_col = refusal_col,
      max_iter = max_iter,
      tol_value = tol_value
    ),
    fit_stats = list(
      iterations = iter,
      converged = (iter < max_iter)
    ),
    loss_value = diff,
    odds_ratio = odds_matrix,
    m_y_x1 = m_y_x1_matrix,
    m_y_x = m_y_x_matrix
  ))
}
