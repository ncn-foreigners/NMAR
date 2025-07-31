run_em_nmar_nonparametric <- function(
    data,
    outcome_cols,
    refusal_col,
    common_covariates,
    instrumental_covariates,
    max_iter = 100,
    tol = 1e-6
) {
  # --- Data Validation and Preparation ---
  required_cols <- c(outcome_cols, refusal_col, common_covariates, instrumental_covariates)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in data frame:", paste(missing_cols, collapse = ", ")))
  }

  # Create a unique key for common_covariates (i1) - base R alternative to unite()
  data_processed <- data
  data_processed$common_covariate_key <- apply(data_processed[, common_covariates, drop = FALSE],
                                               1, paste, collapse = "_")

  # Calculate p_{j|i} (observed outcome probabilities) - base R alternative
  total_observed <- rowSums(data_processed[, outcome_cols], na.rm = TRUE)
  p_cols <- paste0("p_", outcome_cols)

  for (col in outcome_cols) {
    data_processed[[paste0("p_", col)]] <- ifelse(total_observed == 0, 0,
                                                  data_processed[[col]] / total_observed)
  }

  outcome_classes <- outcome_cols
  common_covariate_keys <- unique(data_processed$common_covariate_key)

  # --- Initialize O_values (O_{j, i1}^{(t)}) ---
  O_values <- matrix(1.0, nrow = length(outcome_classes), ncol = length(common_covariate_keys),
                     dimnames = list(outcome_classes, common_covariate_keys))

  cat("--- Starting EM Algorithm ---\n")
  cat("Initial O_values (t=0):\n")
  print(O_values)

  # --- Main EM loop ---
  for (iter in 1:max_iter) {
    # E-step: Estimate m_{ij}^{(t)}
    P_j_i_matrix <- as.matrix(data_processed[, p_cols])
    colnames(P_j_i_matrix) <- outcome_classes

    Refusal_i_vector <- data_processed[[refusal_col]]
    current_keys <- data_processed$common_covariate_key

    # Base R alternative to dplyr operations
    O_for_each_row <- t(O_values[, current_keys, drop = FALSE])
    numerator_components <- P_j_i_matrix * O_for_each_row
    denominator_vector <- rowSums(numerator_components)

    weights_for_refusal_split <- Refusal_i_vector / denominator_vector
    weights_for_refusal_split[is.na(weights_for_refusal_split) | is.infinite(weights_for_refusal_split)] <- 0

    m_ij_t <- numerator_components * weights_for_refusal_split
    m_ij_t[is.na(m_ij_t) | is.infinite(m_ij_t)] <- 0

    # M-step: Update O_{j, i1}^{(t+1)}
    O_values_next <- matrix(0.0, nrow = length(outcome_classes), ncol = length(common_covariate_keys),
                            dimnames = list(outcome_classes, common_covariate_keys))

    # Base R alternative to group_by/summarize
    df_with_m_ij_t <- cbind(data_processed, m_ij_t)
    colnames(df_with_m_ij_t)[(ncol(df_with_m_ij_t)-length(outcome_classes)+1):ncol(df_with_m_ij_t)] <-
      paste0("m_", outcome_classes)

    for (key in common_covariate_keys) {
      subset_data <- df_with_m_ij_t[df_with_m_ij_t$common_covariate_key == key, ]

      for (oc in outcome_classes) {
        sum_m <- sum(subset_data[[paste0("m_", oc)]], na.rm = TRUE)
        sum_n <- sum(subset_data[[oc]], na.rm = TRUE)

        O_values_next[oc, key] <- ifelse(sum_n == 0, 0, sum_m / sum_n)
      }
    }

    # Check for convergence
    if (iter > 1 && max(abs(O_values_next - O_values)) < tol) {
      cat(paste("Algorithm converged at iteration", iter, "\n"))
      O_values <- O_values_next
      break
    }

    O_values <- O_values_next

    if (iter %% 10 == 0 || iter == 1) {
      cat(paste("\nIteration", iter, "\n"))
    }
  }

  cat("\n--- Final EM Results ---\n")
  cat("Final O_values (O_{j, i1}):\n")
  print(O_values)

  final_m_ij_df <- as.data.frame(m_ij_t)
  colnames(final_m_ij_df) <- paste0("m_est_", outcome_classes)
  final_df_with_m_ij_t <- cbind(data, final_m_ij_df)

  cat("\nEstimated refusals (m_ij) in final iteration (first 5 rows):\n")
  print(head(final_df_with_m_ij_t, 5))

  return(list(
    final_O_values = O_values,
    final_m_ij = final_m_ij_df,
    iterations = iter,
    final_data = final_df_with_m_ij_t
  ))
}
