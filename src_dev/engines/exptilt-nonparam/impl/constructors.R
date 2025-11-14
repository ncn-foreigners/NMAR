#' @keywords internal
new_nmar_result_exptilt_nonparam <- function(
    n_y_x_matrix,
    m_x_vec,
    p_hat_y_given_x_matrix,
    x1_cols,
    x2_cols,
    outcome_cols,
    engine_params,
    fit_stats,
    loss_value = NULL,
    data_to_return,
    ...

) {
  result <- list(
    n_y_x_matrix = n_y_x_matrix,
    m_x_vec = m_x_vec,
    p_hat_y_given_x_matrix = p_hat_y_given_x_matrix,
    x1_cols = x1_cols,
    x2_cols = x2_cols,
    outcome_cols = outcome_cols,
    engine_params = engine_params,
    fit_stats = fit_stats,
    loss_value = loss_value,
    data_final = data_to_return

  )
  class(result) <- c("nmar_result_exptilt_nonparam", "nmar_result")
  return(result)
}
