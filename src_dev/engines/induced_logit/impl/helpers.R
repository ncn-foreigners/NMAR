#' Shared induced-logit helpers
#'
#' @keywords internal
#' @noRd
IL_COL_R <- "..nmar_r.."
IL_COL_MU_HAT <- "..nmar_mu_hat.."

il_try_outcome_name <- function(formula) {
  tryCatch(
    {
      f <- tryCatch(Formula::as.Formula(formula), error = function(e) NULL)
      if (is.null(f)) return(NA_character_)
      lhs <- f[[2L]]
      if (is.symbol(lhs)) return(as.character(lhs))
      NA_character_
    },
    error = function(e) NA_character_
  )
}

#' @keywords internal
#' @noRd
il_try_n_respondents <- function(data, outcome_name) {
  if (!is.character(outcome_name) || length(outcome_name) != 1L || is.na(outcome_name)) return(NA_integer_)
  if (!is.data.frame(data) || !outcome_name %in% names(data)) return(NA_integer_)
  y <- data[[outcome_name]]
  if (is.null(y)) return(NA_integer_)
  as.integer(sum(!is.na(y)))
}

#' @keywords internal
#' @noRd
il_capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}
