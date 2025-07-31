#' Print method for nmar_result
#'
#' @param x nmar_result object
#' @param ... Additional parameters

print.nmar_result <- function(x, ...) {
  validate_nmar_result(x)
  cat("NMAR Model Results\n")
  cat("-----------------\n")
  cat("Theta parameters:\n")
  print(x$theta)
  cat("\nEstimated mean:", x$est_mean, "\n")

  # Handle vector loss_value
  if (length(x$loss_value) == 1) {
    cat("Loss value:", x$loss_value, "\n")
  } else {
    cat("Loss values:\n")
    print(x$loss_value)
  }

  invisible(x)
}

#' Summary method for nmar_result
#'
#' @param object nmar_result object
#' @param ... Additional parameters

summary.nmar_result <- function(object, ...) {
  validate_nmar_result(object)

  # Handle vector loss_value
  if (length(object$loss_value) > 1) {
    convergence <- ifelse(all(object$loss_value < 1e-5),
                          "Achieved", "Not achieved")
    loss_value <- min(object$loss_value) # Show smallest loss
  } else {
    convergence <- ifelse(object$loss_value < 1e-5,
                          "Achieved", "Not achieved")
    loss_value <- object$loss_value
  }

  structure(
    list(
      theta_summary = summary(object$theta),
      est_mean = object$est_mean,
      loss_value = loss_value,
      convergence = convergence,
      all_loss_values = if(length(object$loss_value) > 1) object$loss_value else NULL
    ),
    class = "summary_nmar_result"
  )
}

#' Print method for summary.nmar_result
#'
#' @param x summary_nmar_result object
#' @param ... Additional parameters

print.summary_nmar_result <- function(x, ...) {
  cat("NMAR Model Summary\n")
  cat("=================\n")
  cat("Theta parameters summary:\n")
  print(x$theta_summary)
  cat("\nEstimated mean:", x$est_mean, "\n")
  cat("Final loss value:", x$loss_value, "\n")
  if (!is.null(x$all_loss_values)) {
    cat("All loss values:", paste(x$all_loss_values, collapse = ", "), "\n")
  }
  cat("Convergence:", x$convergence, "\n")
}
