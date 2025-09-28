#!/usr/bin/env Rscript


if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Devtools required to run check()!")
}

cat("--- Running devtools::check() before git push... ---\n")

check_result <- tryCatch({
  devtools::check(
    document = FALSE,
    build_args = c("--no-build-vignettes", "--no-manual"),
    error_on = "error"
  )
}, error = function(e) {
  cat("\n======================================================\n")
  cat("ERR: devtools::check() failed.")
  cat("\nPush blocked. Fix your code\n")
  cat("======================================================\n")
  return(1)
})


if (is.list(check_result) && all(check_result$errors == 0)) {
  cat("Check() run succesfully. \u2705\n")
  quit(save = "no", status = 0)
} else if (is.numeric(check_result) && check_result == 1) {
  quit(save = "no", status = 1)
} else {
  cat("Check finished. Some warns but not Errs. \u2714\n")
  quit(save = "no", status = 0)
}
