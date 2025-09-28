# This development-only guard enforces engine boundaries.
# It is skipped on CRAN and when developing from the built R/ tree only.

test_that("engines do not reference EL internals", {
  skip_on_cran()
  if (!dir.exists("src_dev")) skip("boundary test runs only in src_dev context")
# Files outside EL engine impl
  all_r_files <- list.files("src_dev", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  offenders <- grep("/engines/el/", all_r_files, invert = TRUE, value = TRUE)
# Look for calls to el_* helpers
  bad <- character(0)
  for (f in offenders) {
    txt <- readLines(f, warn = FALSE)
    if (any(grepl("\\bel_\\w+\\s*\\(", txt))) bad <- c(bad, f)
  }
  expect(length(bad) == 0, failure_message = paste("Cross-engine references to EL internals found in:\n", paste(bad, collapse = "\n")))
})
