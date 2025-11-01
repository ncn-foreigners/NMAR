## Ensure package is attached when tests are run outside devtools::test()
## testthat::test_dir() does not load the package automatically.
## This is a no-op under devtools::test() where NMAR is already loaded.
if (!("package:NMAR" %in% search())) {
  suppressPackageStartupMessages(library(NMAR))
}
