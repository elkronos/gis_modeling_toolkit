library(testthat)
library(spatialkit)

# CRAN runs the files listed in cran-tests.txt, which keeps the check inside
# its time budget. With NOT_CRAN=true, as in CI and devtools::test(), every
# test file runs.
filter <- NULL
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  files <- readLines("cran-tests.txt")
  files <- files[nzchar(files) & !startsWith(files, "#")]
  filter <- paste0("^(", paste(files, collapse = "|"), ")$")
}

test_check("spatialkit", filter = filter)
