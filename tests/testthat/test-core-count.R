# ===========================================================================
# Worker-count resolution for the parallel paths.
#
# parallel::detectCores() returns NA on some platforms, and NA arithmetic
# propagates rather than erroring, so an unsanitised count reaches
# mclapply() as NA. See test-cv-parallel.R for fold-level reproducibility.
# ===========================================================================

test_that(".sanitize_core_count collapses NA/invalid input to the fallback", {
  f <- spatialkit:::.sanitize_core_count

  expect_identical(f(NA), 1L)
  expect_identical(f(NA_integer_ - 1L), 1L)   # detectCores() NA arithmetic
  expect_identical(f(0), 1L)
  expect_identical(f(-3), 1L)
  expect_identical(f(NULL), 1L)
  expect_identical(f("not a number"), 1L)
  expect_identical(f(7), 7L)
  expect_identical(f(3.9), 3L)
  expect_identical(f(NA, fallback = 2L), 2L)
})


test_that(".resolve_n_cores returns sane values for the sequential paths", {
  f <- spatialkit:::.resolve_n_cores

  expect_identical(f(parallel = FALSE), 1L)
  expect_identical(f(parallel = 1), 1L)
  expect_identical(f(parallel = NA), 1L)

  skip_on_os("windows")
  expect_identical(f(parallel = 3), 3L)
  expect_identical(f(parallel = FALSE, n_cores = 2), 2L)
})
