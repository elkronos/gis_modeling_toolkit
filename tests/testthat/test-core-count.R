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
  expect_identical(f(parallel = FALSE, n_cores = NULL), 1L)
  expect_identical(f(parallel = 1), 1L)
  expect_identical(f(parallel = NA), 1L)
  expect_identical(f(parallel = NULL), 1L)

  skip_on_os("windows")
  expect_identical(f(parallel = 3), 3L)
  expect_identical(f(parallel = 2L), 2L)
  expect_identical(f(parallel = FALSE, n_cores = 2), 2L)
  # n_cores wins over parallel, and is sanitised on the way through.
  expect_identical(f(parallel = 8L, n_cores = 3L), 3L)
  expect_identical(f(parallel = TRUE, n_cores = NA), 1L)
})


test_that(".resolve_n_cores auto-detection returns a usable worker count", {
  # parallel = TRUE reads detectCores(), which is machine-dependent, so the
  # assertable contract is the type and the floor -- an NA reaching mclapply()
  # is the failure this guards.
  skip_on_os("windows")
  cores <- spatialkit:::.resolve_n_cores(TRUE)
  expect_type(cores, "integer")
  expect_length(cores, 1L)
  expect_false(is.na(cores))
  expect_gte(cores, 1L)
  # Never more than the machine has; detectCores() - 1 is the documented rule.
  detected <- parallel::detectCores(logical = FALSE)
  if (!is.na(detected) && detected > 1L) expect_lte(cores, detected - 1L)
})
