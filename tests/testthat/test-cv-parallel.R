# tests/testthat/test-cv-parallel.R
# -------------------------------------------------------------------
# Tests for parallel fold-fitting support in CV functions
# -------------------------------------------------------------------

test_that(".resolve_n_cores returns 1L for FALSE / NULL", {
  expect_equal(spatialkit:::.resolve_n_cores(FALSE), 1L)
  expect_equal(spatialkit:::.resolve_n_cores(FALSE, NULL), 1L)
})

test_that(".resolve_n_cores respects explicit integer", {
  skip_on_os("windows")
  expect_equal(spatialkit:::.resolve_n_cores(parallel = 2L), 2L)
  expect_equal(spatialkit:::.resolve_n_cores(parallel = FALSE, n_cores = 3L), 3L)
})

test_that(".resolve_n_cores auto-detects when TRUE", {
  skip_on_os("windows")
  cores <- spatialkit:::.resolve_n_cores(TRUE)
  expect_true(is.integer(cores))
  expect_true(cores >= 1L)
})

test_that(".resolve_n_cores falls back to 1L on Windows", {
  skip_on_os(c("mac", "linux", "solaris"))
  expect_message(
    cores <- spatialkit:::.resolve_n_cores(TRUE),
    "not available on Windows"
  )
  expect_equal(cores, 1L)
})

test_that("cv_gwr accepts parallel argument without error", {
  # This test just verifies the parameter is accepted in the signature;

  # actual parallel execution requires spatial model backends.
  expect_true("parallel" %in% names(formals(cv_gwr)))
})

test_that("cv_bayes accepts parallel argument without error", {
  expect_true("parallel" %in% names(formals(cv_bayes)))
})

test_that("cv_spatial accepts parallel argument without error", {
  expect_true("parallel" %in% names(formals(cv_spatial)))
})
