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


# ---------------------------------------------------------------------------
# Fold-level reproducibility.
#
# These replace three earlier tests that asserted only
#   expect_true("parallel" %in% names(formals(cv_gwr)))
# which cannot fail for any implementation, correct or not.
#
# All four use boot_fit_fn() from helper-bootlm.R, which resamples inside
# every fold.  A deterministic learner would make these pass regardless of
# whether fold seeding works.
# ---------------------------------------------------------------------------

test_that("parallel and sequential CV agree under a fixed seed", {
  skip_on_os("windows")
  skip_on_cran()

  pts <- make_cv_test_points()
  seq_res <- cv_spatial(pts, "z", "w", fit_fn = boot_fit_fn, k = 4,
                        seed = 99, parallel = FALSE)
  par_res <- cv_spatial(pts, "z", "w", fit_fn = boot_fit_fn, k = 4,
                        seed = 99, parallel = 2L)

  expect_equal(seq_res$overall, par_res$overall)
  expect_equal(seq_res$fold_metrics$RMSE, par_res$fold_metrics$RMSE)
  expect_equal(seq_res$predictions$yhat, par_res$predictions$yhat)
})

test_that("parallel CV is reproducible across repeated runs", {
  skip_on_os("windows")
  skip_on_cran()

  pts <- make_cv_test_points()
  a <- cv_spatial(pts, "z", "w", fit_fn = boot_fit_fn, k = 4,
                  seed = 7, parallel = 2L)
  b <- cv_spatial(pts, "z", "w", fit_fn = boot_fit_fn, k = 4,
                  seed = 7, parallel = 2L)
  expect_equal(a$overall, b$overall)
})

test_that("different seeds give different results", {
  # Guards the two tests above.  Without this, a refactor that accidentally
  # made boot_fit_fn deterministic would turn them green for the wrong reason.
  skip_on_cran()

  pts <- make_cv_test_points()
  a <- cv_spatial(pts, "z", "w", fit_fn = boot_fit_fn, k = 4,
                  seed = 1, parallel = FALSE)
  b <- cv_spatial(pts, "z", "w", fit_fn = boot_fit_fn, k = 4,
                  seed = 2, parallel = FALSE)
  expect_false(isTRUE(all.equal(a$overall$RMSE, b$overall$RMSE)))
})

test_that("CV does not disturb the caller's RNG state", {
  # cran-comments.md tells reviewers that RNG state is saved and restored
  # around all seeded operations.  This is the automated form of that claim.
  skip_on_cran()

  pts <- make_cv_test_points()
  set.seed(5); before <- runif(1)
  set.seed(5); invisible(cv_spatial(pts, "z", "w", fit_fn = boot_fit_fn,
                                    k = 4, seed = 99, parallel = FALSE))
  after <- runif(1)
  expect_equal(before, after)
})
