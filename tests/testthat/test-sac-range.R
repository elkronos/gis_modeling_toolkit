# tests/testthat/test-sac-range.R
# ---------------------------------------------------------------------------
# estimate_sac_range(): the range guard, and the fit it now returns.
#
# gstat::fit.variogram() returns a finite number even when the empirical
# variogram never reaches a sill.  The range is then unidentified and the value
# is a fitting artefact -- on continental-extent data it came back at 1.8x the
# diameter of the data itself.  Returning that silently is worse than returning
# NA, because make_folds(auto_range = TRUE) sizes blocks from it.
# ---------------------------------------------------------------------------

# A Gaussian random field with an EXACTLY KNOWN exponential range, simulated by
# Cholesky factorisation of the covariance matrix.
#
# An earlier version kernel-smoothed white noise at scattered anchors.  That
# produced a field whose correlation length was far longer than the nominal
# bandwidth -- gstat fitted a range of 1352 on a 1401-unit extent, the guard
# (correctly) rejected it, and every success-path test silently skipped.  A
# simulated field is worth the extra lines: the true range is known, so the
# tests can assert that it is recovered rather than merely that something
# finite came back.
sac_test_field <- function(n = 250, extent = 1000, true_range = 80, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, extent); y <- runif(n, 0, extent)
  d <- as.matrix(stats::dist(cbind(x, y)))
  # Exponential covariance exp(-h/a); gstat reports `a`, and .fit_vgm_range()
  # converts to the 95% practical range 3a.  So a = true_range / 3.
  C <- exp(-d / (true_range / 3))
  diag(C) <- diag(C) + 1e-4                     # nugget for numerical PD
  z <- as.numeric(t(chol(C)) %*% rnorm(n)) + rnorm(n, 0, 0.05)
  sf::st_as_sf(data.frame(x = x, y = y, z = z),
               coords = c("x", "y"), crs = 3857)
}

TRUE_RANGE <- 80

test_that("range_frac rejects a range the extent cannot support", {
  # Deterministic regardless of how gstat fits: range_frac = 1e-6 rejects any
  # positive range, so this exercises the guard itself rather than the fit.
  skip_if_not_installed("gstat")
  pts <- sac_test_field()

  expect_true(is.na(estimate_sac_range(pts, "z", range_frac = 1e-6, seed = 1)))

  lines <- capture_spatialkit_log(
    estimate_sac_range(pts, "z", range_frac = 1e-6, seed = 1)
  )
  expect_true(log_has(lines, "exceeds the largest"))
  expect_true(log_has(lines, "unidentified"))
})

test_that("a rejected range is NA but keeps the variogram that rejected it", {
  # The VALUE must stay NA -- make_folds() guards with is.finite() and must
  # never size blocks from an unidentified range.  But the variogram was
  # already computed, and discarding it left plot(type = "variogram") unable to
  # draw exactly the case worth looking at: a curve with no sill.  Attributes
  # do not affect is.na()/is.finite(), so the guards behave identically.
  skip_if_not_installed("gstat")
  r <- estimate_sac_range(sac_test_field(), "z", range_frac = 1e-6, seed = 1)
  expect_true(is.na(r))
  expect_false(is.finite(r))
  expect_identical(as.numeric(r), NA_real_)
  expect_null(attr(r, "class"))            # still a bare numeric, not classed

  # the diagnostics that justify the rejection
  expect_false(is.null(attr(r, "variogram")))
  expect_s3_class(attr(r, "variogram"), "data.frame")
  expect_true(is.finite(attr(r, "rejected_range")))
  expect_true(is.finite(attr(r, "cutoff_dist")))
  expect_match(attr(r, "rejected_reason"), "exceeds the largest lag")

  # and the guard downstream still sees it as unusable
  expect_false(isTRUE(is.finite(r)))
})

test_that("a supportable range is returned and carries its fit", {
  skip_if_not_installed("gstat")
  pts <- sac_test_field()
  r <- estimate_sac_range(pts, "z", seed = 1)
  # No skip: the field is simulated deterministically with a fixed seed, so a
  # failed fit is a defect to surface, not a platform quirk to step around.
  expect_false(is.na(r))

  expect_s3_class(r, "sac_range")
  expect_gt(as.numeric(r), 0)
  expect_true(is.finite(attr(r, "max_dist")))
  expect_true(is.finite(attr(r, "cutoff_dist")))
  expect_lte(as.numeric(r), attr(r, "cutoff_dist"))      # the default guard
  expect_false(is.null(attr(r, "directional")))
  expect_length(attr(r, "directional"), 2L)

  # The field was simulated with a known exponential range, so the estimate
  # should land near it.  A wide band -- variogram estimation on 250 irregular
  # points is noisy -- but tight enough to catch a systematic error such as the
  # Exp 3a conversion being dropped.
  expect_gt(as.numeric(r), TRUE_RANGE / 3)
  expect_lt(as.numeric(r), TRUE_RANGE * 3)
})

test_that("the returned range behaves as an ordinary number", {
  # Callers do arithmetic and comparisons on this; the class must be transparent.
  skip_if_not_installed("gstat")
  r <- estimate_sac_range(sac_test_field(), "z", seed = 1)
  expect_false(is.na(r))

  expect_true(is.numeric(r))
  expect_equal(r * 2, as.numeric(r) * 2, ignore_attr = TRUE)
  expect_true(is.finite(r + 1))
  expect_true(r > 0)
  expect_true(is.finite(max(as.numeric(r), 1)))
})

test_that("print.sac_range leads with a bare number", {
  skip_if_not_installed("gstat")
  r <- estimate_sac_range(sac_test_field(), "z", seed = 1)
  expect_false(is.na(r))
  expect_output(print(r), "^[0-9]")
})

test_that("make_folds(auto_range) uses an identified range as the block size", {
  # Happy path: the field has real structure, so the range is usable and gets
  # applied.  (This test was previously named for the fallback and asserted
  # nothing about it -- with a well-identified field it never fell back.)
  skip_if_not_installed("gstat")
  pts <- sac_test_field()

  f <- make_folds(pts, k = 4, method = "block_kfold", auto_range = TRUE,
                  response_var = "z", seed = 1)
  expect_equal(f$method, "block_kfold")
  expect_gte(length(f$folds), 2L)
  expect_gt(f$params$grid_nx * f$params$grid_ny, 1L)
  expect_true(is.finite(f$params$sac_range))
  expect_gt(f$params$sac_range, 0)
})

test_that("make_folds(auto_range) falls back when the range is unidentified", {
  # The failure the guard exists to prevent: an unusable range must not size a
  # single block over the whole extent.  range_frac forces rejection
  # deterministically, rather than hoping a contrived field fails to fit.
  skip_if_not_installed("gstat")
  pts <- sac_test_field()

  lines <- capture_spatialkit_log(
    f <- make_folds(pts, k = 4, method = "block_kfold", auto_range = TRUE,
                    range_frac = 1e-6, response_var = "z", seed = 1)
  )
  expect_true(log_has(lines, "falling back to geometric blocks"))

  # ... and the grid is still a real grid, not one block.
  expect_gt(f$params$grid_nx * f$params$grid_ny, 1L)
  expect_gte(length(f$folds), 2L)
})
