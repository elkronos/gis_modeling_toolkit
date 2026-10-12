# tests/testthat/test-fix-variogram-directional.R
# ---------------------------------------------------------------------------
# The largest directional range never stands in for the all-pairs fit.
#
# When the all-pairs fit was singular or had not converged and two or more
# directions had fitted, estimate_sac_range() returned the largest
# directional range as the estimate (anisotropy_used = TRUE) and logged
# "notable anisotropy detected ... Using the maximum".  On ISOTROPIC fields
# (150 sites, effective range 30 km) that was 6 of 86 identified ranges, at a
# median 3.9 times the truth against 1.4 for the all-pairs answers, each from
# a quarter of the point pairs.
# ---------------------------------------------------------------------------

# Isotropic exponential fields of the review's reproduction: 150 sites in a
# 500 x 320 km box, sill 1, noise 0.2, effective range 30 km.
fvd_iso <- function(rep, eff = 30e3) {
  set.seed(90000L + rep + eff / 10)
  n <- 150L
  x <- stats::runif(n, 0, 500e3); y <- stats::runif(n, 0, 320e3)
  D <- as.matrix(stats::dist(cbind(x, y)))
  z <- drop(crossprod(chol(exp(-3 * D / eff) + diag(1e-8, n)), stats::rnorm(n))) +
    stats::rnorm(n, 0, sqrt(0.2))
  sf::st_as_sf(data.frame(x = 3e5 + x, y = 4.5e6 + y, z = z), coords = c("x", "y"),
               crs = 32615)
}
# Pairs of points no farther apart than `cutoff`: what an all-pairs variogram
# out to that lag is computed from.
fvd_pairs_within <- function(pts, cutoff) {
  d <- stats::dist(sf::st_coordinates(pts))
  sum(d <= cutoff)
}
# Draws on which the directional maximum used to be returned (12: 162.7 km,
# 17: 271.7 km, 55: 196.5 km, 5: 63.2 km, 34: 68.6 km, 65: 27.0 km) or was
# selected and then refused (50, 94).
FVD_REPS <- c(5, 12, 17, 34, 50, 55, 65, 94)

test_that("an isotropic field never gets the directional maximum, at the default", {
  skip_if_not_installed("gstat")
  for (rep in FVD_REPS) {
    pts <- fvd_iso(rep)
    lines <- capture_spatialkit_log(r <- estimate_sac_range(pts, "z"))
    info <- paste("rep", rep)
    expect_identical(attr(r, "anisotropy_used"), FALSE, info = info)
    expect_false(log_has(lines, "notable anisotropy"), info = info)
    expect_false(log_has(lines, "Using the maximum"), info = info)
    # The variogram attached is the all-pairs one, not a quarter of the pairs.
    expect_equal(sum(attr(r, "variogram")$np), fvd_pairs_within(pts, attr(r, "cutoff_dist")),
                 info = info)
    # The directional sweep is still reported, as the diagnostic it is.
    expect_identical(names(attr(r, "directional")), c("0", "45", "90", "135"))
    expect_identical(names(attr(r, "directional_status")), c("0", "45", "90", "135"))
    expect_identical(names(attr(r, "directional_fitted")), c("0", "45", "90", "135"))
    # An identified range is the all-pairs model's own.
    if (is.finite(r))
      expect_equal(as.numeric(r), 3 * max(attr(r, "variogram_model")$range), info = info)
  }
})

test_that("nor with the structure test off: the all-pairs fit decides, alone", {
  skip_if_not_installed("gstat")
  # structure_alpha = 1 takes away the one refusal the old code did not have,
  # so the same layers reach the fit, and the fit is what used to be replaced.
  for (rep in FVD_REPS) {
    pts <- fvd_iso(rep)
    lines <- capture_spatialkit_log(
      r <- estimate_sac_range(pts, "z", structure_alpha = 1))
    info <- paste("rep", rep)
    expect_identical(attr(r, "anisotropy_used"), FALSE, info = info)
    expect_false(log_has(lines, "Using the maximum"), info = info)
    expect_equal(sum(attr(r, "variogram")$np), fvd_pairs_within(pts, attr(r, "cutoff_dist")),
                 info = info)
    if (is.finite(r)) {
      expect_equal(as.numeric(r), 3 * max(attr(r, "variogram_model")$range), info = info)
    } else {
      expect_false(is.null(attr(r, "rejected_reason")), info = info)
    }
  }
  # The three draws that returned 163, 272 and 197 km for a 30 km field.
  far <- vapply(c(12, 17, 55), function(rep)
    as.numeric(suppressWarnings(estimate_sac_range(fvd_iso(rep), "z", structure_alpha = 1))),
    numeric(1))
  expect_true(all(is.na(far) | far < 100e3))
})

test_that("with both lag settings, and for a refusal, anisotropy_used is FALSE", {
  skip_if_not_installed("gstat")
  pts <- fvd_iso(12)
  for (l in c("adaptive", "equal")) {
    r <- suppressWarnings(estimate_sac_range(pts, "z", lags = l, structure_alpha = 1))
    expect_identical(attr(r, "anisotropy_used"), FALSE)
    rr <- suppressWarnings(estimate_sac_range(pts, "z", lags = l, range_frac = 1e-6,
                                              structure_alpha = 1))
    expect_true(is.na(rr))
    expect_identical(attr(rr, "anisotropy_used"), FALSE)
  }
})

test_that("an all-pairs range past the lags is refused whatever the directions found", {
  skip_if_not_installed("gstat")
  # An exponential field under an east-west trend: the pooled variogram rises
  # past the fitted lags.  Directions across the slope can reach a sill; they
  # are reported and nothing is built from them.
  set.seed(11)
  n <- 300
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 50) + diag(0.2 + 1e-8, n))) %*% rnorm(n)) + 4 * x / 1000
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
  lines <- capture_spatialkit_log(r <- suppressWarnings(estimate_sac_range(pts, "z")))
  expect_true(is.na(r))
  expect_identical(attr(r, "rejected_reason"), "fitted range exceeds the largest lag fitted")
  expect_gt(attr(r, "rejected_range"), attr(r, "cutoff_dist"))
  expect_identical(attr(r, "anisotropy_used"), FALSE)
  st <- attr(r, "directional_status")
  expect_true(all(st %in% c("ok", "over_cutoff", "not_converged", "no_fit")))
  expect_equal(sum(is.finite(attr(r, "directional"))), sum(st == "ok"))
  expect_true(log_has(lines, "exceeds the largest lag"))
  expect_false(log_has(lines, "Using the maximum"))
})
