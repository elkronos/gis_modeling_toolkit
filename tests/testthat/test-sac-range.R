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
  # A real measurement nugget (sd 0.3 on a unit-sill field), not a token one.
  # Nugget-free, the weighted least-squares fit sits against the nugget's zero
  # bound and gstat's exponential fit fails to converge from every start
  # while its spherical fit is flagged singular or not depending on
  # floating-point details: the same field fitted on Linux and came back
  # singular on an arm64 Mac.  With a nugget every start converges to one
  # optimum on both.
  z <- as.numeric(t(chol(C)) %*% rnorm(n)) + rnorm(n, 0, 0.3)
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
  # Classed exactly as the success return, so print.sac_range() fires for a
  # rejected range too.  Before this, printing one dumped the empirical
  # variogram data.frame and the fitted gstat model as raw attributes.
  expect_s3_class(r, "sac_range")
  expect_output(print(r), "^NA")
  printed <- paste(utils::capture.output(print(r)), collapse = "\n")
  expect_false(grepl("np|dist|gamma|psill", printed))   # no variogram dump

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
  # FOUR directions, named by azimuth.  A 0/90 sweep at +/-22.5 degrees covers
  # only 90 of the 180 distinct azimuths, so a field oriented near 45 or 135
  # degrees fell into neither window and had its range halved (measured: 151
  # and 147 against a true 300, versus 255 and 249 at 0 and 90).  c(0, 45, 90,
  # 135) tiles all of them.  An entry is NA when that direction's variogram
  # never reached a sill and was excluded from the maximum.
  d <- attr(r, "directional")
  expect_false(is.null(d))
  expect_length(d, 4L)
  expect_identical(names(d), c("0", "45", "90", "135"))
  # Each direction sees about a quarter of the pairs, so on 250 points some
  # of the four are expected to come back NA; that is a diagnostic, not a
  # failure, since nothing is built from them any more.
  expect_true(all(is.na(d) | d > 0))

  # The returned number is the OMNIDIRECTIONAL fit whenever that fit is
  # usable; the directional ranges are a diagnostic.  Each direction sees
  # about a quarter of the point pairs, so the maximum of four noisy estimates
  # is biased upward (on this isotropic field with a true range of 80 the
  # surviving direction reports 93 while the all-pairs fit lands on 56), and
  # the windows are fixed to the axes, so nothing built from them is rotation
  # invariant.  `anisotropy_used` is TRUE only when the all-pairs fit failed.
  expect_false(is.null(attr(r, "anisotropy_used")))
  expect_false(isTRUE(attr(r, "anisotropy_used")))
  expect_true(is.finite(attr(r, "anisotropy")) || sum(is.finite(d)) < 2L)

  # The field was simulated with a known exponential range, so the estimate
  # should land near it.  A wide band -- variogram estimation on 250 irregular
  # points is noisy -- but tight enough to catch a systematic error such as the
  # Exp 3a conversion being dropped.
  expect_gt(as.numeric(r), TRUE_RANGE / 3)
  expect_lt(as.numeric(r), TRUE_RANGE * 3)
})


test_that("the all-pairs fit is the estimate, not the widest direction", {
  # Regression: the estimate used to be max-of-four directions unconditionally,
  # which on an isotropic field returned the widest of the directions that
  # happened to fit -- three times the truth -- and sized blocks from it.  Each
  # direction sees about a quarter of the point pairs; the omnidirectional fit
  # sees all of them.
  skip_if_not_installed("gstat")
  r <- suppressWarnings(estimate_sac_range(sac_test_field(), "z", seed = 1))
  d <- attr(r, "directional")

  expect_false(isTRUE(attr(r, "anisotropy_used")))
  expect_true(any(is.finite(d)))
  expect_lt(as.numeric(r), max(d, na.rm = TRUE))   # strictly below the max
  expect_lt(abs(as.numeric(r) / TRUE_RANGE - 1), 0.5)
})

test_that("the variogram fit does not depend on gstat's default starting range", {
  # gstat starts the optimiser at a third of the longest lag.  For a field
  # whose range is a small fraction of the extent that is ten times too long,
  # and whether the iteration lands or collapses to a singular model depended
  # on floating-point details: the same 250-point field fitted on Linux and
  # came back singular on an arm64 Mac, where the estimate then fell through
  # to the directional maximum (315 against a true range of 80).  Several
  # starting ranges are tried now, and the best converged fit by gstat's own
  # criterion is kept.  Made deterministic: the default start (range = NA)
  # is forced singular, and the answer must be the one the other starts give.
  skip_if_not_installed("gstat")
  pts <- sac_test_field()
  plain <- estimate_sac_range(pts, "z", seed = 1)
  expect_false(is.na(plain))
  real_fit <- gstat::fit.variogram
  local_mocked_bindings(
    fit.variogram = function(object, model, ...) {
      if (anyNA(model$range)) {
        model$psill <- c(0, 0); model$range <- c(0, max(object$dist) / 3)
        attr(model, "singular") <- TRUE
        return(model)
      }
      real_fit(object, model, ...)
    },
    .package = "gstat")
  retried <- estimate_sac_range(pts, "z", seed = 1)
  expect_false(is.na(retried))
  expect_false(isTRUE(attr(retried, "anisotropy_used")))
  expect_equal(as.numeric(retried), as.numeric(plain), tolerance = 0.1)
  expect_lt(abs(as.numeric(retried) / TRUE_RANGE - 1), 0.5)
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

  # "Use it as the minimum block size" means the block size IS the range and
  # every block edge is at least that long.  Nothing related block_size or
  # the grid to sac_range, so blocks of half the range -- smaller than the
  # correlation length, the leakage blocked CV exists to prevent -- passed
  # (mutation testing, pass 6).
  expect_equal(f$params$block_size, f$params$sac_range)
  bb <- sf::st_bbox(pts)
  edge_x <- as.numeric(bb["xmax"] - bb["xmin"]) / f$params$grid_nx
  edge_y <- as.numeric(bb["ymax"] - bb["ymin"]) / f$params$grid_ny
  expect_gte(min(edge_x, edge_y), f$params$sac_range)
  # ... and the grid is the one that block size implies, not the geometric
  # block_multiplier * k default.
  expect_equal(f$params$grid_nx,
               max(1L, floor(as.numeric(bb["xmax"] - bb["xmin"]) / f$params$sac_range)))
  expect_equal(f$params$grid_ny,
               max(1L, floor(as.numeric(bb["ymax"] - bb["ymin"]) / f$params$sac_range)))
  # The range that sized the blocks is the one estimate_sac_range() reports.
  r <- estimate_sac_range(pts, "z", seed = 1)
  expect_equal(as.numeric(f$params$sac_range), as.numeric(r))
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


# ---------------------------------------------------------------------------
# Input validation and reproducibility
# ---------------------------------------------------------------------------

test_that("a factor response is refused rather than silently coded", {
  # as.numeric() on a factor returns LEVEL CODES -- an arbitrary integer
  # relabelling of the categories -- so a variogram fitted to them changed when
  # the levels were reordered (3700 against 2497 on the same data).
  skip_if_not_installed("gstat")
  pts <- sac_test_field()
  pts$grp <- factor(sample(c("a", "b", "c"), nrow(pts), replace = TRUE))
  expect_error(estimate_sac_range(pts, "grp"), "factor")
  pts$txt <- as.character(pts$grp)
  expect_error(estimate_sac_range(pts, "txt"), "numeric")
  # A missing column is named, not discovered downstream as "too few values".
  expect_error(estimate_sac_range(pts, "nope"), "not found")
  # Logical is fine: 0/1 is a well-defined variogram target.
  pts$flag <- pts$z > stats::median(pts$z)
  expect_false(is.null(estimate_sac_range(pts, "flag")))
})

test_that("the n_max subsample is reproducible and leaves the RNG alone", {
  # `seed` defaults to a constant.  Unseeded, the subsample made the returned
  # range differ between runs on identical input (19531 / 19589 / 19605) and
  # silently advanced the caller's stream -- and make_folds(auto_range = TRUE)
  # sizes its blocks from that number.
  skip_if_not_installed("gstat")
  set.seed(99)
  n   <- 400
  pts <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
    coords = c("x", "y"), crs = 32632)

  a <- suppressWarnings(estimate_sac_range(pts, "z", n_max = 150L))
  b <- suppressWarnings(estimate_sac_range(pts, "z", n_max = 150L))
  expect_equal(as.numeric(a), as.numeric(b))

  # The caller's stream is untouched across the call.
  set.seed(1); before <- runif(3)
  set.seed(1); invisible(suppressWarnings(estimate_sac_range(pts, "z", n_max = 150L)))
  after <- runif(3)
  expect_equal(before, after)

  # And a different seed is still available for a sensitivity check.
  d <- suppressWarnings(estimate_sac_range(pts, "z", n_max = 150L, seed = 999L))
  expect_true(is.na(d) || is.numeric(as.numeric(d)))
})

test_that("all four azimuths are represented in the directional attribute", {
  skip_if_not_installed("gstat")
  r <- estimate_sac_range(sac_test_field(), "z", seed = 1)
  d <- attr(r, "directional")
  expect_identical(names(d), c("0", "45", "90", "135"))
  # 0 and 90 alone at +/-22.5 degrees cover only half the azimuth circle.
  covered <- function(az, tol) {
    vapply(0:179, function(th) any(abs(((th - az + 90) %% 180) - 90) <= tol),
           logical(1))
  }
  expect_equal(sum(covered(c(0, 90), 22.5)), 90L)
  expect_equal(sum(covered(c(0, 45, 90, 135), 22.5)), 180L)
})


test_that("a non-converged variogram fit is refused but stays inspectable", {
  # Two properties that pull against each other, and both matter.
  #
  # (1) The RANGE must be refused. gstat signals non-convergence with a warning
  #     and then returns anyway, so the number it reports is wherever the
  #     optimiser stopped rather than a fitted parameter, and
  #     make_folds(auto_range = TRUE) would size blocks from it.
  # (2) The VARIOGRAM must survive. A curve that never reaches a sill is
  #     exactly the case worth looking at, so plot(type = "variogram") has to
  #     keep working -- discarding the fit made it error with "the residual
  #     variogram could not be fitted".
  #
  # And no bare gstat warning should reach the user: with the sweep at four
  # azimuths each variogram gets about half the pairs, so a direction failing
  # to converge is routine and the caller already handles it.
  skip_if_not_installed("gstat")
  set.seed(21)
  n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  # A pure linear trend: the variogram rises monotonically and never sills.
  pts <- sf::st_as_sf(
    data.frame(x = x, y = y, z = 0.02 * x + 0.01 * y + rnorm(n, 0, 0.5)),
    coords = c("x", "y"), crs = 3857)

  expect_no_warning(r <- estimate_sac_range(pts, "z", seed = 1))

  expect_true(is.na(r))                       # (1) refused
  expect_s3_class(r, "sac_range")
  expect_false(is.null(attr(r, "variogram")))          # (2) still inspectable
  expect_false(is.null(attr(r, "variogram_model")))
  expect_true(is.finite(attr(r, "rejected_range")))
  expect_match(attr(r, "rejected_reason"),
               "did not converge|exceeds the largest lag")
  # It prints as a bare NA rather than dumping its attributes.
  expect_output(print(r), "NA")
})

test_that("estimate_sac_range never emits a raw gstat warning", {
  # The four-azimuth sweep makes a failed directional fit ordinary. Whatever
  # the data, the failure is handled internally -- the direction is excluded,
  # or the whole range is refused -- and never surfaces as gstat's own warning.
  skip_if_not_installed("gstat")
  for (s in 1:6) {
    set.seed(100 + s)
    n <- 120
    pts <- sf::st_as_sf(
      data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
      coords = c("x", "y"), crs = 32632)
    expect_no_warning(estimate_sac_range(pts, "z", seed = 1),
                      message = paste("seed", s))
  }
})


test_that("estimate_sac_range records the CRS its range is measured in", {
  # The range is a length in this CRS; summarize_by_cell(deff = "variogram")
  # transforms its points to it before evaluating within-cell distances.
  skip_if_not_installed("gstat")
  r <- estimate_sac_range(sac_test_field(), "z", seed = 1)
  expect_s3_class(attr(r, "crs"), "crs")
  expect_false(is.na(attr(r, "crs")))
  # Geographic input is projected first, so the recorded CRS is projected too.
  geo <- sf::st_transform(sac_test_field(), 4326)
  rg  <- suppressWarnings(estimate_sac_range(geo, "z", seed = 1))
  expect_false(isTRUE(sf::st_is_longlat(attr(rg, "crs"))))
})


# ---------------------------------------------------------------------------
# Detrending: OLS residual variograms are biased toward a shorter range (Lark,
# Cullis & Welham 2006); detrend = "reml" fits trend and covariance together.
# ---------------------------------------------------------------------------

# A field with a quadratic trend surface in the coordinates -- the case where
# the OLS bias is large (measured 0.75 of the oracle range; see
# ?estimate_sac_range) -- with the trend terms as predictor columns.
sac_trend_field <- function(n = 300, seed = 4001) {
  set.seed(seed)
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  S  <- as.numeric(t(chol(exp(-D / 100) + diag(1e-8, n))) %*% rnorm(n))
  X1 <- xy[, 1] / 1000; X2 <- xy[, 2] / 1000
  d <- data.frame(x = xy[, 1], y = xy[, 2], X1 = X1, X2 = X2,
                  X3 = X1^2, X4 = X2^2, X5 = X1 * X2)
  d$z <- 1 + 2 * X1 + X2 + 1.5 * d$X3 - 1.5 * d$X4 + 2 * d$X5 + S +
    rnorm(n, sd = sqrt(0.2))
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
TREND_VARS <- c("X1", "X2", "X3", "X4", "X5")

test_that("detrend defaults to OLS and records which method was used", {
  skip_if_not_installed("gstat")
  fld <- sac_trend_field()
  raw <- estimate_sac_range(fld, "z")
  expect_true(is.na(attr(raw, "detrend_method")))
  expect_false(attr(raw, "detrended"))
  ols <- estimate_sac_range(fld, "z", TREND_VARS)
  expect_identical(attr(ols, "detrend_method"), "ols")
  expect_true(attr(ols, "detrended"))
  expect_null(attr(ols, "reml"))
  expect_error(estimate_sac_range(fld, "z", TREND_VARS, detrend = "gls"),
               "'arg' should be one of")
  expect_error(estimate_sac_range(fld, "z", TREND_VARS, reml_max_n = 10),
               "at least 30")
})

test_that("detrend = 'reml' returns the REML range with its fit attached", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  fld <- sac_trend_field()
  r <- estimate_sac_range(fld, "z", TREND_VARS, detrend = "reml")
  expect_s3_class(r, "sac_range")
  expect_true(is.finite(r))
  expect_identical(attr(r, "detrend_method"), "reml")
  expect_true(attr(r, "detrended"))
  info <- attr(r, "reml")
  expect_type(info, "list")
  expect_identical(info$n_used, nrow(fld))
  expect_false(info$subsampled)
  expect_true(info$nugget_prop >= 0 && info$nugget_prop <= 1)
  expect_true(info$sigma2 > 0)
  # The model behind the number is the REML model, in gstat's layout, so
  # plot() and summarize_by_cell(deff = "variogram") read it as usual.
  vm <- attr(r, "variogram_model")
  expect_identical(as.character(vm$model), c("Nug", "Exp"))
  expect_equal(sum(vm$psill), info$sigma2, tolerance = 1e-8)
  expect_equal(as.numeric(r), 3 * vm$range[2], tolerance = 1e-8)
  expect_equal(sac_nugget(r), vm$psill[1], tolerance = 1e-12)
  # The empirical variogram of the REML residuals travels with it, and the
  # directional sweep on those residuals is still reported.
  expect_s3_class(attr(r, "variogram"), "data.frame")
  expect_length(attr(r, "directional"), 4L)
  expect_false(attr(r, "anisotropy_used"))
  expect_s3_class(plot(r), "ggplot")
})

test_that("REML detrending lands closer to the oracle than OLS on a trend surface", {
  # A pin of the direction on one draw, not a proof: over 40 draws the
  # median ratios to the oracle were OLS 0.75 and REML 1.06 on this design
  # (?estimate_sac_range).  This seed's draw is one where the ordering holds;
  # a change that flips it here is worth looking at.
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  fld <- sac_trend_field(seed = 4001)
  ols  <- as.numeric(estimate_sac_range(fld, "z", TREND_VARS))
  reml <- as.numeric(estimate_sac_range(fld, "z", TREND_VARS, detrend = "reml"))
  expect_true(is.finite(ols) && is.finite(reml))
  expect_gt(reml, ols)
})

test_that("reml_max_n subsamples the REML fit reproducibly and applies the trend to every point", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  fld <- sac_trend_field()
  r1 <- estimate_sac_range(fld, "z", TREND_VARS, detrend = "reml", reml_max_n = 100)
  r2 <- estimate_sac_range(fld, "z", TREND_VARS, detrend = "reml", reml_max_n = 100)
  expect_identical(attr(r1, "reml")$n_used, 100L)
  expect_true(attr(r1, "reml")$subsampled)
  expect_identical(as.numeric(r1), as.numeric(r2))       # seeded subsample
  # The residual variogram still covers every point, not just the subsample.
  expect_identical(sum(attr(r1, "variogram")$np),
                   sum(attr(estimate_sac_range(fld, "z", TREND_VARS), "variogram")$np))
  # Exact duplicate locations are dropped for the fit and do not break it.
  dup <- rbind(fld[1:60, ], fld[1:20, ])
  rd <- estimate_sac_range(dup, "z", TREND_VARS, detrend = "reml")
  expect_identical(attr(rd, "detrend_method"), "reml")
  expect_identical(attr(rd, "reml")$n_used, 60L)
})

test_that("a REML fit that does not converge falls back to OLS with a warning", {
  skip_if_not_installed("gstat")
  fld <- sac_trend_field()
  local_mocked_bindings(.reml_trend = function(...) NULL, .package = "spatialkit")
  expect_warning(
    r <- estimate_sac_range(fld, "z", TREND_VARS, detrend = "reml"),
    "did not converge.*falling back to OLS")
  expect_identical(attr(r, "detrend_method"), "ols")
  expect_null(attr(r, "reml"))
  expect_equal(as.numeric(r), as.numeric(estimate_sac_range(fld, "z", TREND_VARS)))
})


# ---------------------------------------------------------------------------
# The nugget, on every classed return path.
# ---------------------------------------------------------------------------

test_that("sac_nugget reads the fitted model's nugget and is NA where there is none", {
  skip_if_not_installed("gstat")
  r <- estimate_sac_range(sac_test_field(), "z")
  expect_true(is.finite(r))
  vm <- attr(r, "variogram_model")
  expect_equal(sac_nugget(r), sum(vm$psill[vm$model == "Nug"]))
  expect_identical(sac_nugget(r), attr(r, "nugget"))
  expect_true(sac_nugget(r) >= 0)
  # A rejected range still carries the nugget of the model that was refused.
  rej <- suppressWarnings(estimate_sac_range(sac_test_field(), "z", range_frac = 1e-6))
  expect_true(is.na(rej))
  expect_true(is.finite(sac_nugget(rej)))
  # Nothing fitted: NA, never an error.
  expect_identical(sac_nugget(NA), NA_real_)
  expect_identical(sac_nugget(NA_real_), NA_real_)
  expect_identical(sac_nugget(NULL), NA_real_)
  expect_identical(sac_nugget(42), NA_real_)
  expect_identical(sac_nugget(structure(NA_real_, class = c("sac_range", "numeric"))),
                   NA_real_)
  expect_identical(spatialkit:::.vgm_nugget_of(NULL), NA_real_)
  expect_identical(spatialkit:::.vgm_nugget_of(data.frame(model = "Exp", psill = 1)), 0)
})


# ---------------------------------------------------------------------------
# A variogram that falls with distance identifies no range.
# ---------------------------------------------------------------------------

test_that(".variogram_decreasing flags a net fall over the shorter lags and nothing else", {
  f <- spatialkit:::.variogram_decreasing
  mk <- function(gamma, np = rep(100, length(gamma)))
    data.frame(np = np, dist = seq(10, by = 10, length.out = length(gamma)), gamma = gamma)
  expect_false(f(mk(c(1, 2, 3, 4, 5, 6, 7, 8))))            # rising
  expect_false(f(mk(rep(5, 8))))                            # flat
  expect_true(f(mk(c(8, 7, 6, 5, 4, 3, 2, 1))))             # falling
  # Only the shorter half counts: a fall confined to the long lags is not it.
  expect_false(f(mk(c(1, 2, 3, 4, 8, 6, 4, 2))))
  # Tolerance: a fall of 10% of the mean is noise, 20% is not (tol = 0.15).
  expect_false(f(mk(c(10, 10, 10, 9), np = c(100, 100, 100, 100)), frac = 1))
  expect_true(f(mk(c(10, 10, 10, 8), np = c(100, 100, 100, 100)), frac = 1))
  # Weighted by the pairs supporting each step: a fall carried by one bin
  # with hardly any pairs does not count.
  expect_false(f(mk(c(10, 10, 10, 2), np = c(100, 100, 100, 1)), frac = 1))
  # Too few bins, and malformed input, are never "decreasing".
  expect_false(f(mk(c(3, 2, 1))))
  expect_false(f(NULL))
  expect_false(f(data.frame(dist = 1:5, gamma = 5:1)))       # no np column
  expect_false(f(mk(c(NA, NA, NA, NA, NA))))
})

test_that("a periodic field is refused as 'decreases with distance', with the evidence attached", {
  skip_if_not_installed("gstat")
  # A hole-effect variogram: the semivariance falls again past the first
  # quarter-wavelength.  gstat still fits an exponential to it and reports a
  # finite range; that number is not a correlation length.
  set.seed(5001)
  n <- 250
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  S  <- as.numeric(t(chol(exp(-D / 100) + diag(1e-8, n))) %*% rnorm(n))
  per <- sf::st_as_sf(
    data.frame(x = xy[, 1], y = xy[, 2],
               z = 2 * sin(2 * pi * xy[, 1] / 250) + 0.3 * S + rnorm(n, sd = 0.3)),
    coords = c("x", "y"), crs = 32632)
  lines <- capture_spatialkit_log(r <- estimate_sac_range(per, "z"))
  expect_true(is.na(r))
  expect_s3_class(r, "sac_range")
  expect_identical(attr(r, "rejected_reason"), "empirical variogram decreases with distance")
  expect_true(is.finite(attr(r, "rejected_range")))
  expect_s3_class(attr(r, "variogram"), "data.frame")
  expect_true(spatialkit:::.variogram_decreasing(attr(r, "variogram")))
  expect_true(log_has(lines, "decreases with distance"))
  expect_true(log_has(lines, "periodic"))
  expect_false(log_has(lines, "trend --"))
  # It draws, and says why no range is marked.
  skip_if_not_installed("ggplot2")
  p <- plot(r)
  expect_match(p$labels$subtitle, "^No effective range: the semivariance falls")
  # And an ordinary field is not touched by the check.
  ok <- estimate_sac_range(sac_test_field(), "z")
  expect_true(is.finite(ok))
  expect_false(spatialkit:::.variogram_decreasing(attr(ok, "variogram")))
})
