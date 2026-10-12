# tests/testthat/test-fix-variogram-lags.R
# ---------------------------------------------------------------------------
# The lags start at the spacing of the closest pairs, and the nugget says how
# far it can be trusted.
#
# The variogram was binned in 15 equal lags out to half the extent, whatever
# the spacing of the points.  On a clustered layer whose structure is shorter
# than the first bin (900 points in 60 clusters, effective range 30 against a
# bin 38 wide, nugget 0.4 on a sill of 1) the range came back at a median 77
# with 10 of 28 within a factor 2 of the truth, and the nugget at 0.72, exactly
# 0 in 3 fits: all structure shorter than the first bin was counted as
# nugget, or none of the nugget was.
# ---------------------------------------------------------------------------

# The clustered layout of the review's reproduction, and one draw of a field
# with effective range 30, partial sill 0.6 and nugget 0.4 on it.
fvl_sites <- local({
  val <- NULL
  function() {
    if (is.null(val)) {
      set.seed(11)
      cen <- cbind(runif(60, 0, 1000), runif(60, 0, 1000))
      xy <- cen[rep(1:60, each = 15), ] + matrix(rnorm(2 * 900, sd = 8), ncol = 2)
      val <<- list(xy = xy, L = t(chol(0.6 * exp(-3 * as.matrix(stats::dist(xy)) / 30) +
                                         diag(1e-8, 900))))
    }
    val
  }
})
fvl_clustered <- function(i) {
  s <- fvl_sites()
  set.seed(100 + i)
  z <- as.numeric(s$L %*% rnorm(900)) + sqrt(0.4) * rnorm(900)
  sf::st_as_sf(data.frame(x = 5e5 + s$xy[, 1], y = 5e6 + s$xy[, 2], z = z),
               coords = c("x", "y"), crs = 32632)
}
# An exponential field on uniform sites: effective range R, partial sill ps,
# nugget nug.
fvl_uniform <- function(seed, n = 300L, R = 300, ps = 0.8, nug = 0.2) {
  set.seed(seed)
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  z <- drop(crossprod(chol(ps * exp(-3 * as.matrix(stats::dist(xy)) / R) + diag(1e-8, n)),
                      rnorm(n))) + rnorm(n, 0, sqrt(nug))
  sf::st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2], z = z),
               coords = c("x", "y"), crs = 32632)
}

test_that("a range shorter than an equal first bin is recovered on a clustered layer", {
  skip_if_not_installed("gstat")
  draws <- c(2, 4, 5, 6, 7, 8)
  res <- lapply(draws, function(i) estimate_sac_range(fvl_clustered(i), "z"))
  rng <- vapply(res, as.numeric, numeric(1))
  nug <- vapply(res, sac_nugget, numeric(1))
  # Every draw gets a range, and nearly every one within a factor 2 of 30
  # (the equal lags gave 409, 87, 21, 361, 79 and 28 on these six).
  expect_true(all(is.finite(rng)))
  expect_gte(sum(rng > 15 & rng < 60), 5L)
  # The nugget is 0.4: within 0.2 of it in nearly every draw, and never at
  # its bound (the equal lags gave 0.82, 0.68, 0, 0.76, 0.79 and 0).
  expect_gte(sum(abs(nug - 0.4) <= 0.2), 5L)
  expect_false(any(vapply(res, function(r) isTRUE(attr(r, "nugget_at_bound")), logical(1))))
  # The first lag is at the spacing of the points, far inside the first
  # equal bin (38 wide, mean separation 17).
  expect_true(all(vapply(res, function(r) attr(r, "first_lag"), numeric(1)) < 5))
})

test_that("lags = 'equal' keeps the equal bins, and says where its nugget starts from", {
  skip_if_not_installed("gstat")
  pts <- fvl_clustered(4)
  a <- estimate_sac_range(pts, "z")
  e <- suppressWarnings(estimate_sac_range(pts, "z", lags = "equal"))
  expect_identical(attr(a, "lags"), "adaptive")
  expect_identical(attr(e, "lags"), "equal")
  # Fifteen equal bins out to the cutoff.
  be <- attr(e, "lag_boundaries")
  expect_length(be, 15L)
  expect_equal(be, attr(e, "cutoff_dist") * (1:15) / 15)
  expect_lte(nrow(attr(e, "variogram")), 15L)
  # Geometric bins, the last at the cutoff, the first far below an equal bin.
  ba <- attr(a, "lag_boundaries")
  expect_equal(ba[length(ba)], attr(a, "cutoff_dist"))
  expect_equal(diff(log(ba)), rep(diff(log(ba))[1], length(ba) - 1L))
  expect_lt(ba[1], be[1] / 5)
  expect_lte(length(ba), 16L)
  # first_lag is the mean separation in the first bin of the variogram used.
  for (r in list(a, e)) {
    vg <- attr(r, "variogram")
    expect_equal(attr(r, "first_lag"), min(vg$dist[vg$np > 0]))
  }
  expect_gt(attr(e, "first_lag"), 5 * attr(a, "first_lag"))
  # n_lags sets the number of equal bins.
  e10 <- suppressWarnings(estimate_sac_range(pts, "z", lags = "equal", n_lags = 10))
  expect_length(attr(e10, "lag_boundaries"), 10L)
  expect_lte(nrow(attr(e10, "variogram")), 10L)
  expect_equal(attr(e10, "lag_boundaries")[10], attr(e10, "cutoff_dist"))
})

test_that("with equal lags the fit minimises the criterion gstat's fit minimises", {
  skip_if_not_installed("gstat")
  # A well-determined field, where gstat's iteration converges: the two fits
  # of the same variogram agree.
  r <- estimate_sac_range(fvl_uniform(4201), "z", lags = "equal")
  expect_true(is.finite(r))
  vg <- attr(r, "variogram")
  g <- suppressWarnings(gstat::fit.variogram(vg, gstat::vgm(psill = NA, model = "Exp",
                                                            range = NA, nugget = NA)))
  skip_if(isTRUE(attr(g, "singular")), "gstat's own fit is singular on this platform")
  vm <- attr(r, "variogram_model")
  expect_equal(vm$range[2], g$range[2], tolerance = 0.02)
  expect_equal(vm$psill, g$psill, tolerance = 0.02)
  # And at least as good by that criterion: sum of np / h^2 (gamma - model)^2.
  sse <- function(m) {
    fit <- m$psill[1] + m$psill[2] * (1 - exp(-vg$dist / m$range[2]))
    sum(vg$np / vg$dist^2 * (vg$gamma - fit)^2)
  }
  expect_lte(sse(vm), sse(g) * (1 + 1e-6))
})

test_that("every classed result carries the short-lag evidence", {
  skip_if_not_installed("gstat")
  pts <- fvl_uniform(4202)
  new_attrs <- c("first_lag", "nugget_at_bound", "bins_below_range",
                 "short_lag_semivariance", "short_lag_pairs", "short_lag_dist",
                 "structure_p", "morans_i", "lags", "lag_boundaries", "structures",
                 "short_range", "n_used", "n_missing", "n_subsampled")
  ok  <- estimate_sac_range(pts, "z")
  rej <- suppressWarnings(estimate_sac_range(pts, "z", range_frac = 1e-6))
  set.seed(1); pts$noise <- rnorm(nrow(pts))
  flat <- suppressWarnings(estimate_sac_range(pts, "noise"))
  expect_true(is.finite(ok)); expect_true(is.na(rej)); expect_true(is.na(flat))
  for (r in list(ok, rej, flat)) {
    expect_true(all(new_attrs %in% names(attributes(r))),
                info = paste(setdiff(new_attrs, names(attributes(r))), collapse = ", "))
    # ... and every attribute it carried before.
    expect_true(all(c("directional", "anisotropy", "anisotropy_used", "directional_status",
                      "directional_fitted", "detrended", "detrend_method", "max_dist",
                      "cutoff_dist", "crs", "variogram", "variogram_model", "nugget") %in%
                      names(attributes(r))))
  }
  # The semivariance of nearest-neighbour pairs, against its definition: each
  # point paired with its nearest other point, each pair counted once.
  xy <- sf::st_coordinates(pts); z <- pts$z
  D <- as.matrix(stats::dist(xy)); diag(D) <- Inf
  nn <- apply(D, 1, which.min)
  pr <- unique(t(apply(cbind(seq_along(nn), nn), 1, sort)))
  expect_equal(attr(ok, "short_lag_semivariance"), mean(0.5 * (z[pr[, 1]] - z[pr[, 2]])^2))
  expect_identical(attr(ok, "short_lag_pairs"), nrow(pr))
  expect_equal(attr(ok, "short_lag_dist"), stats::median(D[pr]))
  # It bounds the nugget from above: the field's nugget is 0.2.
  expect_gt(attr(ok, "short_lag_semivariance"), 0.2)
  expect_lt(sac_nugget(ok), attr(ok, "short_lag_semivariance") + 0.1)
  # Bins inside the fitted range, counted on the variogram attached.
  expect_identical(attr(ok, "bins_below_range"),
                   sum(attr(ok, "variogram")$dist < as.numeric(ok)))
  expect_identical(attr(ok, "structures"), 1L)
  expect_true(is.na(attr(ok, "short_range")))
})

test_that("a nugget at its lower bound is marked, not handed on bare", {
  skip_if_not_installed("gstat")
  # A field with no nugget at all: the fit puts the nugget at 0, which is the
  # bound of the search and not a measurement of it.
  lines <- capture_spatialkit_log(
    r <- estimate_sac_range(fvl_uniform(603, ps = 1, nug = 0), "z"))
  expect_true(is.finite(r))
  expect_identical(sac_nugget(r), 0)
  expect_identical(attr(r, "nugget_at_bound"), TRUE)
  expect_true(log_has(lines, "is at its lower bound"))
  expect_true(log_has(lines, "nugget_at_bound"))
  # The closest pairs still say something: their semivariance is positive.
  expect_gt(attr(r, "short_lag_semivariance"), 0)
  # A fitted nugget away from the bound is not marked.
  r2 <- estimate_sac_range(fvl_uniform(4202, nug = 0.5), "z")
  expect_gt(sac_nugget(r2), 0.1)
  expect_identical(attr(r2, "nugget_at_bound"), FALSE)
  # A REML nugget stops short of exactly 0; within 1e-4 of the sill counts.
  hand <- gstat::vgm(psill = 0.99, model = "Exp", range = 100, nugget = 6e-7)
  expect_lte(spatialkit:::.vgm_nugget_of(hand), 1e-4 * sum(hand$psill))
})

test_that("every fitted nugget of zero is marked as the bound it is", {
  skip_if_not_installed("gstat")
  # 600 points over 400 x 300 km, effective range 30 km, noise 0.2 on a sill
  # of 1.2: the range spans a few bins and the nugget is weakly identified.
  # Of 20 such fields the fitted nugget is exactly 0 in 3 (it was in 4 with
  # equal lags, with nothing to tell them from an estimate).
  one <- function(i) {
    set.seed(5200 + i)
    n <- 600L
    x <- runif(n, 0, 400e3); y <- runif(n, 0, 300e3)
    D <- as.matrix(stats::dist(cbind(x, y)))
    z <- drop(crossprod(chol(exp(-3 * D / 30e3) + diag(1e-8, n)), rnorm(n))) +
      rnorm(n, 0, sqrt(0.2))
    estimate_sac_range(sf::st_as_sf(data.frame(x = 5e5 + x, y = 4e6 + y, z = z),
                                    coords = c("x", "y"), crs = 32615), "z")
  }
  res <- lapply(c(1, 2, 3, 5, 7, 14), one)
  at_bound <- vapply(res, function(r) attr(r, "nugget_at_bound"), logical(1))
  nug  <- vapply(res, sac_nugget, numeric(1))
  sill <- vapply(res, function(r) sum(attr(r, "variogram_model")$psill), numeric(1))
  expect_true(all(vapply(res, is.finite, logical(1))))
  # The mark is exactly "at or within 1e-4 of the sill of zero".
  expect_identical(at_bound, nug <= 1e-4 * sill)
  expect_true(any(at_bound)); expect_true(any(!at_bound))
  expect_true(all(nug[at_bound] == 0))
  # Where the nugget is at its bound the closest pairs still give a number:
  # their semivariance, well above 0, is what a nugget cannot exceed.
  sls <- vapply(res, function(r) attr(r, "short_lag_semivariance"), numeric(1))
  expect_true(all(is.finite(sls) & sls > 0.2))
  expect_true(all(nug <= sls))
})

test_that("subsamples of one layer agree on its range", {
  skip_if_not_installed("gstat")
  # 20 subsamples of 80 percent of one clustered layer (range 30).  With
  # equal lags the 20 ranges ran from 63 to 153 and beyond (largest over
  # smallest 5.4, median 94).
  pts <- fvl_clustered(4)
  r <- vapply(1:20, function(b) {
    set.seed(3000 + b)
    as.numeric(suppressWarnings(estimate_sac_range(pts[sort(sample(900, 720)), ], "z")))
  }, numeric(1))
  expect_gte(sum(is.finite(r)), 18L)
  expect_lt(max(r, na.rm = TRUE) / min(r, na.rm = TRUE), 3)
  expect_gt(stats::median(r, na.rm = TRUE), 20)
  expect_lt(stats::median(r, na.rm = TRUE), 45)
})

test_that("bins with too few pairs are merged, pair-weighted, and none is dropped", {
  f <- spatialkit:::.sac_merge_bins
  vg <- data.frame(np = c(10, 10, 15, 40, 5), dist = c(1, 2, 3, 4, 5),
                   gamma = c(1, 2, 3, 4, 5))
  m <- f(vg, min_np = 30)
  # 10 + 10 + 15 make the first bin; 40 the second, and the tail of 5 joins it.
  expect_equal(m$np, c(35, 45))
  expect_equal(sum(m$np), sum(vg$np))
  expect_equal(m$dist, c((10 * 1 + 10 * 2 + 15 * 3) / 35, (40 * 4 + 5 * 5) / 45))
  expect_equal(m$gamma, c((10 * 1 + 10 * 2 + 15 * 3) / 35, (40 * 4 + 5 * 5) / 45))
  # Nothing to do when every bin is large enough, or there is one bin.
  big <- data.frame(np = c(30, 50), dist = 1:2, gamma = 1:2)
  expect_identical(f(big, min_np = 30), big)
  expect_identical(f(vg[1, ], min_np = 30), vg[1, ])
  # Fewer pairs in all than one bin needs: a single bin with all of them.
  expect_equal(f(vg[1:3, ], min_np = 100)$np, 35)
})

test_that("the adaptive boundaries start where the closest pairs end", {
  f <- spatialkit:::.sac_lag_boundaries
  set.seed(41); n <- 400
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  b <- f(xy, cutoff_dist = 700)
  expect_true(all(diff(b) > 0))
  expect_equal(b[length(b)], 700)
  expect_lte(length(b), 16L)                 # the first bin and up to 15 more
  expect_gte(length(b), 6L)
  expect_equal(diff(log(b)), rep(diff(log(b))[1], length(b) - 1L))
  # The first bin holds about max(50, n / 4) = 100 pairs, and no fewer.
  d <- as.numeric(stats::dist(xy))
  expect_gte(sum(d <= b[1]), 100L)
  expect_lte(sum(d <= b[1]), 200L)
  # It is capped at a quarter of the cutoff on a layer too small to fill it.
  small <- f(xy[1:12, ], cutoff_dist = 700)
  expect_lte(small[1], 700 / 4 * (1 + 1e-12))
  # Repeat visits to a location do not break it: zero distances are not a lag.
  dup <- rbind(xy, xy[1:100, ])
  bd <- f(dup, cutoff_dist = 700)
  expect_true(all(is.finite(bd)) && all(bd > 0) && all(diff(bd) > 0))
  # The same boundaries in any unit of length.
  expect_equal(f(xy / 1000, cutoff_dist = 0.7) * 1000, b, tolerance = 1e-12)
})

test_that("a variogram that rises steeply on geometric bins is not called decreasing", {
  f <- spatialkit:::.variogram_decreasing
  # Twelve bins a factor 2 apart whose pair counts grow by 3: a steep rise
  # over the first five, then a sag of a tenth over the long lags.  Weighted
  # by pairs, the last steps spoke for all of them and this was refused.
  vg <- data.frame(np = 60 * 3^(0:11), dist = 2^(0:11),
                   gamma = c(0.40, 0.55, 0.70, 0.85, 0.97, 1.00, 1.00, 0.99, 0.97,
                             0.95, 0.93, 0.90))
  expect_false(f(vg))
  # A variogram that falls from its first bin on is still caught ...
  fall <- vg; fall$gamma <- seq(1.6, 0.8, length.out = 12)
  expect_true(f(fall))
  # ... and so is a hump: up to a peak in the short lags and back down.
  hump <- vg; hump$gamma <- c(0.5, 0.6, 0.8, 1.0, 1.3, 1.6, 1.9, 2.2, 2.4, 1.6, 1.0, 0.9)
  expect_true(f(hump))
  # On a layer: an ordinary field on the clustered layout is not refused for it.
  skip_if_not_installed("gstat")
  r <- estimate_sac_range(fvl_clustered(2), "z")
  expect_false(f(attr(r, "variogram")))
  expect_false(identical(attr(r, "rejected_reason"), "empirical variogram decreases with distance"))
})

test_that("on uniform points the adaptive lags lose nothing to the equal ones", {
  skip_if_not_installed("gstat")
  # 300 points, effective range 300: the reference case of the old help page.
  seeds <- 4301:4312
  a <- vapply(seeds, function(s) as.numeric(suppressWarnings(
    estimate_sac_range(fvl_uniform(s), "z"))), numeric(1))
  e <- vapply(seeds, function(s) as.numeric(suppressWarnings(
    estimate_sac_range(fvl_uniform(s), "z", lags = "equal"))), numeric(1))
  expect_gte(sum(is.finite(a)), 11L)
  expect_gte(sum(is.finite(a)), sum(is.finite(e)) - 1L)
  ma <- stats::median(a, na.rm = TRUE); me <- stats::median(e, na.rm = TRUE)
  expect_gt(ma, 200); expect_lt(ma, 450)
  expect_lt(abs(ma / me - 1), 0.25)
})

# A field with structure at two scales on the clustered layout: effective
# ranges `rs` (partial sill 0.4) and `rl` (0.3) over a nugget of 0.3.
fvl_nested <- function(seed, rs, rl) {
  s <- fvl_sites()
  D <- as.matrix(stats::dist(s$xy))
  set.seed(seed)
  z <- drop(crossprod(chol(0.4 * exp(-3 * D / rs) + 0.3 * exp(-3 * D / rl) +
                             diag(1e-8, 900)), rnorm(900))) + rnorm(900, 0, sqrt(0.3))
  sf::st_as_sf(data.frame(x = 5e5 + s$xy[, 1], y = 5e6 + s$xy[, 2], z = z),
               coords = c("x", "y"), crs = 32632)
}

test_that("a much longer structure is kept in the model and sets the range", {
  skip_if_not_installed("gstat")
  # Structure at two scales on the clustered layout: 15 (partial sill 0.4)
  # and 450 (0.3) over a nugget of 0.3.  One structure fitted to lags that
  # resolve the points follows the short one; blocks have to clear the long.
  pts <- fvl_nested(9211, 15, 450)
  lines <- capture_spatialkit_log(r <- estimate_sac_range(pts, "z"))
  e <- estimate_sac_range(pts, "z", lags = "equal")
  expect_identical(attr(r, "structures"), 2L)
  vm <- attr(r, "variogram_model")
  expect_identical(as.character(vm$model), c("Nug", "Exp", "Exp"))
  expect_true(all(vm$psill[2:3] > 0))
  # The long range is the equal-lag fit's, it is inside the lags, and it is
  # the estimate (277 here; the single adaptive fit gave 24).
  expect_true(is.finite(r)); expect_true(is.finite(e))
  expect_equal(as.numeric(r), as.numeric(e), tolerance = 1e-8)
  expect_equal(as.numeric(r), 3 * vm$range[3])
  expect_lt(as.numeric(r), attr(r, "cutoff_dist"))
  expect_false(attr(r, "lag_fallback"))
  expect_identical(attr(r, "lags"), "adaptive")
  # The short one is at most half as long, and is reported.
  expect_equal(attr(r, "short_range"), 3 * vm$range[2])
  expect_lte(2 * attr(r, "short_range"), as.numeric(r))
  expect_lt(attr(r, "short_range"), 60)
  # The long structure carries its share of the semivariance, and both
  # fits have a minimum, so the model is marked as one.
  expect_gte(attr(vm, "long_share"), 0.03)
  expect_identical(attr(vm, "converged"), TRUE)
  # The nugget comes from the short lags: nearer 0.3 than the equal-lag
  # fit's, which counts the short structure as nugget (0.26 against 0.59).
  expect_lt(abs(sac_nugget(r) - 0.3), abs(sac_nugget(e) - 0.3))
  expect_true(log_has(lines, "structure at two scales"))
  # Bins inside the short structure are what the nugget is told apart from.
  expect_identical(attr(r, "bins_below_range"),
                   sum(attr(r, "variogram")$dist < attr(r, "short_range")))
})

test_that("a long range at the edge of the search marks the two-structure model as no minimum", {
  skip_if_not_installed("gstat")
  # On this draw the equal-lag variogram has no sill inside four cutoffs: its
  # criterion is still falling at the edge of the ranges searched, so the
  # long range (2292 = 4 x 573) and the partial sill solved at it are that
  # edge.  The model built on it was marked converged = TRUE, which the help
  # says of a fitted sill only.
  pts <- fvl_nested(9203, 15, 450)
  lines <- capture_spatialkit_log(r <- suppressWarnings(estimate_sac_range(pts, "z")))
  e <- suppressWarnings(estimate_sac_range(pts, "z", lags = "equal"))
  expect_identical(attr(attr(e, "variogram_model"), "converged"), FALSE)
  expect_true(is.na(r))
  expect_identical(attr(r, "rejected_reason"), "fitted range exceeds the largest lag fitted")
  expect_identical(attr(r, "structures"), 2L)
  expect_equal(attr(r, "rejected_range"), 4 * attr(r, "cutoff_dist"))
  vm <- attr(r, "variogram_model")
  expect_identical(attr(vm, "converged"), FALSE)
  expect_identical(spatialkit:::.vgm_converged(vm), FALSE)
  expect_true(log_has(lines, "the model is marked converged = FALSE"))
  # A model so marked gives Cp no nugget: the flag is one bit for the whole
  # model, and resolution_profile() takes nothing from a model without a
  # minimum.
  out <- suppressWarnings(resolution_profile(pts, "z", sac = r, n_levels = 4, nstart = 2))
  expect_identical(attr(out, "cp_noise")$source, "finest-level residual mean square")
  # The helper passes the flag through, and defaults to a minimum.
  f <- spatialkit:::.sac_fit_nested
  vg <- attr(r, "variogram")
  expect_identical(attr(f(vg, attr(r, "cutoff_dist"), 30, 400), "converged"), TRUE)
  expect_identical(attr(f(vg, attr(r, "cutoff_dist"), 30, 400, converged = FALSE),
                        "converged"), FALSE)
})

test_that("a weaker long structure is kept too: the range does not drop under the equal-lag one", {
  skip_if_not_installed("gstat")
  # Ranges 20 and 200.  A long structure was kept only when the equal-lag
  # range was three times the adaptive one and carried a tenth of the
  # semivariance, and on 30 such layers the range returned was under half
  # the equal-lag one in 8.  Blocks sized from it leak.  The bars are a
  # factor two and 3 percent now.
  for (seed in c(9301, 9309, 9311)) {
    pts <- fvl_nested(seed, 20, 200)
    r <- suppressWarnings(estimate_sac_range(pts, "z"))
    e <- suppressWarnings(estimate_sac_range(pts, "z", lags = "equal"))
    expect_true(is.finite(r) && is.finite(e), info = seed)
    expect_identical(attr(r, "structures"), 2L, info = seed)
    expect_equal(as.numeric(r), as.numeric(e), tolerance = 1e-8, info = seed)
    expect_gte(as.numeric(r), 100)
    # The short structure is still there to read, near its 20.
    expect_lt(attr(r, "short_range"), 50)
  }
  # On 9309 the long structure has 7 percent, under the old bar: the range
  # was 34 where the equal lags give 142.
  r <- suppressWarnings(estimate_sac_range(fvl_nested(9309, 20, 200), "z"))
  sh <- attr(attr(r, "variogram_model"), "long_share")
  expect_gt(sh, 0.03); expect_lt(sh, 0.10)
  # A layer with one short structure keeps its short range where the
  # adaptive variogram gives the equal-lag range nothing (draw 6: 361 with
  # equal lags, and under 3 percent).
  one <- estimate_sac_range(fvl_clustered(6), "z")
  eq1 <- estimate_sac_range(fvl_clustered(6), "z", lags = "equal")
  expect_identical(attr(one, "structures"), 1L)
  expect_gt(as.numeric(eq1), 4 * as.numeric(one))
  expect_lt(as.numeric(one), 60)
})

test_that("a range the adaptive lags leave past the lags falls back to the equal-lag result", {
  skip_if_not_installed("gstat")
  # 400 uniform points, effective range 250, partial sill 0.3 under a nugget
  # of 0.7.  The adaptive lags put little weight on the long ones and their
  # fit runs past the largest lag; the equal lags identify 381.  The default
  # refused 14 of 120 such layers that lags = "equal" answers.
  set.seed(323); n <- 400
  u <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  z <- drop(crossprod(chol(0.3 * exp(-3 * as.matrix(stats::dist(u)) / 250) + diag(1e-8, n)),
                      rnorm(n))) + rnorm(n, 0, sqrt(0.7))
  pts <- sf::st_as_sf(data.frame(x = 5e5 + u[, 1], y = 5e6 + u[, 2], z = z),
                      coords = c("x", "y"), crs = 32632)
  lines <- capture_spatialkit_log(r <- estimate_sac_range(pts, "z"))
  e <- estimate_sac_range(pts, "z", lags = "equal")
  expect_true(is.finite(e))
  expect_true(is.finite(r))
  expect_identical(attr(r, "lag_fallback"), TRUE)
  expect_identical(attr(e, "lag_fallback"), FALSE)
  # The equal-lag result, whole: the value, the variogram, the model, the
  # directional sweep, the lags and the first lag.
  same <- setdiff(names(attributes(e)), "lag_fallback")
  expect_identical(attributes(r)[same], attributes(e)[same])
  expect_identical(as.numeric(r), as.numeric(e))
  expect_identical(attr(r, "lags"), "equal")
  expect_length(attr(r, "lag_boundaries"), 15L)
  expect_true(log_has(lines, "The equal-lag result is returned"))
  expect_true(log_has(lines, "lag_fallback"))
  # The adaptive fit it replaces is past the lags, with nothing to fall
  # back to once the equal-lag range is refused as well.
  tight <- suppressWarnings(estimate_sac_range(pts, "z", range_frac = 0.3))
  expect_true(is.na(tight))
  expect_identical(attr(tight, "rejected_reason"), "fitted range exceeds the largest lag fitted")
  expect_identical(attr(tight, "lag_fallback"), FALSE)
  expect_identical(attr(tight, "lags"), "adaptive")
  # No other refusal falls back.  Two white-noise layers, let through with
  # the test for structure off: the adaptive variogram of one is flat and
  # that of the other falls with distance, the equal lags fit a range to
  # each (181 and 223), and neither is handed on.
  for (case in list(list(70302, "the empirical variogram is flat: no spatial structure at the lags resolved"),
                    list(70329, "empirical variogram decreases with distance"))) {
    set.seed(case[[1]])
    wn <- sf::st_as_sf(data.frame(x = 5e5 + runif(300, 0, 1000), y = 5e6 + runif(300, 0, 1000),
                                  z = rnorm(300)), coords = c("x", "y"), crs = 32632)
    a <- suppressWarnings(estimate_sac_range(wn, "z", structure_alpha = 1))
    expect_true(is.finite(suppressWarnings(
      estimate_sac_range(wn, "z", structure_alpha = 1, lags = "equal"))))
    expect_true(is.na(a))
    expect_identical(attr(a, "rejected_reason"), case[[2]])
    expect_identical(attr(a, "lag_fallback"), FALSE)
    expect_identical(attr(a, "lags"), "adaptive")
  }
})
