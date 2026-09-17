# tests/testthat/test-resolution-profile.R
# ---------------------------------------------------------------------------
# The resolution project: a smooth WSS curve (k-means++ restarts), a bumps
# detector on it, and a per-level profile scored on several criteria with the
# bounds the data impose and the region over which a criterion is flat.
# ---------------------------------------------------------------------------

# An exponential Gaussian field with a known nugget, sampled at random
# locations, plus a white-noise covariate the response partly depends on.
rp_field <- function(n = 400, a = 100, nugget = 0.3, seed = 2, with_pred = TRUE) {
  set.seed(seed)
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  S  <- as.numeric(t(chol(exp(-D / a) + diag(1e-8, n))) %*% rnorm(n))
  xy$w <- rnorm(n)
  xy$z <- (if (with_pred) 2 * xy$w else 0) + S + rnorm(n, sd = sqrt(nugget))
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
}

# Eight tight, well-separated clusters: the layout on which a handful of
# random restarts leaves bumps on the WSS curve.
rp_clustered <- function(n = 1500, seed = 4) {
  set.seed(seed)
  ctr <- cbind(runif(8, 0, 1000), runif(8, 0, 1000))
  g   <- sample(8, n, TRUE)
  ctr[g, ] + matrix(rnorm(2 * n, sd = 60), n)
}


test_that(".kmeanspp_centers draws k distinct points from the data", {
  set.seed(1)
  xy <- cbind(runif(200), runif(200))
  ctr <- spatialkit:::.kmeanspp_centers(xy, 7L)
  expect_identical(dim(ctr), c(7L, 2L))
  expect_false(anyDuplicated(ctr) > 0)
  # Every centre is a data point.
  expect_true(all(apply(ctr, 1L, function(r) any(xy[, 1] == r[1] & xy[, 2] == r[2]))))
  expect_identical(dim(spatialkit:::.kmeanspp_centers(xy, 1L)), c(1L, 2L))
  # Reproducible under the seed.
  set.seed(9); a <- spatialkit:::.kmeanspp_centers(xy, 5L)
  set.seed(9); b <- spatialkit:::.kmeanspp_centers(xy, 5L)
  expect_identical(a, b)
})

test_that(".kmeans_best keeps the best restart and reports the spread", {
  xy <- rp_clustered(400)
  set.seed(3)
  kb <- spatialkit:::.kmeans_best(xy, 8L, nstart = 10L)
  expect_s3_class(kb$km, "kmeans")
  expect_identical(kb$wss, kb$km$tot.withinss)
  expect_gte(kb$spread, 0)
  expect_identical(kb$n_ok, 10L)
  # The best of ten is no worse than any single k-means++ start.
  set.seed(3)
  singles <- replicate(10, stats::kmeans(xy, centers = spatialkit:::.kmeanspp_centers(xy, 8L),
                                         iter.max = 50, nstart = 1)$tot.withinss)
  expect_lte(kb$wss, min(singles) + 1e-6)
})

test_that(".wss_bumps counts the rises on a curve", {
  f <- spatialkit:::.wss_bumps
  expect_identical(f(c(10, 8, 6, 5, 4.5)), 0L)
  expect_identical(f(c(10, 8, 9, 5, 6)), 2L)
  expect_identical(f(c(10, NA, 6)), 0L)
  expect_identical(f(5), 0L)
})

test_that("25 k-means++ restarts leave the WSS curve monotone where 5 random ones did not", {
  # Measured: on this layout stats::kmeans(nstart = 5) rose at two steps
  # of a 1..30 sweep; the restart budget is the fix, and this pins it.
  xy <- rp_clustered()
  set.seed(4)
  wss <- vapply(2:30, function(k) spatialkit:::.kmeans_best(xy, k, nstart = 25L)$wss, numeric(1))
  expect_identical(spatialkit:::.wss_bumps(wss), 0L)
})

test_that("determine_optimal_levels warns before the sweep when no k can clear the floor", {
  pts <- rp_field(120, with_pred = TRUE)
  lines <- capture_spatialkit_log(
    out <- determine_optimal_levels(pts, max_levels = 6, response_var = "z",
                                    predictor_vars = "w", criterion = "morans_i"))
  expect_true(log_has(lines, "nine cells or fewer"))
  expect_true(log_has(lines, "Raise max_levels"))
  expect_type(out, "integer")
  # With room above the floor the warning is not raised.
  quiet <- capture_spatialkit_log(
    determine_optimal_levels(pts, max_levels = 14, response_var = "z",
                             predictor_vars = "w", criterion = "morans_i"))
  expect_false(log_has(quiet, "nine cells or fewer"))
})

test_that("determine_optimal_levels reports the restart budget and bumps in its diagnostics", {
  # Sixteen separated clusters and a ladder to 30: the elbow's evaluation
  # window reaches past the nine-cell floor on this draw, so the model-aware
  # path returns its diagnostics rather than falling back.
  set.seed(5)
  ctr <- cbind(runif(16, 0, 1000), runif(16, 0, 1000))
  g   <- sample(16, 640, TRUE)
  xy  <- ctr[g, ] + matrix(rnorm(2 * 640, sd = 20), 640)
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], w = rnorm(640)),
                      coords = c("x", "y"), crs = 32632)
  pts$z <- pts$w + rnorm(640)
  out <- determine_optimal_levels(pts, max_levels = 30, response_var = "z",
                                  predictor_vars = "w", criterion = "combined")
  d <- attr(out, "diagnostics")
  expect_false(is.null(d))
  expect_identical(d$nstart, 25L)
  expect_true(is.integer(d$wss_bumps) && d$wss_bumps >= 0L)
  expect_length(d$wss_spread, length(d$wss))
  expect_true(all(d$wss_spread >= 0))
  # The geometric path stays a plain integer vector (no attribute), as
  # documented and as its callers pin.
  geo <- determine_optimal_levels(pts, max_levels = 30)
  expect_null(attributes(geo))
})


# ---------------------------------------------------------------------------
# The analytic reliability
# ---------------------------------------------------------------------------

test_that(".rbar_rect is the mean pairwise correlation over a rectangle", {
  rb <- spatialkit:::.rbar_rect
  expect_equal(rb(function(h) rep(0.4, length(h)), 10, 10), 0.4)      # exchangeable
  expect_lt(rb(function(h) exp(-h / 10), 1000, 1000), 0.01)           # far beyond the range
  expect_gt(rb(function(h) exp(-h / 100), 10, 10), 0.9)               # well inside it
  # Larger rectangles have lower mean correlation.
  f <- function(h) exp(-h / 100)
  expect_gt(rb(f, 50, 50), rb(f, 200, 200))
  expect_true(is.na(rb(f, 0, 10)))
})

test_that(".reliability_at has an interior optimum and the measured shape", {
  f <- function(h) exp(-h / 100)                    # effective range 300
  area <- 1e6; N <- 600; c0 <- 0.3; c1 <- 1
  rbV <- spatialkit:::.rbar_rect(f, 1000, 1000)
  L <- c(1, 2, 3, 4, 6, 8, 12, 16, 24, 36, 54, 66)
  r <- vapply(L, spatialkit:::.reliability_at, numeric(1), area = area,
              n_total = N, nugget = c0, psill = c1, cor_fn = f, rbar_V = rbV)
  expect_true(all(r >= 0 & r <= 1))
  expect_identical(r[1L], 0)                         # one cell: nothing between
  # Interior maximum, where the simulation put it (analytic argmax 6,
  # empirical 5 on this design), and a decline toward the support ceiling.
  expect_true(L[which.max(r)] %in% 3:12)
  expect_lt(r[length(r)], max(r))
  # Broad: within 2% of the maximum over a factor of at least three in L.
  flat <- L[r >= 0.98 * max(r)]
  expect_gte(max(flat) / min(flat), 3)
  # More noise lowers reliability everywhere.
  r_noisy <- vapply(L, spatialkit:::.reliability_at, numeric(1), area = area,
                    n_total = N, nugget = 2, psill = c1, cor_fn = f, rbar_V = rbV)
  expect_true(all(r_noisy[-1L] < r[-1L]))
})


# ---------------------------------------------------------------------------
# resolution_profile()
# ---------------------------------------------------------------------------

test_that("a geometry-only profile scores the ladder on WSS, elbow and support", {
  pts <- rp_field(300)
  prof <- resolution_profile(pts, n_levels = 8)
  expect_s3_class(prof, "resolution_profile")
  expect_s3_class(prof, "data.frame")
  expect_named(prof, c("levels", "wss", "wss_spread", "elbow", "cell_n_min",
                       "cell_n_median", "cell_diam_median", "rss", "cp",
                       "moran_i", "moran_z", "reliability"))
  b <- attr(prof, "bounds")
  expect_identical(b$floor, 2L)
  expect_identical(b$ceiling, 33L)                  # floor(300 / 9)
  expect_true(b$supported)
  expect_identical(min(prof$levels), 2L)
  expect_identical(max(prof$levels), 33L)
  expect_true(all(diff(prof$levels) > 0))
  expect_true(all(is.finite(prof$wss)))
  expect_true(all(diff(prof$wss) < 0))              # monotone with 25 restarts
  expect_true(all(prof$wss_spread >= 0))
  expect_true(all(is.na(prof$rss)) && all(is.na(prof$cp)) &&
                all(is.na(prof$moran_z)) && all(is.na(prof$reliability)))
  expect_true(all(prof$cell_n_min >= 1L))
  expect_true(all(diff(prof$cell_diam_median) < 0))  # finer cells are smaller
  expect_null(attr(prof, "variogram"))
  expect_true(is.na(attr(prof, "variable")))
  expect_identical(attr(prof, "nstart"), 25L)
  expect_identical(attr(prof, "wss_bumps"), 0L)
  expect_output(print(prof), "geometry only")
})

test_that("a response adds cp, moran_z and reliability, with the floor from the range", {
  skip_if_not_installed("gstat")
  pts <- rp_field(500)
  prof <- resolution_profile(pts, response_var = "z", predictor_vars = "w", n_levels = 10)
  expect_identical(attr(prof, "variable"), "residuals")
  vg <- attr(prof, "variogram")
  expect_false(is.null(vg))
  expect_true(is.finite(vg$nugget) && vg$nugget >= 0)
  expect_true(vg$psill > 0)
  b <- attr(prof, "bounds")
  expect_true(is.finite(b$range))
  expect_identical(b$floor, max(2L, as.integer(ceiling(b$area / b$range^2))))
  expect_identical(b$ceiling, 55L)                  # floor(500 / 9)
  expect_true(b$supported)
  expect_identical(min(prof$levels), b$floor)
  expect_true(all(is.finite(prof$rss)))
  # Partitions at different levels are not nested, so RSS need not fall at
  # every step; over the ladder as a whole finer cells fit better.
  expect_lt(prof$rss[nrow(prof)], prof$rss[1L])
  expect_equal(prof$cp, prof$rss / 500 + 2 * vg$nugget * prof$levels / 500)
  expect_true(all(is.finite(prof$reliability)))
  expect_true(all(prof$reliability > 0 & prof$reliability < 1))
  # moran_z is NA at nine cells or fewer and finite above.
  expect_true(all(is.na(prof$moran_z[prof$levels <= 9L])))
  expect_true(all(is.finite(prof$moran_z[prof$levels > 9L])))
  expect_s3_class(attr(prof, "sac"), "sac_range")
  expect_output(print(prof), "residuals")
})

test_that("a supplied sac object is used as given", {
  skip_if_not_installed("gstat")
  pts <- rp_field(300)
  vm  <- data.frame(model = c("Nug", "Exp"), psill = c(0.5, 2), range = c(0, 80),
                    stringsAsFactors = FALSE)
  sac <- structure(240, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(pts))
  prof <- resolution_profile(pts, response_var = "z", sac = sac, n_levels = 6)
  vg <- attr(prof, "variogram")
  expect_equal(vg$nugget, 0.5)
  expect_equal(vg$psill, 2)
  expect_equal(vg$range, 240)
  expect_equal(prof$cp, prof$rss / 300 + 2 * 0.5 * prof$levels / 300)
  expect_identical(attr(prof, "bounds")$floor, max(2L, as.integer(ceiling(attr(prof, "bounds")$area / 240^2))))
})

test_that("a floor above the ceiling is reported as a finding, not resolved silently", {
  skip_if_not_installed("gstat")
  # Few points and a long range: cells narrower than the range would need
  # more cells than nine points each can support.
  pts <- rp_field(60)
  vm  <- data.frame(model = c("Nug", "Exp"), psill = c(0.2, 1), range = c(0, 40),
                    stringsAsFactors = FALSE)
  sac <- structure(120, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(pts))
  lines <- capture_spatialkit_log(
    prof <- resolution_profile(pts, response_var = "z", sac = sac, n_levels = 5))
  b <- attr(prof, "bounds")
  expect_false(b$supported)
  expect_gt(b$floor, b$ceiling)
  expect_true(log_has(lines, "cannot support a tessellation"))
  expect_identical(min(prof$levels), 2L)
  expect_identical(max(prof$levels), b$ceiling)
  expect_output(print(prof), "not supported")
})

test_that("an explicit ladder is honoured and filtered", {
  pts <- rp_field(200)
  prof <- resolution_profile(pts, levels = c(3, 5, 5, 8, 1, 400))
  expect_identical(prof$levels, c(3L, 5L, 8L))
  expect_error(resolution_profile(pts, levels = c(1, 1000)), "no usable value")
})

test_that("resolution_profile validates its input", {
  pts <- rp_field(100)
  expect_error(resolution_profile(sf::st_drop_geometry(pts)), "must be an sf object")
  expect_error(resolution_profile(pts, min_cell_n = 0), "min_cell_n")
  expect_error(resolution_profile(pts, nstart = 0), "nstart")
  expect_error(resolution_profile(pts, predictor_vars = "w"), "needs a `response_var`")
  expect_error(resolution_profile(pts, response_var = "nope"), "not found")
  pts$f <- factor(sample(letters[1:3], 100, TRUE))
  expect_error(resolution_profile(pts, response_var = "f"), "numeric or logical")
  expect_error(resolution_profile(pts, response_var = "z", predictor_vars = "f"),
               "numeric or logical")
})

test_that("resolution_profile is reproducible from its seed and subsamples large layers", {
  pts <- rp_field(200)
  a <- resolution_profile(pts, n_levels = 5, seed = 11)
  b <- resolution_profile(pts, n_levels = 5, seed = 11)
  expect_identical(a$wss, b$wss)
  sub <- resolution_profile(pts, n_levels = 5, sample_n = 120)
  expect_identical(attr(sub, "bounds")$n, 120L)
  expect_identical(attr(sub, "bounds")$ceiling, 13L)
})


# ---------------------------------------------------------------------------
# select_resolution()
# ---------------------------------------------------------------------------

test_that("select_resolution reads the optimum and the flat region off each criterion", {
  skip_if_not_installed("gstat")
  pts <- rp_field(500)
  prof <- resolution_profile(pts, response_var = "z", predictor_vars = "w", n_levels = 10)

  s_cp <- select_resolution(prof, "cp")
  expect_s3_class(s_cp, "resolution_selection")
  expect_identical(s_cp$best, prof$levels[which.min(prof$cp)])
  expect_true(all(prof$cp[prof$levels %in% s_cp$flat] <= min(prof$cp) * 1.02))
  expect_true(s_cp$best %in% s_cp$flat)
  expect_identical(s_cp$at_ceiling, s_cp$best == max(prof$levels))

  s_rel <- select_resolution(prof, "reliability")
  expect_identical(s_rel$best, prof$levels[which.max(prof$reliability)])
  expect_true(all(prof$reliability[prof$levels %in% s_rel$flat] >= max(prof$reliability) * 0.98))
  expect_identical(s_rel$at_floor, s_rel$best == min(prof$levels))

  s_el <- select_resolution(prof, "elbow")
  expect_identical(s_el$best, prof$levels[which.max(prof$elbow)])

  s_z <- select_resolution(prof, "moran_z", tol = 0.1)
  ok <- is.finite(prof$moran_z)
  expect_identical(s_z$best, prof$levels[ok][which.min(abs(prof$moran_z[ok]))])
  expect_true(all(s_z$flat %in% prof$levels[ok]))

  # A wider tolerance never narrows the region.
  expect_true(all(s_cp$flat %in% select_resolution(prof, "cp", tol = 0.1)$flat))
  expect_output(print(s_cp), "Resolution by cp")
  if (s_cp$at_ceiling) expect_output(print(s_cp), "support ceiling")
})

test_that("select_resolution refuses a criterion that is NA everywhere, and bad input", {
  pts <- rp_field(200)
  geo <- resolution_profile(pts, n_levels = 5)
  expect_error(select_resolution(geo, "cp"), "NA at every level.*response")
  expect_error(select_resolution(geo, "reliability"), "NA at every level")
  expect_error(select_resolution(geo, "moran_z"), "NA at every level")
  expect_s3_class(select_resolution(geo, "elbow"), "resolution_selection")
  expect_error(select_resolution(geo, "elbow", tol = -1), "non-negative")
  expect_error(select_resolution(data.frame(levels = 1:3), "elbow"), "must come from")
})


test_that("print() survives a subset that no longer carries the ladder", {
  skip_if_not_installed("gstat")
  prof <- suppressWarnings(suppressMessages(
    resolution_profile(rp_field(n = 250), response_var = "z", n_levels = 8)))

  expect_output(print(prof), "^Resolution profile: ")
  # A row subset keeps the bounds attribute and still summarises.  A column
  # subset loses the attribute and the `levels` column both, and used to abort
  # on is.finite(NULL) in the ladder line.
  expect_output(print(prof[1:3, ]), "^Resolution profile: 3 levels")
  expect_error(utils::capture.output(print(prof[, 1:3])), NA)
  expect_output(print(prof[, 1:3]), "subset")
})

