# tests/testthat/test-core-partition-properties.R
# ---------------------------------------------------------------------------
# What the density kind promises: cells of about equal counts on clustered
# points once the cells outnumber the clusters (and counts like plain
# k-means' below that), and a result that depends on the points alone, not
# on the order of their rows, the unit of their coordinates, repeated
# locations, the count asked for, the caller's random-number stream or
# whether FNN is installed.
# ---------------------------------------------------------------------------

.pp_engine <- function(xy, k, ..., seed = 1L) {
  set.seed(seed)
  spatialkit:::.partition_centres(xy, k, partition = "density", ...)
}

.pp_ok <- function(fit, xy, k) {
  # The shape every engine result has.
  expect_identical(dim(fit$centers), c(as.integer(k), 2L))
  expect_true(all(is.finite(fit$centers)))
  expect_length(fit$cluster, nrow(xy))
  expect_true(all(fit$cluster >= 1L & fit$cluster <= k))
  expect_identical(sum(fit$size), nrow(xy))
  expect_identical(fit$size, tabulate(fit$cluster, nbins = k))
  expect_true(all(is.finite(fit$withinss)) && all(fit$withinss >= 0))
  expect_true(is.finite(fit$objective) && fit$objective >= 0)
  expect_true(all(is.finite(fit$weights)) && all(fit$weights > 0))
  invisible(fit)
}


# ---- Counts ------------------------------------------------------------------

test_that("on a clustered layout the density kind's counts are close to equal, k-means' are not", {
  # The measure is the coefficient of variation of the cell counts (sd over
  # mean; 0 when every cell holds the same number) and the fullest cell's
  # count over the emptiest's.  1,200 points, 900 of them in four clusters
  # over a sparse background.  Simulated over 20 such layers, 5 starts,
  # medians (and ranges) of the coefficient of variation:
  #     cells   k-means             density
  #       20    1.36 (1.09-1.53)    0.31 (0.18-0.43)
  #       60    1.20 (1.04-1.40)    0.28 (0.24-0.34)
  #      130    1.18 (1.10-1.36)    0.34 (0.30-0.37)
  # and of fullest over emptiest, 24, 40 and 50 against 3.3, 4.3 and 6.8.
  for (seed in c(41L, 52L)) {
    xy <- partition_clustered_xy(seed = seed)
    for (k in c(20L, 60L, 130L)) {
      set.seed(1)
      km <- spatialkit:::.partition_centres(partition_sorted(xy), k, partition = "kmeans",
                                            init = "kmeans++", nstart = 5L)
      den <- .pp_engine(xy, k, nstart = 5L, iter.max = 100L)
      cv_km <- partition_count_cv(km$size); cv_den <- partition_count_cv(den$size)
      expect_gt(cv_km, 0.9)
      expect_lt(cv_den, 0.5)
      expect_lt(cv_den, cv_km / 2.5)
      expect_gt(max(km$size) / min(km$size), 3 * max(den$size) / min(den$size))
      expect_gt(min(den$size), 0L)
    }
  }
})

test_that("the counts even out only once the cells outnumber the dense places", {
  # 2,000 points with 28 dense places: half of the points in 3 large
  # clusters, 30 percent in 25 small tight ones, the rest spread out.
  # Simulated over 20 such layers, 5 starts, medians (and ranges) of the
  # coefficient of variation of the cell counts, and of the density kind's
  # over k-means':
  #     cells   k-means             density             density over k-means
  #        5    0.63 (0.41-0.93)    0.71 (0.40-0.92)    1.01 (0.69-1.37)
  #       10    0.90 (0.79-1.38)    0.97 (0.80-1.24)    1.00 (0.88-1.49)
  #       20    1.08 (0.85-1.32)    1.01 (0.87-1.47)    0.98 (0.78-1.32)
  #       60    0.95 (0.86-1.27)    0.31 (0.24-0.35)    0.31 (0.21-0.37)
  #      130    1.05 (0.92-1.12)    0.43 (0.37-0.45)    0.41 (0.36-0.47)
  # With fewer cells than dense places a cell takes in whole clusters, large
  # and small, under either kind, and the density kind's counts are as
  # uneven as k-means'.  On the four-cluster layout of the test above every
  # count tried is above the number of clusters.
  for (seed in c(1L, 2L)) {
    xy <- partition_towns_xy(seed = seed)
    cv_at <- function(k) {
      set.seed(1)
      km <- spatialkit:::.partition_centres(partition_sorted(xy), k, partition = "kmeans",
                                            init = "kmeans++", nstart = 5L)
      den <- .pp_engine(xy, k, nstart = 5L, iter.max = 100L)
      c(km = partition_count_cv(km$size), den = partition_count_cv(den$size))
    }
    # Measured on these two layers: 1.12 and 1.05 at 10 cells, 0.87 and 1.17
    # at 20, between 0.89 and 1.31 of k-means'.
    for (k in c(10L, 20L)) {
      v <- cv_at(k)
      expect_gt(v[["den"]], 0.7)
      expect_gt(v[["den"]], 0.75 * v[["km"]])
    }
    # And 0.25 and 0.31 at 60 cells, 0.45 and 0.41 at 130, under half of
    # k-means' 0.88 to 0.97.
    for (k in c(60L, 130L)) {
      v <- cv_at(k)
      expect_gt(v[["km"]], 0.7)
      expect_lt(v[["den"]], 0.55)
      expect_lt(v[["den"]], 0.6 * v[["km"]])
    }
  }
})

test_that("the public functions show both ends on that layout", {
  xy  <- partition_towns_xy(seed = 1L)
  pts <- partition_points_sf(xy)
  count_cv <- function(seeds) partition_count_cv(
    tabulate(partition_nearest(xy, partition_xy(seeds)), nbins = nrow(seeds)))
  # Measured: 0.82 at 10 cells under the density kind, as under k-means from
  # either start; 0.31 at 60 cells, against 0.88 for k-means from k-means++
  # starts.
  expect_gt(count_cv(voronoi_seeds_kmeans(pts, 10, partition = "density")), 0.7)
  expect_gt(count_cv(voronoi_seeds_kmeans(pts, 10, init = "kmeans++")), 0.7)
  expect_lt(count_cv(voronoi_seeds_kmeans(pts, 60, partition = "density")), 0.5)
  expect_gt(count_cv(voronoi_seeds_kmeans(pts, 60, init = "kmeans++")), 0.7)
  # In a profile the median cell over the mean cell (1 for equal counts) was
  # 0.42 and 0.42 at 10 and 20 cells under the density kind, where k-means
  # gave 0.51 and 0.42, and 0.81 and 0.72 at 60 and 130 cells.
  lv <- c(10L, 20L, 60L, 130L)
  p_den <- resolution_profile(pts, levels = lv, nstart = 4, partition = "density")
  share <- p_den$cell_n_median / (nrow(pts) / lv)
  expect_true(all(share[1:2] < 0.55))
  expect_true(all(share[3:4] > 0.65))
})

test_that("the public functions give those counts", {
  xy  <- partition_clustered_xy(seed = 41L)
  pts <- partition_points_sf(xy)
  k <- 60L
  count_cv <- function(seeds) partition_count_cv(
    tabulate(partition_nearest(xy, partition_xy(seeds)), nbins = nrow(seeds)))
  # Measured on this layer: 0.74 for the function's own random starts, 1.12
  # for k-means++ starts (the lower sum of squares, and the less even
  # counts), 0.32 for the density kind and 0.34 for it from random starts.
  cv_default <- count_cv(voronoi_seeds_kmeans(pts, k))
  cv_pp      <- count_cv(voronoi_seeds_kmeans(pts, k, init = "kmeans++"))
  cv_den     <- count_cv(voronoi_seeds_kmeans(pts, k, partition = "density"))
  cv_den_rnd <- count_cv(voronoi_seeds_kmeans(pts, k, partition = "density", init = "random"))
  expect_gt(cv_default, 0.6)
  expect_gt(cv_pp, 0.9)
  expect_lt(cv_den, 0.5)
  expect_lt(cv_den_rnd, 0.5)
  expect_lt(cv_den, cv_default / 2)
  expect_lt(count_cv(get_voronoi_seeds(method = "kmeans", n = k, sample_points = pts,
                                       partition = "density", set_seed = 2)), 0.5)
  # In a profile.  With equal counts the median cell holds the mean count,
  # n / L; under k-means a few large cells of sparse ground hold most of the
  # points and the median cell holds far fewer.  Measured: 0.91 and 0.95 of
  # the mean under the density kind at 20 and 60 cells, 0.34 and 0.40 under
  # k-means; smallest cells of 35 and 8 points against 10 and 2.
  lv <- c(20L, 60L)
  p_km  <- resolution_profile(pts, levels = lv, nstart = 4)
  p_den <- resolution_profile(pts, levels = lv, nstart = 4, partition = "density")
  expect_true(all(p_den$cell_n_median / (nrow(pts) / lv) > 0.8))
  expect_true(all(p_km$cell_n_median / (nrow(pts) / lv) < 0.5))
  expect_true(all(p_den$cell_n_min >= 2L * p_km$cell_n_min))
})

test_that("density_power runs from plain k-means' counts at 0 to equal counts at 1", {
  xy <- partition_clustered_xy(seed = 41L)
  cv_at <- function(p) partition_count_cv(.pp_engine(xy, 60L, nstart = 4L, iter.max = 100L,
                                                     density_power = p)$size)
  cv <- vapply(c(0, 0.5, 1), cv_at, numeric(1))
  expect_gt(cv[1L], 0.9)                 # every weight is one
  expect_lt(cv[2L], cv[1L])
  expect_lt(cv[3L], cv[2L])
  # Power 0 is plain k-means by Lloyd's iteration: the weights are all one
  # and the weighted sum is the within-cell sum of squares about the centres,
  # which are then the cell means.
  f0 <- .pp_engine(xy, 25L, nstart = 3L, iter.max = 200L, density_power = 0)
  expect_identical(f0$weights, rep(1, nrow(xy)))
  expect_true(f0$converged)
  expect_equal(f0$objective, f0$wss, tolerance = 1e-10)
  expect_equal(f0$centers, rowsum(xy, f0$cluster) / f0$size, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("withinss is the unweighted sum of squares of each cell about its own mean", {
  xy  <- partition_clustered_xy(seed = 41L)
  fit <- .pp_ok(.pp_engine(xy, 18L, nstart = 3L), xy, 18L)
  mu  <- rowsum(xy, fit$cluster) / fit$size
  ss  <- as.numeric(rowsum(rowSums((xy - mu[fit$cluster, ])^2), fit$cluster))
  expect_equal(fit$withinss, ss, tolerance = 1e-10)
  expect_equal(fit$wss, sum(ss), tolerance = 1e-10)
  # Between and within add up to the total, as for any partition into cells.
  tss <- sum(sweep(xy, 2L, colMeans(xy))^2)
  bss <- sum(fit$size * rowSums(sweep(mu, 2L, colMeans(xy))^2))
  expect_equal(fit$wss + bss, tss, tolerance = 1e-10)
  # The centres are weighted means, so they are not the cell means, and the
  # cells are those of the centres: every point is with its nearest one.
  expect_gt(max(abs(fit$centers - mu)), 1e-3)
  expect_identical(fit$cluster, as.integer(partition_nearest(xy, fit$centers)))
  wm <- rowsum(xy * fit$weights, fit$cluster) / as.numeric(rowsum(fit$weights, fit$cluster))
  expect_true(fit$converged)
  expect_equal(fit$centers, wm, tolerance = 1e-9, ignore_attr = TRUE)
})


# ---- Row order ----------------------------------------------------------------

test_that("the density kind does not depend on the order of the rows", {
  xy <- partition_clustered_xy(n_bg = 150L, n_cl = 350L, seed = 7L)
  set.seed(99); perm <- sample.int(nrow(xy))
  for (init in c("kmeans++", "random")) for (k in c(1L, 6L, 40L)) {
    a <- .pp_engine(xy, k, nstart = 3L, init = init)
    b <- .pp_engine(xy[perm, ], k, nstart = 3L, init = init)
    expect_identical(b$centers, a$centers)
    expect_identical(b$cluster, a$cluster[perm])
    expect_identical(b$weights, a$weights[perm])
    expect_identical(b[c("size", "withinss", "wss", "objective", "spread", "iter")],
                     a[c("size", "withinss", "wss", "objective", "spread", "iter")])
  }
  # Through the public functions.
  pts <- partition_points_sf(xy)
  s1 <- voronoi_seeds_kmeans(pts, 25, partition = "density")
  s2 <- voronoi_seeds_kmeans(pts[perm, ], 25, partition = "density")
  expect_identical(partition_xy(s2), partition_xy(s1))
  g1 <- get_voronoi_seeds(method = "kmeans", n = 25, sample_points = pts,
                          partition = "density", set_seed = 3)
  g2 <- get_voronoi_seeds(method = "kmeans", n = 25, sample_points = pts[perm, ],
                          partition = "density", set_seed = 3)
  expect_identical(partition_xy(g2), partition_xy(g1))
  expect_identical(attr(g2, "kmeans")$cluster, attr(g1, "kmeans")$cluster[perm])
  p1 <- resolution_profile(pts, n_levels = 5, nstart = 3, partition = "density")
  p2 <- resolution_profile(pts[perm, ], n_levels = 5, nstart = 3, partition = "density")
  expect_identical(as.data.frame(p2), as.data.frame(p1))
  expect_identical(attr(p2, "centres"), attr(p1, "centres"))
  # k-means++ starts in a seeding function are order-free as well; the
  # function's own random starts cluster the rows as they come.
  expect_identical(partition_xy(voronoi_seeds_kmeans(pts[perm, ], 25, init = "kmeans++")),
                   partition_xy(voronoi_seeds_kmeans(pts, 25, init = "kmeans++")))
})


# ---- Units ----------------------------------------------------------------------

test_that("the density kind does not depend on the unit of the coordinates", {
  xy <- partition_clustered_xy(n_bg = 150L, n_cl = 350L, seed = 7L)
  for (k in c(5L, 30L)) {
    m <- .pp_engine(xy, k, nstart = 3L)
    for (f in c(1e-3, 3937 / 1200, 1e6)) {          # kilometres, US survey feet, micrometres
      u <- .pp_engine(xy * f, k, nstart = 3L)
      expect_identical(u$cluster, m$cluster)
      expect_identical(u$iter, m$iter)
      expect_lt(max(abs(u$centers / f - m$centers)), 1e-9 * 1000)
      expect_equal(u$weights, m$weights, tolerance = 1e-9)
      expect_equal(u$wss / f^2, m$wss, tolerance = 1e-9)
      expect_equal(u$objective / f^2, m$objective, tolerance = 1e-9)
      expect_equal(u$spread, m$spread, tolerance = 1e-9)
    }
    # A shift of origin moves the centres with the points.
    s <- .pp_engine(sweep(xy, 2L, c(5e5, 5e6), "+"), k, nstart = 3L)
    expect_identical(s$cluster, m$cluster)
    expect_lt(max(abs(sweep(s$centers, 2L, c(5e5, 5e6)) - m$centers)), 1e-6)
  }
})

test_that("seeds are the same places whether the layer is in metres or in feet", {
  # One projection in two units: NAD83 / North Carolina in metres (EPSG:32119)
  # and in US survey feet (EPSG:2264).
  xy <- partition_clustered_xy(n_bg = 120L, n_cl = 280L, side = 20000,
                               origin = c(6e5, 2e5), seed = 9L)
  m  <- partition_points_sf(xy, crs = 32119)
  ft <- tryCatch(sf::st_transform(m, 2264), error = function(e) NULL)
  skip_if(is.null(ft), "EPSG:2264 is not available in this PROJ installation")
  ratio <- stats::median(sf::st_coordinates(ft)[, 1L] / sf::st_coordinates(m)[, 1L])
  skip_if_not(isTRUE(all.equal(ratio, 3937 / 1200, tolerance = 1e-6)),
              "the two CRSs do not differ by the unit alone here")
  s_m  <- voronoi_seeds_kmeans(m, 20, partition = "density")
  s_ft <- voronoi_seeds_kmeans(ft, 20, partition = "density")
  back <- partition_xy(sf::st_transform(s_ft, 32119))
  expect_lt(max(abs(back - partition_xy(s_m))), 1e-4)       # a tenth of a millimetre
})


# ---- Repeated locations ----------------------------------------------------------

test_that("repeated locations are counted, share a cell and never give a zero distance", {
  # Five stations visited thirty times each, and 200 single visits.
  cleanup <- spatialkit:::.with_seed(8L)
  st <- cbind(stats::runif(5, 100, 900), stats::runif(5, 100, 900))
  sc <- cbind(stats::runif(200, 0, 1000), stats::runif(200, 0, 1000))
  cleanup()
  xy <- rbind(st[rep(1:5, each = 30L), ], sc)
  station <- c(rep(1:5, each = 30L), rep(NA_integer_, 200L))
  for (k in c(1L, 5L, 12L, 30L, 205L)) {
    fit <- .pp_ok(.pp_engine(xy, k, nstart = 3L), xy, k)
    # The thirty visits to a station are in one cell.
    per_station <- tapply(fit$cluster[1:150], station[1:150], function(v) length(unique(v)))
    expect_true(all(per_station == 1L))
    # And carry one weight, which is finite and of ordinary size: a station
    # is thirty points in a disc that reaches its eighth other location, not
    # a point with eight neighbours at distance zero.
    w_station <- tapply(fit$weights[1:150], station[1:150], function(v) length(unique(v)))
    expect_true(all(w_station == 1L))
    expect_lt(max(fit$weights) / min(fit$weights), 1e4)
  }
  # One centre per distinct location is the most there can be.
  expect_identical(nrow(unique(xy)), 205L)
  full <- .pp_engine(xy, 205L, nstart = 2L)
  expect_identical(sort(full$size, decreasing = TRUE)[1:6], c(30L, 30L, 30L, 30L, 30L, 1L))
  expect_equal(full$wss, 0, tolerance = 1e-12)
  expect_null(.pp_engine(xy, 206L, nstart = 2L))
  # The rows in another order give the same cells.
  set.seed(5); perm <- sample.int(nrow(xy))
  expect_identical(.pp_engine(xy[perm, ], 12L, nstart = 3L)$centers,
                   .pp_engine(xy, 12L, nstart = 3L)$centers)
})

test_that("layers that are all repeats, one location or a line are handled", {
  st <- cbind(c(0, 10, 10, 0, 5), c(0, 0, 10, 10, 5))
  xy <- st[rep(1:5, times = c(30L, 30L, 30L, 30L, 30L)), ]
  for (k in 1:5) {
    fit <- .pp_ok(.pp_engine(xy, k, nstart = 3L), xy, k)
    expect_true(all(fit$size %% 30L == 0L))
  }
  expect_identical(sort(.pp_engine(xy, 5L)$size), rep(30L, 5L))
  expect_null(.pp_engine(xy, 6L))
  # Every point at one location: one centre, there.
  one <- matrix(c(3, 4), 20L, 2L, byrow = TRUE)
  f1 <- .pp_ok(.pp_engine(one, 1L), one, 1L)
  expect_identical(f1$centers, matrix(c(3, 4), 1L))
  expect_identical(f1$wss, 0)
  expect_null(.pp_engine(one, 2L))
  # Points on a line, and two points.
  line <- cbind(seq(0, 100, length.out = 40L), 7)
  fl <- .pp_ok(.pp_engine(line, 6L, nstart = 3L), line, 6L)
  expect_equal(fl$centers[, 2L], rep(7, 6L), tolerance = 1e-12)
  two <- rbind(c(0, 0), c(1, 1))
  expect_identical(.pp_ok(.pp_engine(two, 2L), two, 2L)$size, c(1L, 1L))
  expect_equal(.pp_engine(two, 1L)$centers, matrix(c(0.5, 0.5), 1L), tolerance = 1e-12)
  # A non-finite coordinate is an error, not a non-finite centre.
  expect_error(.pp_engine(rbind(two, c(NA, 1)), 1L), "non-finite coordinates")
})


# ---- The count asked for ---------------------------------------------------------

test_that("k = 1, k = n - 1 and k = n give that many centres, and k above n none", {
  xy <- partition_clustered_xy(n_bg = 30L, n_cl = 50L, seed = 12L)
  n  <- nrow(xy)
  for (init in c("kmeans++", "random")) {
    f1 <- .pp_ok(.pp_engine(xy, 1L, init = init), xy, 1L)
    # One centre: the weighted mean of the points.
    expect_equal(as.numeric(f1$centers), colSums(xy * f1$weights) / sum(f1$weights),
                 tolerance = 1e-12)
    expect_identical(f1$size, n)
    fn1 <- .pp_ok(.pp_engine(xy, n - 1L, init = init, iter.max = 100L), xy, n - 1L)
    expect_identical(sort(fn1$size), c(rep(1L, n - 2L), 2L))
    fn <- .pp_ok(.pp_engine(xy, n, init = init), xy, n)
    expect_identical(fn$size, rep(1L, n))
    expect_equal(fn$wss, 0, tolerance = 1e-12)
    expect_equal(fn$objective, 0, tolerance = 1e-12)
    expect_null(.pp_engine(xy, n + 1L, init = init))
    expect_null(.pp_engine(xy, 0L, init = init))
  }
})

test_that("the public functions clamp a count at or above the number of points", {
  xy  <- partition_clustered_xy(n_bg = 30L, n_cl = 50L, seed = 12L)
  pts <- partition_points_sf(xy)
  n   <- nrow(pts)
  expect_identical(nrow(voronoi_seeds_kmeans(pts, 1, partition = "density")), 1L)
  expect_identical(nrow(voronoi_seeds_kmeans(pts, n - 1L, partition = "density")), n - 1L)
  for (k in c(n, n + 50L)) {
    expect_warning(s <- voronoi_seeds_kmeans(pts, k, partition = "density"), "clamping to 79")
    expect_identical(nrow(s), n - 1L)
    expect_true(all(is.finite(partition_xy(s))))
    expect_warning(g <- get_voronoi_seeds(method = "kmeans", n = k, sample_points = pts,
                                          partition = "density", set_seed = 1),
                   "clamping")
    expect_identical(nrow(g), n - 1L)
    expect_identical(sum(attr(g, "kmeans")$size), n)
  }
  # A layer of repeats is held to its distinct locations.
  dup <- pts[rep(1:10, each = 6L), ]
  expect_warning(s <- voronoi_seeds_kmeans(dup, 25, partition = "density"),
                 "only 10 unique positions")
  expect_identical(nrow(s), 10L)
  # A density profile's ladder stops where a k-means one does.
  p <- resolution_profile(pts, n_levels = 4, nstart = 2, min_cell_n = 1, partition = "density")
  expect_identical(max(p$levels), n - 1L)
  expect_identical(p$levels, resolution_profile(pts, n_levels = 4, nstart = 2,
                                                min_cell_n = 1)$levels)
  expect_true(all(is.finite(p$wss)))
  expect_true(all(vapply(attr(p, "centres"), function(m) all(is.finite(m)), logical(1))))
})


# ---- Seeds and the caller's stream -------------------------------------------------

test_that("a seed fixes the result, and the caller's random-number stream is untouched", {
  xy  <- partition_clustered_xy(n_bg = 100L, n_cl = 200L, seed = 13L)
  pts <- partition_points_sf(xy)
  set.seed(2024); before <- .Random.seed; expected_next <- stats::runif(3)

  check <- function(call) {
    set.seed(2024)
    a <- call()
    expect_identical(.Random.seed, before)
    expect_identical(stats::runif(3), expected_next)
    set.seed(77)                                   # another stream, the same answer
    expect_identical(call(), a)
    a
  }
  v <- check(function() voronoi_seeds_kmeans(pts, 15, partition = "density"))
  check(function() voronoi_seeds_kmeans(pts, 15, init = "kmeans++"))
  check(function() get_voronoi_seeds(method = "kmeans", n = 15, sample_points = pts,
                                     partition = "density", set_seed = 5))
  check(function() resolution_profile(pts, n_levels = 4, nstart = 2, partition = "density"))
  check(function() resolution_profile(pts, n_levels = 4, nstart = 2, init = "random"))
  # Another seed is another partition.
  expect_false(identical(partition_xy(voronoi_seeds_kmeans(pts, 15, partition = "density",
                                                           set_seed = 1)),
                         partition_xy(v)))
  # set_seed = NULL draws from the session's stream, as it always did.
  set.seed(31); a <- get_voronoi_seeds(method = "kmeans", n = 15, sample_points = pts,
                                       partition = "density")
  set.seed(31); b <- get_voronoi_seeds(method = "kmeans", n = 15, sample_points = pts,
                                       partition = "density")
  expect_identical(a, b)
  # The engine itself draws from the stream it is given, and nowhere else.
  set.seed(8); e1 <- spatialkit:::.partition_centres(xy, 9L, partition = "density")
  set.seed(8); e2 <- spatialkit:::.partition_centres(xy, 9L, partition = "density")
  expect_identical(e1, e2)
})


# ---- Nothing non-finite ---------------------------------------------------------------

test_that("no input gives a non-finite centre", {
  cleanup <- spatialkit:::.with_seed(3L)
  awkward <- list(
    lattice   = as.matrix(expand.grid(x = 1:15, y = 1:15)),
    far       = cbind(1e9 + stats::runif(200), -4e8 + stats::runif(200)),
    tiny      = cbind(stats::runif(200), stats::runif(200)) * 1e-9,
    huge      = cbind(stats::runif(200), stats::runif(200)) * 1e150,
    sliver    = cbind(stats::runif(200, 0, 1e6), stats::runif(200, 0, 1e-3)),
    two_scale = rbind(cbind(stats::runif(150), stats::runif(150)) * 1e-4,
                      cbind(stats::runif(50, 0, 1e4), stats::runif(50, 0, 1e4))),
    near_dup  = rbind(cbind(stats::runif(100), stats::runif(100)),
                      matrix(0.5, 30L, 2L) + stats::rnorm(60, sd = 1e-13)))
  cleanup()
  for (nm in names(awkward)) {
    xy <- unname(awkward[[nm]])
    for (k in c(1L, 3L, 20L, 90L)) for (p in c(0, 1, 2)) for (init in c("kmeans++", "random")) {
      fit <- .pp_engine(xy, k, nstart = 2L, density_power = p, init = init)
      expect_false(is.null(fit))
      .pp_ok(fit, xy, k)
      bb <- apply(xy, 2L, range)
      # The centres lie within the points' bounding box.
      slack <- 1e-9 * max(bb[2L, ] - bb[1L, ])
      expect_true(all(fit$centers[, 1L] >= bb[1L, 1L] - slack &
                        fit$centers[, 1L] <= bb[2L, 1L] + slack))
      expect_true(all(fit$centers[, 2L] >= bb[1L, 2L] - slack &
                        fit$centers[, 2L] <= bb[2L, 2L] + slack))
    }
  }
  # density_k from one neighbour to more than there are points.
  xy <- awkward$lattice
  for (kn in c(1L, 3L, 8L, 224L, 5000L))
    .pp_ok(.pp_engine(xy, 12L, nstart = 2L, density_k = kn), xy, 12L)
})


# ---- FNN ------------------------------------------------------------------------------

test_that("the result is the same with and without FNN, ties included", {
  skip_if_not_installed("FNN")
  cleanup <- spatialkit:::.with_seed(4L)
  lattice <- as.matrix(expand.grid(x = 1:20, y = 1:20))
  sets <- list(
    clustered = partition_clustered_xy(n_bg = 150L, n_cl = 350L, seed = 7L),
    # On a lattice a point is often the same distance from two centres, and
    # several locations are the same distance from a point.
    lattice   = lattice,
    repeats   = rbind(lattice, lattice[sample.int(400L, 150L), ]))
  cleanup()
  for (nm in names(sets)) {
    xy <- unname(sets[[nm]])
    for (k in c(4L, 16L, 60L)) for (kn in c(5L, 8L)) for (init in c("kmeans++", "random")) {
      with_fnn <- .pp_engine(xy, k, nstart = 2L, density_k = kn, init = init, use_fnn = TRUE)
      without  <- .pp_engine(xy, k, nstart = 2L, density_k = kn, init = init, use_fnn = FALSE)
      expect_identical(without, with_fnn)
    }
  }
  # The pieces on their own: the neighbour sets and distances, and the
  # nearest centre under ties.
  u <- spatialkit:::.partition_locations(sets$repeats)$u
  a <- spatialkit:::.knn_locations(u, 5L, use_fnn = TRUE)
  b <- spatialkit:::.knn_locations(u, 5L, use_fnn = FALSE)
  expect_identical(a$d, b$d)
  expect_identical(t(apply(a$idx, 1L, sort)), t(apply(b$idx, 1L, sort)))
  cx <- c(5, 10, 15, 10); cy <- c(10, 5, 10, 15)
  expect_identical(spatialkit:::.nearest_centre(lattice[, 1L], lattice[, 2L], cx, cy, use_fnn = TRUE),
                   spatialkit:::.nearest_exact(lattice[, 1L], lattice[, 2L], cx, cy))
  # The point in the middle is the same distance from all four: the lowest
  # centre number has it.
  mid <- which(lattice[, 1L] == 10 & lattice[, 2L] == 10)
  expect_identical(spatialkit:::.nearest_centre(lattice[, 1L], lattice[, 2L], cx, cy)[mid], 1L)
})

test_that("tied neighbour distances are settled by the kd-tree as the exact search settles them", {
  skip_if_not_installed("FNN")
  same <- function(u, k) {
    a <- spatialkit:::.knn_locations(u, k, use_fnn = TRUE)
    b <- spatialkit:::.knn_locations(u, k, use_fnn = FALSE)
    expect_identical(a$d, b$d)
    expect_identical(t(apply(a$idx, 1L, sort)), t(apply(b$idx, 1L, sort)))
    # The tie search on its own, asked for every row.
    tied <- spatialkit:::.knn_tied_rows(u, seq_len(nrow(u)), k)
    expect_identical(t(apply(tied, 1L, sort)), t(apply(b$idx, 1L, sort)))
    invisible(a)
  }
  # Layers on which most rows tie: a square lattice at projected magnitudes,
  # a triangular one, a line, and coordinates rounded to a coarse grid.
  square <- as.matrix(expand.grid(x = 1:24, y = 1:24)) * 100 + 4e5
  tri <- cbind(rep(1:20, 20L) + rep(c(0, 0.5), each = 20L, times = 10L),
               rep(1:20, each = 20L) * sqrt(3) / 2)
  line <- cbind(1:200, 0)
  rounded <- spatialkit:::.partition_locations(
    round(partition_clustered_xy(n_bg = 300L, n_cl = 900L, seed = 3L) / 10) * 10)$u
  for (k in c(1L, 2L, 3L, 6L, 8L, 13L, 40L)) {
    same(square, k); same(tri, k); same(line, k); same(rounded, k)
  }
  same(square, nrow(square) - 2L)
  same(square, nrow(square) - 1L)
  # Twelve locations at distance exactly 5 from the origin (3-4-5 triangles)
  # and a far ring.  For one neighbour the first shortlist, four long, lies
  # inside the tie, and so does the second; the origin's k nearest are the
  # first k of the twelve in coordinate order.
  ring <- rbind(c(5, 0), c(-5, 0), c(0, 5), c(0, -5), c(3, 4), c(3, -4), c(-3, 4), c(-3, -4),
                c(4, 3), c(4, -3), c(-4, 3), c(-4, -3))
  u <- spatialkit:::.partition_locations(rbind(c(0, 0), ring, ring * 7 + 0.25))$u
  origin <- which(u[, 1L] == 0 & u[, 2L] == 0)
  on_ring <- which(rowSums(u^2) == 25)
  expect_length(on_ring, 12L)
  for (k in c(1L, 2L, 3L, 5L, 11L, 12L, 13L, 24L)) {
    a <- same(u, k)
    if (k <= 12L) {
      expect_identical(sort(a$idx[origin, ]), on_ring[seq_len(k)])
      expect_identical(a$d[origin], 5)
    }
  }
  # In any row order: the rule is the lower row number, whatever order the
  # kd-tree meets the tied locations in.  36 locations at distance exactly
  # 65 from the origin, in shuffled rows, so that a shortlist that ends
  # inside the tie seldom holds the lowest row numbers and has to be
  # lengthened until it holds every tied location.
  p65  <- rbind(c(16, 63), c(25, 60), c(33, 56), c(39, 52))
  oct  <- rbind(p65, p65[, 2:1])
  r65  <- rbind(oct, cbind(-oct[, 1L], oct[, 2L]), cbind(oct[, 1L], -oct[, 2L]), -oct,
                c(65, 0), c(-65, 0), c(0, 65), c(0, -65))
  expect_identical(unique(rowSums(r65^2)), 4225)
  all65 <- rbind(c(0, 0), r65, r65 * 3 + 0.5)
  for (seed in 1:4) {
    cleanup <- spatialkit:::.with_seed(seed)
    u <- all65[sample.int(nrow(all65)), ]
    cleanup()
    origin <- which(u[, 1L] == 0 & u[, 2L] == 0)
    tied <- which(rowSums(u^2) == 4225)
    expect_length(tied, 36L)
    for (k in c(1L, 2L, 3L, 7L, 20L)) {
      a <- same(u, k)
      expect_identical(sort(a$idx[origin, ]), tied[seq_len(k)])
    }
  }
  # With FNN no row is compared with every location, tied or not: the exact
  # search is stood in for by an error, and the tied lattice still goes
  # through.
  expected <- spatialkit:::.knn_locations(square, 6L, use_fnn = FALSE)
  local_mocked_bindings(.knn_exact_rows = function(...) stop("the exact search was called"))
  got <- spatialkit:::.knn_locations(square, 6L, use_fnn = TRUE)
  expect_identical(got$d, expected$d)
  expect_error(spatialkit:::.knn_locations(square, 6L, use_fnn = FALSE),
               "the exact search was called")
})

test_that("the floor on the neighbour distance is a millionth of the extent, in any unit", {
  # A lattice of spacing 100 with two close pairs of locations: one pair a
  # ten-thousandth of the extent apart, which is above the floor, and one a
  # hundred-millionth apart, which is read as a millionth.  With one
  # neighbour the local density at a location is 1 / (pi d^2), d the distance
  # to the nearest other location.
  lat <- as.matrix(expand.grid(x = 0:5, y = 0:5)) * 100
  extent <- sqrt(2) * 500
  xy <- unname(rbind(lat,
                     c(150.5, 250.5), c(150.5 + 1e-4 * extent, 250.5),
                     c(350.5, 250.5), c(350.5 + 1e-8 * extent, 250.5)))
  for (f in c(1, 1000, 1e-3)) {
    fit <- .pp_engine(xy * f, 3L, density_k = 1L, nstart = 1L)
    w <- fit$weights
    # The corner of the lattice has its nearest neighbour at 100.
    expect_equal(w[37:38] / w[1L], rep((100 / (1e-4 * extent))^2, 2L), tolerance = 1e-6)
    expect_equal(w[39:40] / w[1L], rep((100 / (1e-6 * extent))^2, 2L), tolerance = 1e-9)
  }
  # The function that reads the densities, on coordinates already divided by
  # the extent.
  u <- cbind(c(0, 0.3, 0.3 + 1e-4, 0.7, 0.7 + 1e-8, 1), c(0, 0.5, 0.5, 0.5, 0.5, 1) / 1)
  w <- spatialkit:::.density_weights(u, rep(1L, 6L), power = 1, k_nn = 1L, use_fnn = FALSE)
  d_corner <- sqrt(0.3^2 + 0.5^2)
  expect_equal(w[2:3] / w[1L], rep((d_corner / 1e-4)^2, 2L), tolerance = 1e-6)
  expect_equal(w[4:5] / w[1L], rep((d_corner / 1e-6)^2, 2L), tolerance = 1e-9)
  expect_identical(spatialkit:::.DENSITY_FLOOR_REL, 1e-6)
})

test_that("the local density is the usual nearest-neighbour estimate where no location repeats", {
  cleanup <- spatialkit:::.with_seed(6L)
  xy <- cbind(stats::runif(120), stats::runif(120))
  cleanup()
  loc <- spatialkit:::.partition_locations(xy)
  expect_identical(loc$m, rep(1L, 120L))
  d <- as.matrix(stats::dist(loc$u)); diag(d) <- Inf
  d8 <- unname(apply(d, 1L, function(v) sort(v)[8L]))
  w <- spatialkit:::.density_weights(loc$u, loc$m, power = 1, k_nn = 8L, use_fnn = FALSE)
  expect_equal(w, (8 / (pi * d8^2)) / mean(8 / (pi * d8^2)), tolerance = 1e-12)
  # With repeats: the points at the location and at its 8 nearest other
  # locations, over the disc that reaches the furthest of those.
  m2 <- rep(1L, 120L); m2[c(3L, 50L)] <- c(30L, 4L)
  w2 <- spatialkit:::.density_weights(loc$u, m2, power = 1, k_nn = 8L, use_fnn = FALSE)
  nb <- unname(t(apply(d, 1L, function(v) order(v)[1:8])))
  cnt <- (m2 - 1) + rowSums(matrix(m2[nb], 120L, 8L))
  rho <- cnt / (pi * d8^2)
  expect_equal(w2, rho / mean(rep.int(rho, m2)), tolerance = 1e-12)
  expect_equal(mean(rep.int(w2, m2)), 1, tolerance = 1e-12)
  # Power 0, and one location: every weight is one.
  expect_identical(spatialkit:::.density_weights(loc$u, loc$m, 0, 8L, FALSE), rep(1, 120L))
  expect_identical(spatialkit:::.density_weights(loc$u[1L, , drop = FALSE], 5L, 1, 8L, FALSE), 1)
})
