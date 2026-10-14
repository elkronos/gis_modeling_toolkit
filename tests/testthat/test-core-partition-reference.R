# tests/testthat/test-core-partition-reference.R
# ---------------------------------------------------------------------------
# The package's density kind and its random starts reproduce the experiment
# that tested them.
#
# den_centres() and rnd_centres() (helper-partition-reference.R) are that
# experiment's functions.  Two of den_centres()'s guards depend on the unit
# of the coordinates: the neighbour distance is not taken below 1, and Lloyd's
# iteration stops when no centre moves by 1e-6.  The package states both as
# fractions of the extent.  On point sets where neither guard binds (every
# 8th-neighbour distance is far above 1, and a step either moves a centre by
# far more than 1e-6 or not at all), given the same point order and seed,
# the package's centres are the reference's.
#
# den_centres() takes the rows as they come; the package works in coordinate
# order, so the reference is given the rows in that order.
# ---------------------------------------------------------------------------

.pr_sets <- function() list(
  # 1,200 points, 900 of them in four clusters, on a square of side 1000.
  clustered = partition_clustered_xy(seed = 11L),
  # 800 points spread over a 2000 x 700 rectangle.
  uniform = {
    cleanup <- spatialkit:::.with_seed(5L)
    xy <- cbind(stats::runif(800, 0, 2000), stats::runif(800, 0, 700))
    cleanup()
    xy
  })

.pr_density <- function(xy, k, seed, ...) {
  set.seed(seed)
  spatialkit:::.partition_centres(xy, k, partition = "density", init = "kmeans++",
                                  nstart = DEN_NSTART, iter.max = DEN_ITER,
                                  density_power = DEN_POWER, density_k = DEN_NN, ...)
}


test_that("neither unit-dependent guard of the reference binds on the test point sets", {
  skip_if_not_installed("FNN")
  for (nm in names(.pr_sets())) {
    xy <- .pr_sets()[[nm]]
    d8 <- FNN::get.knn(xy, k = DEN_NN)$nn.dist[, DEN_NN]
    # The reference's floor of 1 unit, and the package's of 1e-6 of the extent.
    extent <- sqrt(sum((apply(xy, 2L, max) - apply(xy, 2L, min))^2))
    expect_gt(min(d8), 1)
    expect_gt(min(d8) / extent, spatialkit:::.DENSITY_FLOOR_REL)
    # With no floor in play the reference's density is k / (pi d^2) exactly.
    expect_identical(local_density(xy), DEN_NN / (pi * d8^2))
  }
})

test_that("the density kind's centres are den_centres()'s, to 1e-8", {
  skip_if_not_installed("FNN")
  worst <- 0
  for (nm in names(.pr_sets())) {
    xy <- .pr_sets()[[nm]]
    xs <- partition_sorted(xy)
    for (k in c(2L, 7L, 30L, 120L)) for (seed in c(1L, 77L)) {
      ref <- den_centres(xs, k, seed)
      got <- .pr_density(xy, k, seed)
      expect_identical(dim(got$centers), dim(ref))
      d <- max(abs(got$centers - ref))
      worst <- max(worst, d)
      expect_lt(d, 1e-8)
      # The cells are the reference's as well: every point to its nearest
      # centre.
      expect_identical(got$cluster, as.integer(partition_nearest(xy, ref)))
    }
  }
  # Measured: 2e-12 at most on these coordinates, which run to 2000.
  expect_lt(worst, 1e-10)
})

test_that("the weights are the reference's: local density to the power 1, mean 1", {
  skip_if_not_installed("FNN")
  xy <- .pr_sets()$clustered
  w_ref <- local_density(xy)^DEN_POWER
  w_ref <- w_ref / mean(w_ref)
  got <- .pr_density(xy, 12L, 3L)
  expect_equal(got$weights, w_ref, tolerance = 1e-12)
  expect_equal(mean(got$weights), 1, tolerance = 1e-12)
  # The objective is the reference's weighted sum at the centres returned.
  cl <- partition_nearest(xy, got$centers)
  expect_equal(got$objective, sum(w_ref * rowSums((xy - got$centers[cl, ])^2)),
               tolerance = 1e-10)
})

test_that("at projected magnitudes the two agree to the reference's own rounding", {
  # The same comparison at coordinates of 5e5 and 5e6, where one unit in the
  # last place of a coordinate is 9.3e-10.  The reference adds up its
  # weighted means at that magnitude; the package centres the coordinates on
  # their bounding box first and adds up numbers no larger than the extent.
  # The two then differ by a few units in the last place of the coordinates:
  # 2.8e-9 to 9.3e-9 on this layer.  The size depends on the layout and
  # grows with the points per cell, so other layers pass the 1e-8 that holds
  # at small coordinates.  It is the reference's rounding and not a
  # difference of method: the cells are identical, and against the weighted
  # mean of each cell computed on centred coordinates the package is within
  # 1.2e-10 where the reference is up to 9.3e-9 away.
  skip_if_not_installed("FNN")
  xy <- partition_clustered_xy(n_bg = 500L, n_cl = 2500L, clusters = 6L, side = 1e5,
                               origin = c(5e5, 5e6), seed = 3L)
  xs <- partition_sorted(xy)
  mid <- colMeans(xy)
  xc  <- sweep(xy, 2L, mid)
  for (k in c(7L, 30L, 120L)) {
    ref <- den_centres(xs, k, 77L)
    got <- .pr_density(xy, k, 77L)
    expect_lt(max(abs(got$centers - ref)), 1e-6)
    expect_identical(got$cluster, as.integer(partition_nearest(xy, ref)))
    wm <- sweep(rowsum(xc * got$weights, got$cluster) /
                  as.numeric(rowsum(got$weights, got$cluster)), 2L, mid, "+")
    expect_lt(max(abs(got$centers - wm)), 1e-9)
    expect_lte(max(abs(got$centers - wm)), max(abs(ref - wm)))
  }
})

test_that("where the reference's floor of one unit binds, the package's does not", {
  # The clustered layout on a square of side 25: the 8th neighbour of most
  # points is nearer than 1, so the reference reads the same density off all
  # of them.  It is the same layout as on a square of side 1000, where the
  # floor does not bind, and the reference gives a different partition of it.
  # The package gives the same one at either size.
  skip_if_not_installed("FNN")
  big   <- partition_clustered_xy(n_bg = 100L, n_cl = 200L, clusters = 3L, side = 1000, seed = 8L)
  small <- big / 40
  d8 <- FNN::get.knn(small, k = DEN_NN)$nn.dist[, DEN_NN]
  expect_gt(mean(d8 < 1), 0.5)
  expect_gt(min(FNN::get.knn(big, k = DEN_NN)$nn.dist[, DEN_NN]), 1)
  k <- 12L
  ref_big   <- den_centres(partition_sorted(big), k, 4L)
  ref_small <- den_centres(partition_sorted(small), k, 4L)
  expect_gt(max(abs(ref_small * 40 - ref_big)), 1)          # the unit changed the answer
  got_big   <- .pr_density(big, k, 4L)
  got_small <- .pr_density(small, k, 4L)
  expect_lt(max(abs(got_small$centers * 40 - got_big$centers)), 1e-8)
  expect_identical(got_small$cluster, got_big$cluster)
  # And the package's answer is the reference's where the reference is
  # unit-free.
  expect_lt(max(abs(got_big$centers - ref_big)), 1e-8)
})

test_that("get_voronoi_seeds(partition = \"density\") returns den_centres()'s centres", {
  # Through the public function: 3 starts, 40 Lloyd steps, 8 neighbours,
  # power 1, the settings of the experiment.
  skip_if_not_installed("FNN")
  xy  <- .pr_sets()$clustered
  pts <- partition_points_sf(xy)
  for (k in c(9L, 40L)) {
    ref <- den_centres(partition_sorted(xy), k, 21L)
    s <- get_voronoi_seeds(method = "kmeans", n = k, sample_points = pts,
                           partition = "density", kmeans_nstart = DEN_NSTART,
                           kmeans_iter = DEN_ITER, set_seed = 21L)
    expect_lt(max(abs(partition_xy(s) - ref)), 1e-8)
    expect_identical(attr(s, "partition")$kind, "density")
    expect_identical(attr(s, "partition")$init, "kmeans++")
  }
})

test_that("a density profile's centres are the engine's at the profile's settings", {
  # resolution_profile() runs `nstart` starts of at most 100 Lloyd steps.
  # Geometry only, so nothing draws from the stream before the sweep.
  xy  <- .pr_sets()$uniform
  pts <- partition_points_sf(xy)
  prof <- resolution_profile(pts, levels = 15L, nstart = 4, seed = 9, partition = "density")
  set.seed(9)
  eng <- spatialkit:::.partition_centres(xy, 15L, partition = "density", nstart = 4L,
                                         iter.max = 100L)
  expect_identical(attr(prof, "centres")[["15"]], eng$centers)
  expect_identical(prof$wss, eng$wss)
  expect_identical(unname(attr(prof, "partition")$objective), eng$objective)
})

test_that("with the reference's 40 steps a density profile's centres are den_centres()'s", {
  # Every start of this case stops moving before step 40, so the limit of
  # 100 the profile runs to makes no difference.
  skip_if_not_installed("FNN")
  xy  <- .pr_sets()$clustered
  pts <- partition_points_sf(xy)
  k <- 30L
  set.seed(1)
  eng <- spatialkit:::.partition_centres(xy, k, partition = "density", nstart = DEN_NSTART,
                                         iter.max = DEN_ITER)
  expect_true(eng$converged)
  prof <- resolution_profile(pts, levels = k, nstart = DEN_NSTART, seed = 1,
                             partition = "density")
  expect_lt(max(abs(attr(prof, "centres")[[as.character(k)]] -
                      den_centres(partition_sorted(xy), k, 1L))), 1e-8)
})


# ---- Random starts ----------------------------------------------------------

test_that("random starts are rnd_centres()'s: stats::kmeans() from 5 starts, 100 iterations", {
  xy <- .pr_sets()$clustered
  for (k in c(3L, 20L, 150L)) for (seed in c(4L, 9L)) {
    ref <- rnd_centres(xy, k, seed)
    expect_identical(attr(ref, "warned"), 0L)
    expect_false(attr(ref, "fell"))
    ref <- matrix(as.numeric(ref), ncol = 2L)
    # One start at a time, as the profile runs them ...
    set.seed(seed)
    each <- spatialkit:::.partition_centres(xy, k, partition = "kmeans", init = "random",
                                            nstart = RND_NSTART, iter.max = RND_ITER)
    # ... and as one stats::kmeans() call, as the seeding functions run them.
    set.seed(seed)
    one <- spatialkit:::.partition_centres(xy, k, partition = "kmeans", init = "random",
                                           nstart = RND_NSTART, iter.max = RND_ITER,
                                           one_call = TRUE)
    expect_identical(each$centers, ref)
    expect_identical(one$centers, ref)
    expect_identical(each$cluster, one$cluster)
    expect_identical(each$wss, one$wss)
    # The start-by-start run also has the spread of the starts.
    expect_true(is.finite(each$spread) && each$spread >= 0)
    expect_true(is.na(one$spread))
  }
})

test_that("a single random start is the one stats::kmeans(nstart = 1) draws", {
  # With one start stats::kmeans() draws among all the rows and turns to the
  # distinct rows only when two of those drawn coincide.
  xy  <- .pr_sets()$uniform
  dup <- rbind(xy[1:60, ], xy[1:60, ], xy[1:30, ])
  for (x in list(xy, dup)) for (seed in 1:3) {
    set.seed(seed)
    km <- suppressWarnings(stats::kmeans(x, centers = 12L, nstart = 1L, iter.max = 100L))
    set.seed(seed)
    got <- spatialkit:::.partition_centres(x, 12L, partition = "kmeans", init = "random",
                                           nstart = 1L, iter.max = 100L)
    expect_identical(got$centers, unname(km$centers))
  }
})

test_that("get_voronoi_seeds() and a random-start profile return rnd_centres()'s centres", {
  xy  <- .pr_sets()$clustered
  pts <- partition_points_sf(xy)
  k <- 20L
  # The seeding function clusters the rows as they come, as rnd_centres() does.
  ref <- matrix(as.numeric(rnd_centres(xy, k, 4L)), ncol = 2L)
  s <- get_voronoi_seeds(method = "kmeans", n = k, sample_points = pts,
                         kmeans_nstart = RND_NSTART, kmeans_iter = RND_ITER, set_seed = 4L)
  expect_identical(partition_xy(s), ref)
  # The profile puts the rows in coordinate order first.
  ref_sorted <- matrix(as.numeric(rnd_centres(partition_sorted(xy), k, 4L)), ncol = 2L)
  prof <- resolution_profile(pts, levels = k, nstart = RND_NSTART, seed = 4, init = "random")
  expect_identical(attr(prof, "centres")[[as.character(k)]], ref_sorted)
  expect_identical(attr(prof, "partition")[c("kind", "init")],
                   list(kind = "kmeans", init = "random"))
})
