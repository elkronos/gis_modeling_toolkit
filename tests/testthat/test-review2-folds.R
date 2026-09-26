# tests/testthat/test-review2-folds.R
# ---------------------------------------------------------------------------
# Second review round: make_folds().  Each test fails on the code before
# the fix.  Helpers are in helper-review2-folds.R.
# ---------------------------------------------------------------------------

# ---- block_kfold: k and the blocks that hold points ------------------------

test_that("drop_empty_blocks = FALSE lowers k to the blocks that hold points", {
  # Two clusters on a 4 x 4 grid: 2 occupied blocks, but the highest occupied
  # block id is 16, so k = 5 used to be kept and three folds came back with
  # no test points and no warning.
  set.seed(1)
  xy <- rbind(cbind(rnorm(30, 100, 20), rnorm(30, 100, 20)),
              cbind(rnorm(30, 900, 20), rnorm(30, 900, 20)))
  p <- r2_pts(xy[, 1], xy[, 2])
  lines <- capture_spatialkit_log(
    f <- make_folds(p, k = 5, method = "block_kfold", block_nx = 4, block_ny = 4,
                    drop_empty_blocks = FALSE, seed = 1))
  expect_equal(f$k, 2L)
  expect_length(f$folds, 2L)
  expect_true(all(lengths(lapply(f$folds, `[[`, "test")) > 0L))
  expect_true(is.finite(f$params$balance_ratio))
  expect_true(log_has(lines, "only 2 blocks hold points"))
  # The empty blocks are still kept and packed.
  expect_equal(f$params$blocks_used, 16L)
  expect_identical(sort(unlist(f$params$fold_blocks)), 1:16)

  # All points in one of several blocks is the single-block case, not a fold
  # scheme with an empty training set.
  bnd <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(4000, 0), c(4000, 4000),
                                              c(0, 4000), c(0, 0)))), crs = 32632)
  one <- r2_pts(runif(20, 100, 300), runif(20, 100, 300))
  expect_error(make_folds(one, k = 3, method = "block_kfold", boundary = bnd,
                          block_nx = 4, block_ny = 4, drop_empty_blocks = FALSE),
               "single block")
})


# ---- block_kfold: the automatic grid ---------------------------------------

test_that("the automatic grid is the same for a corridor and for it turned on its side", {
  set.seed(1)
  u <- runif(150, 0, 10000); v <- runif(150, 0, 100)
  ew <- make_folds(r2_pts(5e5 + u, 5e6 + v), k = 5, method = "block_kfold", seed = 1)
  ns <- make_folds(r2_pts(5e5 + v, 5e6 + u), k = 5, method = "block_kfold", seed = 1)
  # block_multiplier * k = 15 blocks either way; 39 x 1 before.
  expect_equal(c(ew$params$grid_nx, ew$params$grid_ny), c(15, 1))
  expect_equal(c(ns$params$grid_nx, ns$params$grid_ny), c(1, 15))
  expect_equal(ew$k, ns$k)
  # A square keeps its aspect-preserving grid.
  sq <- make_folds(r2_pts(runif(150, 0, 1000), runif(150, 0, 1000)), k = 5,
                   method = "block_kfold", seed = 1)
  expect_equal(c(sq$params$grid_nx, sq$params$grid_ny), c(4, 4))
})

test_that("points on one horizontal line get a row of blocks, not a collapsed square grid", {
  set.seed(2)
  h <- r2_pts(5e5 + runif(150, 0, 6000), rep(5e6, 150))
  f <- make_folds(h, k = 5, method = "block_kfold", seed = 1)
  expect_equal(c(f$params$grid_nx, f$params$grid_ny), c(15, 1))
  expect_equal(f$k, 5L)                        # was lowered to 4
})


# ---- block_kfold: block_nx / block_ny --------------------------------------

test_that("one grid dimension is honoured and the other derived; bad ones are refused", {
  set.seed(3)
  d <- r2_pts(runif(200, 0, 1000), runif(200, 0, 2000))
  f <- make_folds(d, k = 4, method = "block_kfold", block_nx = 10, seed = 1)
  expect_equal(f$params$grid_nx, 10)           # was ignored: a 3 x 4 grid
  expect_equal(f$params$grid_ny, 20)
  f <- make_folds(d, k = 4, method = "block_kfold", block_ny = 4, seed = 1)
  expect_equal(c(f$params$grid_nx, f$params$grid_ny), c(2, 4))
  for (v in list(0, -2, NA, c(2, 3), 2.7, "3"))
    expect_error(make_folds(d, k = 4, method = "block_kfold", block_nx = v,
                            block_ny = 3),
                 "`block_nx` must be a single whole number >= 1")
  expect_error(make_folds(d, k = 4, method = "block_kfold", block_nx = 3,
                          block_ny = NA),
               "`block_ny` must be a single whole number >= 1")
})

test_that("a units object is refused by name for block_size", {
  d <- r2_pts(runif(50, 0, 1000), runif(50, 0, 1000))
  expect_error(make_folds(d, k = 3, method = "block_kfold",
                          block_size = units::set_units(1, km)),
               "`block_size` must be a single positive number")
})


# ---- block_kfold: boundary clipping ----------------------------------------

test_that("with a boundary, source_row indexes the full grid and zero-area pieces are not blocks", {
  tri <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1000, 0), c(0, 1000),
                                              c(0, 0)))), crs = 32632)
  set.seed(1)
  pts <- sf::st_sf(geometry = sf::st_sample(tri, 400))
  f <- make_folds(pts, k = 3, method = "block_kfold", block_nx = 4, block_ny = 4,
                  boundary = tri, seed = 1, drop_empty_blocks = FALSE)
  pr <- f$params
  expect_equal(pr$n_blocks, 16L)
  expect_true(all(as.numeric(sf::st_area(pr$blocks)) > 0))
  expect_equal(pr$blocks_used, 10L)
  # Cells 8, 11 and 12 touch the hypotenuse only at a corner, and 14-16 lie
  # outside it; the clipped list numbered them 1..13 instead.
  expect_identical(pr$blocks$source_row, c(1:7, 9L, 10L, 13L))
  full <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(tri)), n = c(4, 4))
  inside <- sf::st_within(sf::st_centroid(sf::st_geometry(pr$blocks)), full)
  expect_identical(vapply(inside, `[`, integer(1), 1L), pr$blocks$source_row)

  # A point exactly on a grid vertex that lies on the boundary joins an
  # areal block, not a one-point POINT block of its own.
  tri2 <- sf::st_sfc(sf::st_polygon(list(rbind(c(1000, 0), c(1000, 1000), c(0, 1000),
                                               c(1000, 0)))), crs = 32632)
  p2 <- sf::st_sf(geometry = c(sf::st_sample(tri2, 200),
                               sf::st_sfc(sf::st_point(c(500, 500)), crs = 32632)))
  f2 <- make_folds(p2, k = 3, method = "block_kfold", block_nx = 4, block_ny = 4,
                   boundary = tri2, seed = 1)
  blk <- f2$params$blocks[f2$assignment$block_id[nrow(p2)], ]
  expect_gt(as.numeric(sf::st_area(blk)), 0)
  expect_true(all(sf::st_dimension(f2$params$blocks) == 2L))
})


# ---- argument validation ---------------------------------------------------

test_that("k is required by the k-fold methods and named when missing", {
  d <- r2_pts(runif(40, 0, 1000), runif(40, 0, 1000), g = rep(1:4, 10))
  expect_error(make_folds(d, method = "block_kfold"),
               "`k` \\(the number of folds\\) is required for method = \"block_kfold\"")
  expect_error(make_folds(d, k = NULL, method = "random_kfold"),
               "`k` .* is required for method = \"random_kfold\"")
  expect_error(make_folds(d, k = NULL, method = "leave_location_out", group_var = "g"),
               "`k` .* is required")
  # The leave-one-out methods never read it.
  expect_equal(make_folds(d, method = "buffered_loo", buffer = 50)$k, 40L)
})

test_that("an invalid buffer is refused by name", {
  d <- r2_pts(runif(40, 0, 1000), runif(40, 0, 1000))
  for (b in list(NA_real_, numeric(0), c(100, 200), NULL, "100", -1,
                 units::set_units(1, km)))
    expect_error(make_folds(d, k = 1, method = "buffered_loo", buffer = b),
                 "`buffer` must be a single positive number")
  # An NA -- what estimate_sac_range() returns when nothing is identified --
  # says so.
  expect_error(make_folds(d, k = 1, method = "buffered_loo", buffer = NA_real_),
               "no range was identified")
})

test_that("a buffer that excludes no neighbour is warned about", {
  # 0.1 on lon/lat input is 0.1 m once projected: plain LOO.
  set.seed(1)
  ll <- r2_pts(10 + runif(80, 0, 1), 50 + runif(80, 0, 1), crs = 4326)
  expect_warning(
    f <- r2_quiet(make_folds(ll, k = 1, method = "buffered_loo", buffer = 0.1)),
    "excludes no neighbour from any fold, so this is plain leave-one-out")
  expect_true(all(lengths(lapply(f$folds, `[[`, "train")) == 79L))
  # A buffer that does exclude something is not.
  expect_no_warning(r2_quiet(make_folds(ll, k = 1, method = "buffered_loo",
                                        buffer = 20000)))
})


# ---- auto_range fallback ---------------------------------------------------

test_that("auto_range falling back to geometric blocks is a warning naming the reason", {
  skip_if_not_installed("gstat")
  set.seed(4); n <- 120
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- r2_pts(x, y, z = 0.01 * x + rnorm(n, 0, 0.1))     # an unremoved trend
  expect_warning(
    f <- r2_quiet(make_folds(d, k = 4, method = "block_kfold", auto_range = TRUE,
                             range_frac = 1e-6, response_var = "z", seed = 1)),
    "no autocorrelation range was identified \\(fitted range exceeds the largest lag fitted\\); falling back to geometric blocks")
  expect_null(f$params$block_size)
})


# ---- nndm ------------------------------------------------------------------

test_that("nndm warns when min_train leaves the folds more optimistic than the target", {
  # One cluster predicted onto a 20 km grid: 96 of 100 folds are held at the
  # floor with a training point far closer than the prediction distances.
  set.seed(1)
  p <- r2_pts(rnorm(100, 10000, 800), rnorm(100, 10000, 800))
  g <- r2_pts(rep(seq(0, 20000, by = 1000), 21), rep(seq(0, 20000, by = 1000), each = 21))
  expect_warning(
    f <- r2_quiet(make_folds(p, method = "nndm", prediction_points = g)),
    "min_train = 0.5 stopped the distance matching in 96 of 100 folds")
  expect_equal(f$params$n_at_min_train, 96L)
  expect_gt(f$params$max_ecdf_excess, 0.5)
  # Limiting the matching with phi is a choice, not the floor: no warning.
  expect_no_warning(
    f_phi <- r2_quiet(make_folds(p, method = "nndm", prediction_points = g, phi = 500)))
  expect_equal(f_phi$params$n_at_min_train, 0L)
  # Nor where the target is reachable.
  set.seed(2)
  u <- r2_pts(runif(100, 0, 10000), runif(100, 0, 10000))
  g2 <- r2_pts(rep(seq(0, 10000, by = 500), 21), rep(seq(0, 10000, by = 500), each = 21))
  expect_no_warning(fu <- r2_quiet(make_folds(u, method = "nndm", prediction_points = g2)))
  expect_equal(fu$params$n_at_min_train, 0L)
})

test_that("the nndm size guard states the worst-case cost", {
  set.seed(1)
  big <- r2_pts(runif(5001, 0, 1e4), runif(5001, 0, 1e4))
  expect_error(make_folds(big, method = "nndm", prediction_points = big[1:10, ]),
               "worst case is O\\(n\\^3\\) time")
})
