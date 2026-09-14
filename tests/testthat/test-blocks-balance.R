# tests/testthat/test-blocks-balance.R
# ---------------------------------------------------------------------------
# make_folds(block_kfold): user-supplied blocks and the fold-balance tolerance.
#
# The grid used to be the only block design the fold builder could consume,
# while build_tessellation() produced every shape it could not.  `blocks`
# closes that gap: a polygon layer supplies the blocks and the points are
# assigned to folds exactly as grid cells are.  `balance_tol` exposes the
# 3:1 residual-imbalance check, which is now a warning a pipeline can catch
# rather than a log line only.
# ---------------------------------------------------------------------------

blk_points <- function(n = 400, seed = 1, extent = 1000) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = runif(n, 0, extent), y = runif(n, 0, extent),
                          z = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}

# Six clusters of unequal weight: the geometric grid packs them unevenly.
blk_clustered <- function(n = 300, seed = 3, extent = 1000) {
  set.seed(seed)
  cx <- runif(6, 0, extent); cy <- runif(6, 0, extent)
  pr <- rexp(6); pr <- pr / sum(pr)
  i <- sample(6, n, replace = TRUE, prob = pr)
  x <- pmin(pmax(cx[i] + rnorm(n, 0, extent / 15), 0), extent)
  y <- pmin(pmax(cy[i] + rnorm(n, 0, extent / 15), 0), extent)
  sf::st_as_sf(data.frame(x = x, y = y, z = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}

# The leakage fixture of test-leakage-diagnostic.R: an exponential field
# with range parameter 70 (effective range ~210, estimated ~200) on a
# 1000-unit extent.
blk_leak_points <- function(n = 300, seed = 8) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 70) + diag(1e-8, n))) %*% rnorm(n)) +
    rnorm(n, sd = 0.3)
  sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}

square_bnd <- function(extent = 1000, crs = 32632) {
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(extent, 0), c(extent, extent), c(0, extent), c(0, 0)))),
    crs = crs))
}

rect_poly <- function(x0, x1, y0, y1) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
}

# Which supplied block each point falls in, computed independently of
# make_folds(): first hit, as the function documents.
block_of <- function(pts, blocks) {
  h <- sf::st_intersects(pts, blocks)
  vapply(h, function(ix) if (length(ix)) ix[1] else NA_integer_, 1L)
}

hex_blocks <- function(pts, n_cells = 30) {
  build_tessellation(pts, boundary = square_bnd(), method = "hex",
                     approx_n_cells = n_cells, quiet = TRUE)$cells
}

catch_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}


# ---- 7.1: supplied blocks ---------------------------------------------------

test_that("supplied polygons are the blocks: every block's points share a fold", {
  pts <- blk_points()
  hex <- hex_blocks(pts)
  expect_gt(nrow(hex), 20L)

  f <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, blocks = hex)
  expect_equal(f$k, 5L)
  expect_length(f$folds, 5L)
  expect_equal(sort(f$assignment$row_id), seq_len(nrow(pts)))
  expect_true(all(f$assignment$fold %in% 1:5))

  bid <- block_of(pts, hex)
  expect_false(anyNA(bid))
  folds_per_block <- tapply(f$assignment$fold, bid, function(v) length(unique(v)))
  expect_true(all(folds_per_block == 1L))

  # Each fold's train and test partition the rows.
  for (fd in f$folds) {
    expect_length(intersect(fd$train, fd$test), 0L)
    expect_equal(sort(c(fd$train, fd$test)), seq_len(nrow(pts)))
  }
})

test_that("params record the supplied design and no grid", {
  pts <- blk_points()
  hex <- hex_blocks(pts)
  f <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, blocks = hex)
  p <- f$params
  expect_true(isTRUE(p$blocks_supplied))
  expect_true(is.na(p$grid_nx)); expect_true(is.na(p$grid_ny))
  expect_equal(p$n_blocks, nrow(hex))
  expect_equal(p$blocks_used, length(unique(block_of(pts, hex))))
  expect_lte(p$blocks_used, p$n_blocks)
  expect_null(p$block_size)
  # Median side of the equal-area square over blocks holding points: for
  # ~30 hexagons on a 1000 x 1000 extent, about 180 units.
  expect_true(is.finite(p$block_scale))
  expect_gt(p$block_scale, 100); expect_lt(p$block_scale, 300)
  expect_equal(p$balance_ratio,
               max(table(f$assignment$fold)) / min(table(f$assignment$fold)))
  expect_equal(p$balance_tol, 3)
  expect_equal(p$crs, "EPSG:32632")

  # A grid run reports the counterpart fields.
  g <- make_folds(pts, k = 5, method = "block_kfold", seed = 1)
  expect_false(g$params$blocks_supplied)
  expect_true(is.na(g$params$block_scale))
  expect_gte(g$params$n_blocks, g$params$blocks_used)
  expect_equal(g$params$blocks_used, g$params$grid_nx * g$params$grid_ny)
})

test_that("an sfc layer and its sf wrapper give the same folds", {
  pts <- blk_points()
  hex <- hex_blocks(pts)
  f1 <- make_folds(pts, k = 4, method = "block_kfold", seed = 2, blocks = hex)
  f2 <- make_folds(pts, k = 4, method = "block_kfold", seed = 2,
                   blocks = sf::st_geometry(hex))
  expect_identical(f1$assignment, f2$assignment)
})

test_that("grid-sizing arguments and boundary are ignored, and say so", {
  pts <- blk_points()
  hex <- hex_blocks(pts)
  ref <- make_folds(pts, k = 4, method = "block_kfold", seed = 2, blocks = hex)
  lines <- capture_spatialkit_log(
    f <- make_folds(pts, k = 4, method = "block_kfold", seed = 2, blocks = hex,
                    block_size = 50, block_nx = 20, block_ny = 20,
                    boundary = square_bnd()))
  expect_identical(f$assignment, ref$assignment)
  expect_null(f$params$block_size)
  expect_false(f$params$boundary_supplied)
  expect_true(log_has(lines, "block_size.*ignored when `blocks` are supplied"))
  expect_true(log_has(lines, "boundary.*ignored when `blocks` are supplied"))
})

test_that("blocks are refused for every other method, and non-polygon input is refused", {
  pts <- blk_points(n = 60)
  hex <- hex_blocks(pts)
  for (m in c("random_kfold", "buffered_loo", "leave_location_out", "nndm"))
    expect_error(make_folds(pts, k = 3, method = m, blocks = hex),
                 "only used by method = \"block_kfold\"")
  tess <- build_tessellation(pts, boundary = square_bnd(), method = "hex",
                             approx_n_cells = 20, quiet = TRUE)
  expect_error(make_folds(pts, k = 3, method = "block_kfold", blocks = tess),
               "pass its `\\$cells`")
  expect_error(make_folds(pts, k = 3, method = "block_kfold", blocks = pts),
               "POLYGON, MULTIPOLYGON")
  expect_error(make_folds(pts, k = 3, method = "block_kfold", blocks = square_bnd()),
               "at least 2 polygons")
})

test_that("points outside every block go to the nearest one, with a warning that counts them", {
  pts <- blk_points()
  hex <- hex_blocks(pts)
  cx <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(hex))))[, 1]
  west <- hex[cx < 600, ]
  bid_west <- block_of(pts, west)
  n_out <- sum(is.na(bid_west))
  expect_gt(n_out, 50L)

  res <- catch_warnings(make_folds(pts, k = 4, method = "block_kfold", seed = 1,
                                   blocks = west))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, sprintf("%d of %d points fall outside every supplied block", n_out, nrow(pts)))
  expect_match(res$warnings, "nearest block")
  f <- res$value
  expect_equal(sort(f$assignment$row_id), seq_len(nrow(pts)))
  expect_true(all(f$assignment$fold %in% seq_len(f$k)))

  # Nearest by distance to the polygon: the outside points and the block
  # nearest to each must share a fold.
  near <- sf::st_nearest_feature(pts[is.na(bid_west), ], west)
  fold_of_block <- tapply(f$assignment$fold[!is.na(bid_west)], bid_west[!is.na(bid_west)],
                          function(v) unique(v)[1])
  expect_equal(unname(f$assignment$fold[is.na(bid_west)]),
               as.vector(fold_of_block[as.character(near)]))
})

test_that("a rounding-error edge is not 'outside': a lon/lat round trip reproduces the folds silently", {
  pts <- blk_points()
  hex <- hex_blocks(pts)
  ref <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, blocks = hex)
  res <- catch_warnings(make_folds(sf::st_transform(pts, 4326), k = 5,
                                   method = "block_kfold", seed = 1,
                                   blocks = sf::st_transform(hex, 4326)))
  expect_length(res$warnings, 0L)
  expect_identical(res$value$assignment$fold, ref$assignment$fold)
})

test_that("overlapping blocks are warned about; blocks that only share edges are not", {
  pts <- blk_points()
  ov <- sf::st_sf(geometry = sf::st_sfc(rect_poly(0, 600, 0, 1000),
                                        rect_poly(400, 1000, 0, 1000)),
                  crs = 32632)
  n_both <- sum(lengths(sf::st_intersects(pts, ov)) > 1L)
  expect_gt(n_both, 0L)
  res <- catch_warnings(make_folds(pts, k = 2, method = "block_kfold", seed = 1,
                                   blocks = ov))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, sprintf("%d point\\(s\\) fall inside more than one", n_both))
  expect_match(res$warnings, "overlap")
  # First block wins: every point of the overlap sits in block 1's fold.
  f <- res$value
  in_first <- lengths(sf::st_intersects(pts, ov[1, ])) > 0L
  expect_length(unique(f$assignment$fold[in_first]), 1L)

  # A grid of touching squares: points on shared edges hit two cells, and
  # the design is a partition, so nothing is said.
  grid <- sf::st_sf(geometry = sf::st_make_grid(square_bnd(), n = c(4, 4)), crs = 32632)
  on_edge <- sf::st_as_sf(data.frame(x = c(250, 500, 750, 100), y = c(100, 500, 300, 750),
                                     z = 0),
                          coords = c("x", "y"), crs = 32632)
  pts2 <- rbind(pts, on_edge)
  expect_gt(sum(lengths(sf::st_intersects(pts2, grid)) > 1L), 0L)
  res2 <- catch_warnings(make_folds(pts2, k = 4, method = "block_kfold", seed = 1,
                                    blocks = grid))
  expect_length(res2$warnings, 0L)
})

test_that("a design that gives every point the same block, or none, is an error", {
  pts <- blk_points(n = 80)
  one <- sf::st_sf(geometry = sf::st_sfc(rect_poly(-1, 2000, -1, 2000),
                                         rect_poly(5000, 6000, 5000, 6000)),
                   crs = 32632)
  expect_error(make_folds(pts, k = 3, method = "block_kfold", blocks = one),
               "all 80 points fall in the same one of the supplied `blocks`")
  far <- sf::st_sf(geometry = sf::st_sfc(rect_poly(5000, 6000, 5000, 6000),
                                         rect_poly(7000, 8000, 5000, 6000)),
                   crs = 32632)
  expect_error(make_folds(pts, k = 3, method = "block_kfold", blocks = far),
               "none of the 80 points fall inside any of the supplied `blocks`")
})

test_that("k is lowered to the number of blocks that hold points", {
  pts <- blk_points(n = 80)
  three <- sf::st_sf(geometry = sf::st_sfc(rect_poly(0, 400, 0, 1000),
                                           rect_poly(400, 700, 0, 1000),
                                           rect_poly(700, 1000, 0, 1000),
                                           rect_poly(3000, 4000, 0, 1000)),
                     crs = 32632)
  lines <- capture_spatialkit_log(
    f <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, blocks = three))
  expect_equal(f$k, 3L)
  expect_equal(f$params$blocks_used, 3L)
  expect_equal(f$params$n_blocks, 4L)
  expect_true(log_has(lines, "blocks < k; reducing k"))
  # Without dropping empties the fourth (empty) block is still absent from
  # every fold, and the folds are the same.
  f2 <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, blocks = three,
                   drop_empty_blocks = FALSE)
  expect_identical(f2$assignment, f$assignment)
})

test_that("CRS-less points are aligned to blocks that carry one", {
  pts <- blk_points()
  hex <- hex_blocks(pts)
  ref <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, blocks = hex)
  p0 <- sf::st_set_crs(pts, NA)
  f <- suppressWarnings(make_folds(p0, k = 5, method = "block_kfold", seed = 1,
                                   blocks = hex))
  expect_identical(f$assignment$fold, ref$assignment$fold)
  expect_equal(f$params$crs, "EPSG:32632")
  # Both CRS-less: the same unnamed space, no reprojection.
  f0 <- suppressWarnings(make_folds(p0, k = 5, method = "block_kfold", seed = 1,
                                    blocks = sf::st_set_crs(hex, NA)))
  expect_identical(f0$assignment$fold, ref$assignment$fold)
  expect_true(is.na(f0$params$crs))
})

test_that("auto_range with supplied blocks compares the range against them and resizes nothing", {
  skip_if_not_installed("gstat")
  pts <- blk_leak_points()      # range ~200 on a 1000-unit extent
  # Small hexagons (~100 units across): well below the range.
  small <- build_tessellation(pts, boundary = square_bnd(), method = "hex",
                              approx_n_cells = 100, quiet = TRUE)$cells
  expect_lt(median(sqrt(as.numeric(sf::st_area(small)))), 150)
  res <- catch_warnings(
    expect_message(
      f <- make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                      blocks = small, auto_range = TRUE, response_var = "z"),
      "used as given, so the range is only compared against them"))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "supplied blocks \\(median scale [0-9.]+\\) are smaller than the autocorrelation range")
  expect_null(f$params$block_size)
  expect_true(is.finite(f$params$sac_range))
  expect_lt(f$params$block_scale, f$params$sac_range)
  expect_equal(f$params$n_blocks, nrow(small))

  # Large blocks: the same call is silent.
  big <- build_tessellation(pts, boundary = square_bnd(), method = "hex",
                            approx_n_cells = 6, quiet = TRUE)$cells
  res2 <- catch_warnings(suppressMessages(
    f2 <- make_folds(pts, k = 3, method = "block_kfold", seed = 1,
                     blocks = big, auto_range = TRUE, response_var = "z")))
  expect_length(res2$warnings, 0L)
  expect_gt(f2$params$block_scale, f2$params$sac_range)
})

test_that("the leakage diagnostic also runs on the default path with supplied blocks", {
  skip_if_not_installed("gstat")
  pts <- blk_leak_points()
  small <- build_tessellation(pts, boundary = square_bnd(), method = "hex",
                              approx_n_cells = 100, quiet = TRUE)$cells
  res <- catch_warnings(make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                                   blocks = small, response_var = "z"))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "smaller than the autocorrelation range")
})


# ---- 7.2: fold balance ------------------------------------------------------

test_that("balance_tol is validated", {
  pts <- blk_points(n = 60)
  expect_error(make_folds(pts, k = 3, method = "block_kfold", balance_tol = 0.5),
               "`balance_tol` must be a single number >= 1")
  expect_error(make_folds(pts, k = 3, method = "block_kfold", balance_tol = c(2, 3)),
               "`balance_tol`")
  expect_error(make_folds(pts, k = 3, method = "block_kfold", balance_tol = NA),
               "`balance_tol`")
  expect_error(make_folds(pts, k = 3, method = "block_kfold", balance_tol = "3"),
               "`balance_tol`")
})

test_that("the imbalance warning is an R condition, fires past the tolerance and not below it", {
  pts <- blk_clustered()
  ref <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, balance_tol = Inf)
  ratio <- ref$params$balance_ratio
  expect_equal(ratio, max(table(ref$assignment$fold)) / min(table(ref$assignment$fold)))
  expect_gt(ratio, 1.2)

  # Tolerance just below the achieved ratio: warns, names both counts, logs.
  lines <- capture_spatialkit_log(
    res <- catch_warnings(make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                                     balance_tol = ratio - 0.05)))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "fold size imbalance")
  expect_match(res$warnings, sprintf("largest fold has %d obs vs %d in smallest",
                                     max(table(ref$assignment$fold)),
                                     min(table(ref$assignment$fold))))
  expect_match(res$warnings, "points per block")
  expect_true(log_has(lines, "fold size imbalance"))
  expect_equal(res$value$params$balance_tol, ratio - 0.05)
  # The tolerance changes what is reported, never the folds.
  expect_identical(res$value$assignment, ref$assignment)

  # Tolerance just above: silent.
  res2 <- catch_warnings(make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                                    balance_tol = ratio + 0.05))
  expect_length(res2$warnings, 0L)
  expect_identical(res2$value$assignment, ref$assignment)
})

test_that("the default tolerance warns when one block holds most of the points", {
  # 85% of the points in one corner block: no packing can balance that.
  set.seed(9)
  n <- 200
  x <- c(runif(170, 0, 240), runif(30, 260, 1000))
  y <- c(runif(170, 0, 240), runif(30, 0, 1000))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = rnorm(n)),
                      coords = c("x", "y"), crs = 32632)
  res <- catch_warnings(make_folds(pts, k = 4, method = "block_kfold", seed = 1,
                                   block_size = 250))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "fold size imbalance")
  expect_match(res$warnings, "tolerance 3")
  expect_gt(res$value$params$balance_ratio, 3)
  # And Inf silences it without changing the folds.
  res2 <- catch_warnings(make_folds(pts, k = 4, method = "block_kfold", seed = 1,
                                    block_size = 250, balance_tol = Inf))
  expect_length(res2$warnings, 0L)
  expect_identical(res2$value$assignment, res$value$assignment)
})

test_that("density-adaptive blocks through `blocks` balance what the grid cannot", {
  pts <- blk_clustered()
  g <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, balance_tol = Inf)
  seeds <- suppressWarnings(get_voronoi_seeds(square_bnd(), method = "kmeans", n = 15,
                                              sample_points = pts, set_seed = 1))
  vor <- build_tessellation(seeds, boundary = square_bnd(), method = "voronoi",
                            quiet = TRUE)$cells
  v <- make_folds(pts, k = 5, method = "block_kfold", seed = 1, blocks = vor,
                  balance_tol = Inf)
  expect_lt(v$params$balance_ratio, g$params$balance_ratio)
  expect_lt(v$params$balance_ratio, 1.5)
})
