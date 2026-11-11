# tests/testthat/test-fix-edge-margin-offset.R
# ---------------------------------------------------------------------------
# BUG-225: the margin within which a point counts as on an edge
# (.edge_tol()) was one part in 1e10 of the largest coordinate, with no cap.
# On a layer 4 km across shifted to 1e12 that is 100 m, and at 1e15 100 km:
# block folds, the join to polygons, the Voronoi cells, cell_values_at() and
# the block-size sweep changed with the origin, and at 1e15 make_folds()
# refused with a false "single block".  The margin is now at most a
# millionth of the diagonal of the coordinates' box (.EDGE_REACH) and
# at least four units in the last place of the largest coordinate.
#
# The layers are built from fixed sequences, not the random stream, and
# every coordinate is a multiple of 1/8, so the shifts by 1e8 and 1e12 are
# exact: the layers are the same geometry in every frame, and the results
# are compared exactly.
# ---------------------------------------------------------------------------

.em_x0 <- 600000
.em_y0 <- 3990000

# 60 points over 4000 m x 3000 m (a Kronecker sequence), none on an edge of
# the polygons below, none repeated, with a response and two predictors.
.em_points <- function(by = 0) {
  i <- seq_len(60)
  px <- .em_x0 + round(8 * 4000 * ((i * 0.6180339887498949) %% 1)) / 8
  py <- .em_y0 + round(8 * 3000 * ((i * 0.7548776662466927) %% 1)) / 8
  on_x <- ((px - (.em_x0 - 250)) %% 1500) == 0
  px[on_x] <- px[on_x] + 0.125
  on_y <- ((py - (.em_y0 - 500)) %% 1000) == 0
  py[on_y] <- py[on_y] + 0.125
  u <- (px - .em_x0) / 1000
  v <- (py - .em_y0) / 1000
  x1 <- round(8192 * (0.8 * u + cos(3 * i))) / 8192
  x2 <- round(8192 * sin(5 * i)) / 8192
  y <- round(8192 * (2 + 1.5 * sin(1.2 * u) + cos(1.4 * v) + 0.6 * x1 - 0.4 * x2)) / 8192
  sf::st_as_sf(data.frame(pid = i, y = y, x1 = x1, x2 = x2, X = px + by, Y = py + by),
               coords = c("X", "Y"), crs = 32617)
}

# 8 x 8 points at 250 m: points on the polygons' edges, and tied distances.
.em_lattice <- function(by = 0) {
  ij <- expand.grid(i = 0:7, j = 0:7)
  sf::st_as_sf(data.frame(pid = seq_len(nrow(ij)), y = sin(ij$i) + cos(ij$j),
                          X = .em_x0 + 250 * ij$i + by, Y = .em_y0 + 250 * ij$j + by),
               coords = c("X", "Y"), crs = 32617)
}

# 12 rectangles (3 x 4, 1500 m x 1000 m) tiling 4500 m x 4000 m.
.em_polygons <- function(by = 0) {
  cells <- list()
  for (j in 0:3) for (i in 0:2) {
    xa <- .em_x0 - 250 + 1500 * i + by
    ya <- .em_y0 - 500 + 1000 * j + by
    cells[[length(cells) + 1L]] <- sf::st_polygon(list(matrix(
      c(xa, ya, xa + 1500, ya, xa + 1500, ya + 1000, xa, ya + 1000, xa, ya),
      ncol = 2L, byrow = TRUE)))
  }
  sf::st_sf(poly_id = sprintf("P%02d", seq_along(cells)),
            geometry = sf::st_sfc(cells, crs = 32617))
}

# The value of `expr`, with messages muffled and warnings collected.
.em_run <- function(expr) {
  w <- list()
  value <- withCallingHandlers(suppressMessages(expr), warning = function(c) {
    w[[length(w) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

.em_imprecise <- function(r)
  sum(vapply(r$warnings, inherits, logical(1), "spatialkit_xy_imprecise"))

.em_offsets <- c(0, 1e8, 1e12)

test_that(".edge_tol() is unchanged on a layer near its CRS's origin, and capped by the extent", {
  bb <- sf::st_bbox(.em_polygons())
  m <- max(abs(as.numeric(bb)))
  # The diagonal is more than a ten-thousandth of the largest coordinate:
  # one part in 1e10 of it, to the bit, in each form a caller passes.
  expect_identical(.edge_tol(bb), 1e-10 * m)
  expect_identical(.edge_tol(c(bb, sf::st_bbox(.em_points()))), 1e-10 * m)
  xy <- sf::st_coordinates(.em_points())
  expect_identical(.edge_tol(.xy_box(xy[, 1], xy[, 2])), 1e-10 * max(abs(xy)))
  xy <- sf::st_coordinates(.em_lattice())
  expect_identical(.edge_tol(.xy_box(xy[, 1], xy[, 2])), 1e-10 * max(abs(xy)))
  # A layer 200 m across at that northing: a millionth of its diagonal.
  expect_equal(.edge_tol(c(xmin = 6e5, ymin = 4e6, xmax = 6e5 + 200, ymax = 4e6)), 2e-4,
               tolerance = 1e-12)
  # Numbers without names give no diagonal: one part in 1e10 of the largest.
  expect_identical(.edge_tol(unname(as.numeric(bb))), 1e-10 * m)
  expect_identical(.edge_tol(c(0, 0, 0, 0)), 0)
  expect_identical(.edge_tol(c(xmin = NA_real_, ymin = NA, xmax = NA, ymax = NA)), 0)
  # Shifted to 1e12 the layer keeps its 4500 m x 4000 m diagonal, and the
  # margin is a millionth of it (6 mm), not 100 m.
  d <- sqrt(4500^2 + 4000^2)
  expect_equal(.edge_tol(sf::st_bbox(.em_polygons(1e12))), 1e-6 * d, tolerance = 1e-12)
  # At 1e15 that is below four units in the last place there (0.5 m).
  expect_identical(.edge_tol(sf::st_bbox(.em_polygons(1e15))), 4 * 0.125)
})

test_that("block folds are the same at offsets 0, 1e8 and 1e12", {
  for (mk in list(points = .em_points, lattice = .em_lattice)) {
    ref <- .em_run(make_folds(mk(0), k = 4, method = "block_kfold", seed = 1))
    for (by in .em_offsets[-1]) {
      r <- .em_run(make_folds(mk(by), k = 4, method = "block_kfold", seed = 1))
      expect_identical(r$value$assignment$fold, ref$value$assignment$fold,
                       info = sprintf("offset %g", by))
      expect_identical(r$value$params$blocks_used, ref$value$params$blocks_used)
      expect_identical(.em_imprecise(r), 0L)
    }
    # A block size given in map units, and a grid given by block_nx.
    for (args in list(list(block_size = 700), list(block_nx = 5L, block_ny = 3L))) {
      ref <- .em_run(do.call(make_folds, c(list(mk(0), k = 3, method = "block_kfold",
                                                seed = 1), args)))
      for (by in .em_offsets[-1]) {
        r <- .em_run(do.call(make_folds, c(list(mk(by), k = 3, method = "block_kfold",
                                                seed = 1), args)))
        expect_identical(r$value$assignment$fold, ref$value$assignment$fold,
                         info = sprintf("offset %g, %s", by, names(args)[1]))
      }
    }
  }
})

test_that("assign_features_to_polygons() gives each point the same polygon at offsets 0, 1e8 and 1e12", {
  for (mk in list(points = .em_points, lattice = .em_lattice)) {
    ref <- .em_run(assign_features_to_polygons(mk(0), .em_polygons(0), keep_unassigned = TRUE))
    for (by in .em_offsets[-1]) {
      for (largest in c(TRUE, FALSE)) {
        r <- .em_run(assign_features_to_polygons(mk(by), .em_polygons(by),
                                                 keep_unassigned = TRUE, largest = largest))
        expect_identical(r$value$poly_id, ref$value$poly_id, info = sprintf("offset %g", by))
        expect_identical(attr(r$value, "ties")$which, attr(ref$value, "ties")$which)
        expect_identical(.em_imprecise(r), 0L)
      }
    }
  }
})

test_that("Voronoi cells hold one point each, the same, at offsets 0, 1e8 and 1e12", {
  # The cells' IDs are the order GEOS returns them in, which the origin can
  # change; each point's own cell is compared by its area.  Areas computed
  # at 1e12, where doubles are 1.2e-4 m apart, differ from those at the
  # origin by about that spacing times a perimeter of a few km, a few parts
  # in ten million of a cell of 0.2 km2; 1e-5 is well above it and far
  # below the area of any other cell.
  own_area <- function(v) {
    a <- as.numeric(sf::st_area(v$cells))
    a[match(v$index, v$cells$cell_id)]
  }
  # Point 60 moved 40 m east of point 1: each lies 20 m from the edge of
  # their two cells, within a margin of 100 m, which put both in one cell.
  pts <- function(by) {
    p <- .em_points(by)
    xy <- sf::st_coordinates(p)
    sf::st_geometry(p)[[60]] <- sf::st_point(xy[1, ] + c(40, 0))
    p
  }
  for (build in c("create_voronoi_polygons", "build_tessellation")) {
    call <- function(by) .em_run(if (build == "create_voronoi_polygons")
      create_voronoi_polygons(pts(by), quiet = TRUE)
      else build_tessellation(pts(by), method = "voronoi", quiet = TRUE))$value
    ref <- call(0)
    expect_identical(length(unique(ref$index)), 60L)
    expect_false(anyNA(own_area(ref)))
    for (by in .em_offsets[-1]) {
      v <- call(by)
      expect_identical(nrow(v$cells), 60L)
      expect_identical(length(unique(v$index)), 60L, info = sprintf("%s, offset %g", build, by))
      expect_equal(own_area(v), own_area(ref), tolerance = 1e-5,
                   info = sprintf("%s, offset %g", build, by))
    }
  }
})

test_that(".seeds_as_centres() reads the same seeds as the centres at offsets 0, 1e8 and 1e12", {
  # Six centres 1000 m apart.  Seeds 50 m off them are not the centres at
  # any offset: with the margin of one part in 1e10 of the largest
  # coordinate, 100 m at 1e12, they were.  Seeds a few units in the last
  # place off them are the centres at every offset.
  cen <- cbind(.em_x0 + c(0, 1000, 2000, 0, 1000, 2000),
               .em_y0 + c(0, 0, 0, 1000, 1000, 1000))
  for (by in .em_offsets) {
    c_by <- cen + by
    for (fnn in unique(c(FALSE, spatialkit:::.has_backend("FNN")))) {
      expect_null(spatialkit:::.seeds_as_centres(c_by + 50, c_by, use_fnn = fnn),
                  info = sprintf("offset %g", by))
      near <- c_by + 2 * spatialkit:::.ulp_of(as.numeric(c_by))
      expect_identical(spatialkit:::.seeds_as_centres(near, c_by, use_fnn = fnn), 1:6,
                       info = sprintf("offset %g", by))
    }
  }
})

test_that("cell_values_at() gives the same cells at offsets 0, 1e8 and 1e12", {
  # The locations asked about: the lattice, on the polygons' edges and
  # 250 m from them, and the points moved 30 m east, some of them metres
  # from an edge.
  at <- function(by) {
    q <- .em_points(by)
    sf::st_geometry(q) <- sf::st_geometry(q) + c(30, 0)
    sf::st_crs(q) <- 32617
    rbind(.em_lattice(by)[, "pid"], q[, "pid"])
  }
  one <- function(by) {
    pts <- suppressWarnings(assign_features_to_polygons(.em_points(by), .em_polygons(by),
                                                        quiet = TRUE))
    vals <- .em_run(summarize_by_cell(pts, "y", id_col = "poly_id",
                                      cells_sf = .em_polygons(by), estimate = "mean"))$value
    .em_run(cell_values_at(at(by), vals, pts, "y", id_col = "poly_id"))$value
  }
  ref <- one(0)
  for (by in .em_offsets[-1]) {
    got <- one(by)
    expect_identical(got$poly_id, ref$poly_id, info = sprintf("offset %g", by))
    expect_identical(got$filled_from, ref$filled_from, info = sprintf("offset %g", by))
  }
})

test_that("cv_block_size_sweep() uses the same blocks at offsets 0, 1e8 and 1e12", {
  # lm() on the response and the two predictors (helper-lmfit.R): the fit
  # does not read the coordinates.
  lm_fit <- function(train_sf) lm_spatial_fit(train_sf, "y", c("x1", "x2"))
  sweep_at <- function(by) .em_run(cv_block_size_sweep(
    .em_points(by), "y", c("x1", "x2"), fit_fn = lm_fit, block_sizes = c(400, 1200),
    k = 3, quiet = TRUE, seed = 1))$value
  ref <- sweep_at(0)
  for (by in .em_offsets[-1]) {
    got <- sweep_at(by)
    expect_identical(got$blocks_used, ref$blocks_used, info = sprintf("offset %g", by))
  }
})
