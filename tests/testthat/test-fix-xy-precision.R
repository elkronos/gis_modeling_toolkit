# tests/testthat/test-fix-xy-precision.R
# ---------------------------------------------------------------------------
# Coordinates so far from their origin that doubles there do not hold a
# millionth of the layer's extent.
#
# BUG-226: at 1e15 doubles are 0.125 apart, and make_folds() lays block
# edges at that size, a corner plus a multiple of the block's side, so a
# point near an edge could fall in the next block with nothing said.
# make_folds() (block_kfold, buffered_loo), assign_features_to_polygons()
# and create_grid_polygons() now warn once, with class
# spatialkit_xy_imprecise, where that spacing is more than a millionth of
# the extent.
# BUG-270: clip_target_for() with no boundary drew the buffer around points
# on one line or at one location at that rounding (a 1 m disk of area
# 2.6875 at 1e15) and returned an empty LINESTRING where the radius was
# below it.  It now warns with the same class, and refuses an empty target
# with class spatialkit_xy_imprecise_error.
#
# Every coordinate is a multiple of 1/8, so the shifts are exact: at 1e15
# the layers are the same geometry, and only the arithmetic at that size
# differs.
# ---------------------------------------------------------------------------

.xp_x0 <- 600000
.xp_y0 <- 3990000

# 60 points over 4000 m x 3000 m, from a fixed sequence.
.xp_points <- function(by = 0, crs = 32617) {
  i <- seq_len(60)
  px <- .xp_x0 + round(8 * 4000 * ((i * 0.6180339887498949) %% 1)) / 8
  py <- .xp_y0 + round(8 * 3000 * ((i * 0.7548776662466927) %% 1)) / 8
  sf::st_as_sf(data.frame(pid = i, y = sin(i), x1 = cos(3 * i), X = px + by, Y = py + by),
               coords = c("X", "Y"), crs = crs)
}

# 12 rectangles of 1500 m x 1000 m tiling 4500 m x 4000 m around them.
.xp_polygons <- function(by = 0, crs = 32617) {
  cells <- list()
  for (j in 0:3) for (i in 0:2) {
    xa <- .xp_x0 - 250 + 1500 * i + by
    ya <- .xp_y0 - 500 + 1000 * j + by
    cells[[length(cells) + 1L]] <- sf::st_polygon(list(matrix(
      c(xa, ya, xa + 1500, ya, xa + 1500, ya + 1000, xa, ya + 1000, xa, ya),
      ncol = 2L, byrow = TRUE)))
  }
  sf::st_sf(poly_id = sprintf("P%02d", seq_along(cells)),
            geometry = sf::st_sfc(cells, crs = crs))
}

.xp_boundary <- function(by = 0) sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(
  xmin = .xp_x0 - 250 + by, ymin = .xp_y0 - 500 + by,
  xmax = .xp_x0 + 4250 + by, ymax = .xp_y0 + 3500 + by), crs = sf::st_crs(32617))))

.xp_run <- function(expr) {
  w <- list()
  value <- withCallingHandlers(suppressMessages(expr), warning = function(c) {
    w[[length(w) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

.xp_of <- function(r) Filter(function(c) inherits(c, "spatialkit_xy_imprecise"), r$warnings)

test_that("make_folds(), the join and the grid warn once where doubles do not hold the layer", {
  r <- .xp_run(make_folds(.xp_points(1e15), k = 4, method = "block_kfold", seed = 1))
  w <- .xp_of(r)
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1]]), paste0(
    "make_folds(): `points_sf` lies so far from the origin of its CRS that its ",
    "coordinates lose their precision: they are rounded to 0.125, more than a ",
    "millionth of its extent (", format(signif(w[[1]]$extent, 3)), "), so edges and ",
    "ties are read at that rounding. Shift the layer nearer the origin, or project it."))
  expect_identical(w[[1]]$spacing, 0.125)
  # The extent is the diagonal of the points' box, a few kilometres.
  xy <- sf::st_coordinates(.xp_points())
  expect_equal(w[[1]]$extent, sqrt(diff(range(xy[, 1]))^2 + diff(range(xy[, 2]))^2),
               tolerance = 1e-12)
  r <- .xp_run(make_folds(.xp_points(1e15), method = "buffered_loo", buffer = 400))
  expect_length(.xp_of(r), 1L)
  r <- .xp_run(assign_features_to_polygons(.xp_points(1e15), .xp_polygons(1e15)))
  w <- .xp_of(r)
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1]]), paste0(
    "assign_features_to_polygons(): `features_sf` and `polygons_sf` lie so far from the ",
    "origin of their CRS that their coordinates lose their precision: they are rounded ",
    "to 0.125, more than a millionth of their extent ("), fixed = TRUE)
  r <- .xp_run(create_grid_polygons(.xp_boundary(1e15), target_cells = 12, quiet = TRUE))
  w <- .xp_of(r)
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1]]), paste0(
    "create_grid_polygons(): `boundary` lies so far from the origin of its CRS that its ",
    "coordinates lose their precision: they are rounded to 0.125"), fixed = TRUE)
  expect_gt(nrow(r$value), 0L)
  # With no CRS the message names no CRS.
  r <- .xp_run(make_folds(.xp_points(1e15, crs = NA), k = 4, method = "block_kfold", seed = 1))
  expect_match(conditionMessage(.xp_of(r)[[1]]),
               "`points_sf` lies so far from the origin that its coordinates", fixed = TRUE)
})

test_that("no such warning where doubles hold the layer, and the results are the origin's", {
  f0 <- .xp_run(make_folds(.xp_points(0), k = 4, method = "block_kfold", seed = 1))$value
  a0 <- .xp_run(assign_features_to_polygons(.xp_points(0), .xp_polygons(0)))$value
  g0 <- .xp_run(create_grid_polygons(.xp_boundary(0), target_cells = 12, quiet = TRUE))$value
  for (by in c(0, 1e8, 1e12)) {
    r <- .xp_run(make_folds(.xp_points(by), k = 4, method = "block_kfold", seed = 1))
    expect_length(.xp_of(r), 0L)
    expect_identical(r$value$assignment$fold, f0$assignment$fold, info = sprintf("offset %g", by))
    r <- .xp_run(make_folds(.xp_points(by), method = "buffered_loo", buffer = 400))
    expect_length(.xp_of(r), 0L)
    r <- .xp_run(assign_features_to_polygons(.xp_points(by), .xp_polygons(by)))
    expect_length(.xp_of(r), 0L)
    expect_identical(r$value$poly_id, a0$poly_id, info = sprintf("offset %g", by))
    r <- .xp_run(create_grid_polygons(.xp_boundary(by), target_cells = 12, quiet = TRUE))
    expect_length(.xp_of(r), 0L)
    expect_identical(nrow(r$value), nrow(g0))
  }
  # A layer at one location has no extent to resolve: nothing is said of
  # it, at any offset.
  one <- .xp_points(1e15)[c(1, 1), ]
  r <- .xp_run(spatialkit:::.warn_xy_imprecise(one, "f", "x"))
  expect_length(.xp_of(r), 0L)
  expect_false(r$value$warned)
})

test_that("a block edge rounded onto a point at 1e15 is said", {
  # Points on x in [0, 3000.25], block_nx = 3: the first edge lies at
  # 1000.0833 near the origin and is held as 1000.125 at 1e15, on the point
  # at x = 1000.125, which then falls in the block to its left.
  x <- c(0, 3000.25, 1000.125, (1:37) * 81)
  y <- c(0, 1000, 500, ((1:37) * 27) %% 1000)
  mk <- function(by) sf::st_as_sf(data.frame(pid = seq_along(x), X = x + by, Y = y + by),
                                  coords = c("X", "Y"), crs = 32617)
  r <- .xp_run(make_folds(mk(1e15), k = 3, method = "block_kfold", seed = 1,
                          block_nx = 3, block_ny = 1))
  expect_length(.xp_of(r), 1L)
  r <- .xp_run(make_folds(mk(1e12), k = 3, method = "block_kfold", seed = 1,
                          block_nx = 3, block_ny = 1))
  expect_length(.xp_of(r), 0L)
  # At 1e12 doubles hold the edge to 1.2e-4 m, and the point 0.042 m right
  # of it is in the block it is in near the origin: the margin there is a
  # millionth of the layer's diagonal (3 mm), not 100 m (BUG-225).
  r0 <- .xp_run(make_folds(mk(0), k = 3, method = "block_kfold", seed = 1,
                           block_nx = 3, block_ny = 1))
  expect_identical(r$value$assignment$fold, r0$value$assignment$fold)
})

test_that("a run of folds on one layer says it once", {
  r <- .xp_run(cv_block_size_sweep(
    .xp_points(1e15), "y", "x1",
    fit_fn = function(tr) lm_spatial_fit(tr, "y", "x1"), block_sizes = c(800, 1600), k = 3,
    quiet = TRUE, seed = 1))
  expect_length(.xp_of(r), 1L)
})

test_that("the warning's bound is a millionth of the extent", {
  # The 60 points span a diagonal of about 4730 m, a millionth of which is
  # 4.7 mm.  Doubles are 2^-8 (3.9 mm) apart at 2e13 and 2^-7 (7.8 mm) apart
  # at 4e13, so the layer is said at 4e13 and not at 2e13.  Every
  # coordinate is a multiple of 1/8, so both shifts are exact.
  r <- .xp_run(make_folds(.xp_points(2e13), k = 4, method = "block_kfold", seed = 1))
  expect_length(.xp_of(r), 0L)
  r <- .xp_run(make_folds(.xp_points(4e13), k = 4, method = "block_kfold", seed = 1))
  w <- .xp_of(r)
  expect_length(w, 1L)
  expect_identical(w[[1]]$spacing, 2^-7)
  # On a box: at 1e15 doubles are 0.125 apart; an extent of 0.125 / 0.9e-6
  # holds them, one of 0.125 / 1.1e-6 does not.
  box <- function(d) c(xmin = 1e15, ymin = 1e15, xmax = 1e15 + d, ymax = 1e15)
  r <- .xp_run(spatialkit:::.warn_xy_imprecise(box(0.125 / 0.9e-6), "f", "x"))
  expect_length(.xp_of(r), 0L)
  expect_false(r$value$warned)
  r <- .xp_run(spatialkit:::.warn_xy_imprecise(box(0.125 / 1.1e-6), "f", "x"))
  expect_length(.xp_of(r), 1L)
  expect_true(r$value$warned)
})

test_that("the Voronoi cells and triangles warn where doubles do not hold the points", {
  # Point 40 lies 0.25 m east of point 1.  At 1e15 doubles are 0.125 apart
  # and the margin of an edge is 0.5 m (.edge_tol()), so both points fall
  # in one cell; near the origin each has its own.
  mk <- function(by) {
    p <- .xp_points(by)[1:40, ]
    xy <- sf::st_coordinates(p)
    xy[40, ] <- xy[1, ] + c(0.25, 0)
    sf::st_set_geometry(p, sf::st_sfc(lapply(seq_len(40), function(i) sf::st_point(xy[i, ])),
                                      crs = 32617))
  }
  for (by in c(0, 1e8, 1e12)) {
    r <- .xp_run(create_voronoi_polygons(mk(by), quiet = TRUE))
    expect_length(.xp_of(r), 0L)
    expect_length(unique(r$value$index), 40L)
    r <- .xp_run(build_tessellation(mk(by), method = "voronoi", quiet = TRUE))
    expect_length(.xp_of(r), 0L)
    expect_length(unique(r$value$index), 40L)
    r <- .xp_run(build_tessellation(mk(by), method = "triangles", quiet = TRUE))
    expect_length(.xp_of(r), 0L)
  }
  r <- .xp_run(create_voronoi_polygons(mk(1e15), quiet = TRUE))
  w <- .xp_of(r)
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1]]), paste0(
    "create_voronoi_polygons(): `points_sf` lies so far from the origin of its CRS that its ",
    "coordinates lose their precision: they are rounded to 0.125, more than a millionth of ",
    "its extent (", format(signif(w[[1]]$extent, 3)), "), so the cells' edges and the ",
    "points' cells are read at that rounding. Shift the layer nearer the origin, or project it."))
  for (m in c("voronoi", "triangles")) {
    r <- .xp_run(build_tessellation(mk(1e15), method = m, quiet = TRUE))
    w <- .xp_of(r)
    expect_length(w, 1L)
    expect_match(conditionMessage(w[[1]]), "^build_tessellation\\(\\): `points_sf` lies so far")
  }
  # The grids say it once, of the boundary, through create_grid_polygons().
  for (m in c("hex", "square")) {
    r <- .xp_run(build_tessellation(mk(1e15), method = m, boundary = .xp_boundary(1e15),
                                    approx_n_cells = 12, quiet = TRUE))
    expect_length(.xp_of(r), 1L)
  }
})

# Points at one location, and two on a line 20 m and 0.5 m apart.
.xp_line <- function(by, d) sf::st_as_sf(data.frame(x = c(by, by + d), y = c(by, by)),
                                         coords = c("x", "y"), crs = 32617)

# The target clip_target_for() draws for points on one line: the buffer of
# their union by `r`, made valid, which the tests compare to the bit.
.xp_disk <- function(p, r) spatialkit:::.safe_make_valid(suppressWarnings(
  sf::st_buffer(sf::st_union(sf::st_geometry(p)), dist = r)))

test_that("clip_target_for() says where its buffer is drawn at the rounding of the coordinates", {
  cases <- list(list(d = 0, r = 1), list(d = 20, r = 1))
  for (cs in cases) {
    for (by in c(0, 1e8)) {
      p <- .xp_line(by, cs$d)
      r <- .xp_run(clip_target_for(p))
      expect_length(.xp_of(r), 0L)
      expect_identical(sf::st_geometry(r$value), .xp_disk(p, cs$r),
                       info = sprintf("d %g, offset %g", cs$d, by))
    }
    for (by in c(1e12, 1e15)) {
      p <- .xp_line(by, cs$d)
      r <- .xp_run(clip_target_for(p, quiet = TRUE))
      w <- .xp_of(r)
      expect_length(w, 1L)
      u <- 2^(floor(log2(by)) - 52)
      expect_identical(conditionMessage(w[[1]]), paste0(
        "clip_target_for(): `points_sf` lies so far from the origin of its CRS that its ",
        "coordinates lose their precision: they are rounded to ", format(signif(u, 3)),
        ", more than a millionth of the target's extent (", format(cs$d + 2 * cs$r),
        "), so the buffer of radius ", cs$r, " around the points is drawn at that ",
        "rounding. Shift the layer nearer the origin, or pass `boundary`."),
        info = sprintf("d %g, offset %g", cs$d, by))
      # The target is still the buffer, drawn at that rounding.
      expect_identical(sf::st_geometry(r$value), .xp_disk(p, cs$r),
                       info = sprintf("d %g, offset %g", cs$d, by))
    }
  }
})

test_that("clip_target_for() refuses a buffer whose radius is below the rounding", {
  # Two points 0.5 m apart: a radius of 0.025, under the 0.125 spacing of
  # doubles at 1e15, where the buffer has no area.
  err <- tryCatch(suppressWarnings(clip_target_for(.xp_line(1e15, 0.5))),
                  error = function(e) e)
  expect_s3_class(err, "spatialkit_xy_imprecise_error")
  expect_match(conditionMessage(err), paste0(
    "clip_target_for(): `points_sf` lies so far from the origin of its CRS that its ",
    "coordinates lose their precision: they are rounded to 0.125, more than the radius ",
    "0.025 of the buffer around the points, which comes back "), fixed = TRUE)
  expect_match(conditionMessage(err), paste0(
    "with no area, so no clip target can be drawn. Shift the layer nearer the origin, ",
    "or pass `boundary`."), fixed = TRUE)
  # Near the origin the same points have their buffer, and nothing is said.
  p <- .xp_line(0, 0.5)
  r <- .xp_run(clip_target_for(p))
  expect_length(.xp_of(r), 0L)
  expect_identical(sf::st_geometry(r$value), .xp_disk(p, 0.025))
})

test_that("clip_target_for() keeps a buffer below the rounding that has an area", {
  # Two points 2 m apart at 1e15: a radius of 0.1, under the 0.125 spacing.
  # The buffer is drawn at that rounding; the call stops only where it comes
  # back empty or with no area, as the help says.
  r <- tryCatch(.xp_run(clip_target_for(.xp_line(1e15, 2))), error = function(e) e)
  if (inherits(r, "error")) {
    expect_s3_class(r, "spatialkit_xy_imprecise_error")
  } else {
    expect_length(.xp_of(r), 1L)
    expect_true(all(as.character(sf::st_geometry_type(r$value)) %in%
                      c("POLYGON", "MULTIPOLYGON")))
    expect_gt(as.numeric(sum(sf::st_area(r$value))), 0)
  }
})

test_that("clip_target_for() says nothing of a box drawn from the points' own corners", {
  r <- .xp_run(clip_target_for(.xp_points(1e15)))
  expect_length(.xp_of(r), 0L)
  expect_identical(as.character(sf::st_geometry_type(r$value)), "POLYGON")
})

test_that("the help pages and NEWS say so", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  for (page in c("make_folds.Rd", "assign_features_to_polygons.Rd", "create_grid_polygons.Rd",
                 "clip_target_for.Rd"))
    expect_match(read("man", page), "spatialkit_xy_imprecise", fixed = TRUE, info = page)
  expect_match(read("man", "clip_target_for.Rd"), "spatialkit_xy_imprecise_error", fixed = TRUE)
  # The refusal is of a buffer that comes back empty or with no area, not of
  # every radius below the spacing (a radius of 0.1 at 1e15 has an area).
  expect_match(read("man", "clip_target_for.Rd"), paste0(
    "Where the buffer comes back empty or with no area, as a radius below that ",
    "spacing can make it, the call stops"), fixed = TRUE)
  expect_no_match(read("man", "clip_target_for.Rd"), "where the radius is below the spacing",
                  fixed = TRUE)
  for (page in c("create_voronoi_polygons.Rd", "build_tessellation.Rd"))
    expect_match(read("man", page), "spatialkit_xy_imprecise", fixed = TRUE, info = page)
  news <- read("NEWS.md")
  expect_match(news, "Coordinates so far from their origin that doubles cannot hold the layer",
               fixed = TRUE)
  expect_match(news, paste0("`assign_features_to_polygons()`, `create_grid_polygons()`, ",
                            "`create_voronoi_polygons()` and `build_tessellation()`"), fixed = TRUE)
  expect_match(news, paste0("`clip_target_for()` drew its buffer around points on one line or ",
                            "at one location at the rounding of their coordinates"), fixed = TRUE)
})
