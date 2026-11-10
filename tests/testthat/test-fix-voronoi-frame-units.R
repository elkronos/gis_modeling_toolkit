# tests/testthat/test-fix-voronoi-frame-units.R
# ---------------------------------------------------------------------------
# Two frames that depended on where the layer sits or on its unit.
#
# create_voronoi_polygons() and build_tessellation(method = "voronoi") with
# no boundary clip the cells to the points' convex hull grown by 2 percent of
# its diagonal, and a cover the same way, as their help says.  The margin was
# 0.02 * max(diagonal, 1): a layer whose diagonal is under 1 CRS unit (a
# 500 m layer in a kilometre CRS) grew by 0.02 units, 20 m, where in metres
# it grew by 15.9 m, and its edge cells took another shape.  The floor is
# gone (.voronoi_grow()); a diagonal of 0 still grows by 0.02 units.
#
# clip_target_for() returns a buffer around the points in place of their
# bounding box when the box's short side is under a millionth of its long
# side.  It also did when all.equal() called the two ends of a side equal,
# which compares them relative to their size: a 3 km x 2 km box at
# x = 1e12, or a 5 cm x 3 cm one at y = 4e6, got disks around its points.
# The side is now read against the span only.
# ---------------------------------------------------------------------------

# 9 points spread over a 500 m square in EPSG:32617, off a regular lattice so
# the Voronoi edges are not tied.
.vfu_pts <- function() {
  d <- data.frame(x = 5e5 + c(0, 250, 500, 0, 250, 500, 0, 250, 500) + c(0, 3, 7, 1, 0, 2, 5, 4, 6),
                  y = 4e6 + c(0, 0, 0, 250, 250, 250, 500, 500, 500) + c(2, 0, 1, 6, 3, 0, 0, 2, 4))
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
}

# The 8 x 8 lattice at 250 m, every distance tied.
.vfu_lattice <- function() {
  g <- expand.grid(x = 5e5 + 250 * 0:7, y = 4e6 + 250 * 0:7)
  sf::st_as_sf(g, coords = c("x", "y"), crs = 32617)
}

# The layer with every coordinate multiplied by `by`, in a CRS whose unit is
# 1 / by metres.
.vfu_scale <- function(x, by) {
  p4  <- trimws(gsub("\\+units=[^ ]+", "", sf::st_crs(x)$proj4string))
  crs <- sf::st_crs(paste(p4, sprintf("+to_meter=%.17g", 1 / by)))
  g <- sf::st_geometry(x) * by
  sf::st_geometry(x) <- sf::st_set_crs(g, crs)
  x
}

# The same ground in UTM 17N in kilometres, a unit sf names.
.vfu_km <- function(x) sf::st_transform(x, "+proj=utm +zone=17 +datum=WGS84 +units=km +no_defs")

# Each cell's area and centroid, in metres, in the order of cell_id.
.vfu_cells <- function(cells, by) {
  cells <- cells[order(cells$cell_id), ]
  g <- sf::st_set_crs(sf::st_geometry(cells), NA)
  list(area = as.numeric(sf::st_area(g)) / by^2,
       xy = unname(sf::st_coordinates(sf::st_centroid(g))) / by)
}

# The relative tolerance of the comparison with the metre layer: the scaled
# layer differs by one rounding of the multiplication by `by` (about 1e-16 of
# coordinates near 4e6 m) and GEOS's Voronoi and clip work from those
# coordinates, so 1e-8 is far above the rounding and far below the 3 percent
# the floor moved the edge cells' areas by.
.vfu_expect_same_cells <- function(got, ref, by) {
  a <- .vfu_cells(got, by); b <- .vfu_cells(ref, 1)
  expect_equal(a$area, b$area, tolerance = 1e-8)
  expect_equal(a$xy, b$xy, tolerance = 1e-8)
}

# The margin the frame was grown by beyond the points' box, in metres.  The
# buffer draws a round corner with straight segments (30 to a quarter
# circle), so where the box's left edge meets a corner of the hull the
# margin falls short of the buffer distance by up to 1 - cos(pi / 120),
# about 3.4e-4 of it: margins are compared to 1e-3 relative.  The floor in
# map units made them differ by a quarter and more.
.vfu_margin <- function(boundary, pts, by) {
  bb <- sf::st_bbox(boundary); pb <- sf::st_bbox(pts)
  as.numeric(pb[["xmin"]] - bb[["xmin"]]) / by
}

test_that("the frame of Voronoi cells grows by 2% of the diagonal in a unit sf does not name", {
  for (fx in list(.vfu_pts(), .vfu_lattice())) {
    ref <- create_voronoi_polygons(fx, quiet = TRUE)
    pb <- sf::st_bbox(fx)
    diag_m <- sqrt(as.numeric(pb[["xmax"]] - pb[["xmin"]])^2 + as.numeric(pb[["ymax"]] - pb[["ymin"]])^2)
    expect_equal(.vfu_margin(ref$boundary, fx, 1), 0.02 * diag_m, tolerance = 1e-3)
    for (by in c(1e-6, 1e-3, 0.304801)) {
      s <- .vfu_scale(fx, by)
      got <- create_voronoi_polygons(s, quiet = TRUE)
      expect_equal(.vfu_margin(got$boundary, s, by), 0.02 * diag_m, tolerance = 1e-3)
      .vfu_expect_same_cells(got$cells, ref$cells, by)
      bt <- build_tessellation(s, method = "voronoi", quiet = TRUE)
      .vfu_expect_same_cells(bt$cells, ref$cells, by)
    }
  }
})

test_that("a 500 m layer in kilometres gets the cells it gets in metres", {
  pm <- .vfu_pts()
  ref <- create_voronoi_polygons(pm, quiet = TRUE)
  km  <- create_voronoi_polygons(.vfu_km(pm), quiet = TRUE)
  expect_equal(sum(.vfu_cells(km$cells, 1e-3)$area), sum(.vfu_cells(ref$cells, 1)$area),
               tolerance = 1e-6)
  expect_equal(.vfu_margin(km$boundary, .vfu_km(pm), 1e-3),
               .vfu_margin(ref$boundary, pm, 1), tolerance = 1e-3)
})

test_that("a cover grows by 2% of its diagonal in a unit sf does not name", {
  pts <- .vfu_pts()
  cover <- sf::st_as_sf(data.frame(x = 5e5 + c(-100, 650, -100, 650), y = 4e6 + c(-80, -80, 600, 600)),
                        coords = c("x", "y"), crs = 32617)
  ref <- create_voronoi_polygons(pts, cover = cover, quiet = TRUE)
  diag_m <- sqrt(750^2 + 680^2)
  expect_equal(.vfu_margin(ref$boundary, cover, 1), 0.02 * diag_m, tolerance = 1e-3)
  for (by in c(1e-6, 0.304801)) {
    got <- create_voronoi_polygons(.vfu_scale(pts, by), cover = .vfu_scale(cover, by), quiet = TRUE)
    expect_equal(.vfu_margin(got$boundary, .vfu_scale(cover, by), by), 0.02 * diag_m, tolerance = 1e-3)
    .vfu_expect_same_cells(got$cells, ref$cells, by)
  }
})

test_that("a frame whose diagonal is over 1 unit is buffered by 0.02 times it, as before", {
  # 0.02 * max(diag, 1) and 0.02 * diag are the same double when diag >= 1.
  pts <- .vfu_pts()
  got <- create_voronoi_polygons(pts, quiet = TRUE)
  hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))
  bb <- sf::st_bbox(hull)
  diag <- sqrt((bb$xmax - bb$xmin)^2 + (bb$ymax - bb$ymin)^2)
  want <- suppressWarnings(sf::st_buffer(hull, dist = 0.02 * max(diag, 1)))
  expect_identical(sf::st_bbox(got$boundary), sf::st_bbox(want))
  expect_identical(.voronoi_grow(as.numeric(diag)), 0.02 * max(as.numeric(diag), 1))
})

test_that("a diagonal of 0 still grows by 0.02 units, and a diagonal that is not finite stays so", {
  expect_identical(.voronoi_grow(0), 0.02)
  expect_identical(.voronoi_grow(0.5), 0.01)
  expect_true(is.na(.voronoi_grow(NA_real_)))
  expect_identical(.voronoi_grow(Inf), Inf)
})

test_that("clip_target_for() gives the points' box at any offset", {
  base <- data.frame(x = c(0, 3000, 1500), y = c(0, 500, 2000))
  for (expand in c(0, 0.1, 500)) {
    ref <- clip_target_for(sf::st_as_sf(base + 0, coords = c("x", "y"), crs = 32617),
                           expand = expand, quiet = TRUE)
    ref_bb <- as.numeric(sf::st_bbox(ref))
    ref_a  <- as.numeric(sf::st_area(ref))
    for (off in c(1e8, 1e12, 1e15)) {
      p <- sf::st_as_sf(data.frame(x = base$x + off, y = base$y + off), coords = c("x", "y"),
                        crs = 32617)
      got <- clip_target_for(p, expand = expand, quiet = TRUE)
      # The spacing of doubles at 1e15 is 0.125, so the box is held to 1 m
      # (the defect shrank it by 186 m) and the area to 1e-4 of it; at 1e12
      # and 1e8 both come out far closer.
      expect_lt(max(abs(as.numeric(sf::st_bbox(got)) - off - ref_bb)), 1)
      expect_equal(as.numeric(sf::st_area(got)), ref_a, tolerance = 1e-4,
                   info = sprintf("offset %g, expand %g", off, expand))
    }
  }
})

test_that("clip_target_for() gives the box of a 5 cm x 3 cm layer at ordinary coordinates", {
  p <- sf::st_as_sf(data.frame(x = 5e5 + c(0, 0.05, 0.02), y = 4e6 + c(0, 0.01, 0.03)),
                    coords = c("x", "y"), crs = 32617)
  got <- clip_target_for(p, expand = 0, quiet = TRUE)
  expect_equal(as.numeric(sf::st_area(got)), 0.05 * 0.03, tolerance = 1e-6)
  expect_identical(sf::st_bbox(got), sf::st_bbox(p))
})

test_that("clip_target_for() still buffers points on a line, on a near line and at one place", {
  line <- sf::st_as_sf(data.frame(x = 1e12 + c(0, 1000, 2000), y = c(5, 5, 5)),
                       coords = c("x", "y"), crs = 32617)
  near <- sf::st_as_sf(data.frame(x = c(0, 1000, 2000), y = c(0, 1e-4, 0)),
                       coords = c("x", "y"), crs = 32617)
  one  <- sf::st_as_sf(data.frame(x = c(10, 10), y = c(20, 20)), coords = c("x", "y"), crs = 32617)
  for (p in list(line, near, one)) {
    expect_message(got <- clip_target_for(p), "degenerate bbox")
    expect_gt(as.numeric(sf::st_area(got)), 0)
    expect_identical(as.character(sf::st_geometry_type(got)) %in% c("POLYGON", "MULTIPOLYGON"), TRUE)
  }
})

test_that("clip_target_for() reads a box as degenerate below a millionth of its long side, at any offset", {
  # A short side of 1e-5 of the long side is a box; one of 5e-7 is a line.
  # At 1e12 the spacing of doubles is 1.2e-4, which moves the short sides
  # (0.02 and 0.001) by under 3 percent, so each ratio stays on its side of
  # 1e-6.
  for (off in c(0, 1e12)) {
    box <- sf::st_as_sf(data.frame(x = off + c(0, 2000, 1000), y = off + c(0, 0, 0.02)),
                        coords = c("x", "y"), crs = 32617)
    expect_no_message(got <- clip_target_for(box), message = "degenerate")
    expect_identical(sf::st_bbox(got), sf::st_bbox(box), info = sprintf("offset %g", off))
    thin <- sf::st_as_sf(data.frame(x = off + c(0, 2000, 1000), y = off + c(0, 0, 0.001)),
                         coords = c("x", "y"), crs = 32617)
    expect_message(got <- clip_target_for(thin), "degenerate bbox", info = sprintf("offset %g", off))
    # The buffer is 5 percent of the 2000 m extent, 100 m, around the points,
    # so the target is 200 m high where the box is 0.001 m.
    expect_gt(as.numeric(sf::st_bbox(got)$ymax - sf::st_bbox(got)$ymin), 199)
  }
})

test_that("clip_target_for() with no boundary refuses points none of which has a location", {
  # The box of such points is NA; the call stopped with R's "missing value
  # where TRUE/FALSE needed".
  why_given <- "no point in `points_sf` has a location \\(every row is empty, or has a coordinate"
  empty <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(), sf::st_point(), crs = 32617))
  expect_error(clip_target_for(empty, quiet = TRUE), why_given)
  infs <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(Inf, 1)), sf::st_point(c(1, NaN)),
                                          sf::st_point(c(2e200, 0)), crs = 32617))
  expect_error(clip_target_for(infs, quiet = TRUE), why_given)
  expect_error(clip_target_for(sf::st_set_crs(infs, NA), quiet = TRUE), why_given)
  ll_empty <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(), sf::st_point(), crs = 4326))
  expect_error(clip_target_for(ll_empty, quiet = TRUE), why_given)
  # In degrees beyond the pole a point has a location as given, and the move
  # to a projected CRS empties it.
  pole <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(10, 95)), sf::st_point(c(11, 96)),
                                          crs = 4326))
  expect_error(suppressWarnings(clip_target_for(pole, quiet = TRUE)),
               "no point in `points_sf` has a location \\(every row is empty after its move")
  # One located point is a degenerate box and gets its buffer, as before;
  # with a boundary the points take no part, and the boundary comes back.
  one <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(), sf::st_point(c(5e5, 4e6)), crs = 32617))
  expect_gt(as.numeric(sf::st_area(clip_target_for(one, quiet = TRUE))), 0)
  b <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10),
                                                                    c(0, 10), c(0, 0)))),
                                       crs = 32617))
  expect_equal(as.numeric(sf::st_area(clip_target_for(empty, boundary = b, quiet = TRUE))), 100)
})

test_that("create_voronoi_polygons() on one location frames a disk of 0.02 units", {
  p <- sf::st_as_sf(data.frame(x = 5e5, y = 4e6), coords = c("x", "y"), crs = 32617)
  got <- create_voronoi_polygons(p, quiet = TRUE)
  # The disk is drawn with 30 segments a quarter circle, whose area is short
  # of pi r^2 by 4.6e-4 of it.
  expect_equal(as.numeric(sf::st_area(got$cells)), pi * 0.02^2, tolerance = 1e-3)
  expect_equal(as.numeric(sf::st_area(got$boundary)), pi * 0.02^2, tolerance = 1e-3)
})
