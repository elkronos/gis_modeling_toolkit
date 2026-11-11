# tests/testthat/test-fix-line-midpoint-units.R
# ---------------------------------------------------------------------------
# coerce_to_points() modes "auto" and "line_midpoint" (and prep_model_data(),
# which pointizes lines with them) took a line's midpoint with
# sf::st_line_sample(sample = 0.5).  sf measures the line with st_length(),
# in metres when the CRS has a +to_meter whose unit sf cannot name, and GDAL
# reads half that length as a distance along the line in CRS units.  In a
# unit of 2 m the point was the line's end, in a unit of 0.5 m a quarter of
# the way along.  The midpoint is now taken at half the planar length in the
# CRS's own unit (.planar_line_midpoints()); lines in metres and in a named
# unit (EPSG:2264, US survey feet) get the same numbers as before, to the bit.
#
# The scaled copies multiply every coordinate by `by` and state the unit as
# +to_meter = 1 / by, so they describe the same ground: by = 1e-6 (a unit of
# 1000 km), by = 0.304801 (a unit of about 3.28 m), by = 0.5 (a unit of 2 m)
# and by = 2 (a unit of 0.5 m).
# ---------------------------------------------------------------------------

# Bent lines with segments of unequal length, so half the length falls inside
# a segment and not on a vertex, and a MULTILINESTRING whose longest part is
# not its first.
.lmu_lines <- function(crs = 32617) {
  x0 <- 5e5; y0 <- 4e6
  l1 <- sf::st_linestring(cbind(x0 + c(0, 120, 400, 410), y0 + c(0, 90, 60, 700)))
  l2 <- sf::st_linestring(cbind(x0 + c(1000, 1300, 1310), y0 + c(50, 450, 455)))
  l3 <- sf::st_linestring(cbind(x0 + c(2000, 2003, 2600), y0 + c(0, 1, 30)))
  sf::st_sf(id = 1:3, geometry = sf::st_sfc(l1, l2, l3, crs = crs))
}

.lmu_mlines <- function(crs = 32617) {
  x0 <- 5e5; y0 <- 4e6
  m1 <- sf::st_multilinestring(list(
    cbind(x0 + c(0, 50), y0 + c(0, 0)),
    cbind(x0 + c(100, 400, 450), y0 + c(100, 500, 520))))
  m2 <- sf::st_multilinestring(list(
    cbind(x0 + c(900, 905, 1500), y0 + c(0, 2, 40)),
    cbind(x0 + c(900, 930), y0 + c(100, 100))))
  sf::st_sf(id = 1:2, geometry = sf::st_sfc(m1, m2, crs = crs))
}

# The layer with every coordinate multiplied by `by`, in a CRS whose unit is
# 1 / by metres.
.lmu_scale <- function(x, by) {
  p4  <- trimws(gsub("\\+units=[^ ]+", "", sf::st_crs(x)$proj4string))
  crs <- sf::st_crs(paste(p4, sprintf("+to_meter=%.17g", 1 / by)))
  g <- sf::st_geometry(x) * by
  sf::st_geometry(x) <- sf::st_set_crs(g, crs)
  x
}

.lmu_xy <- function(x) unname(sf::st_coordinates(x)[, 1:2, drop = FALSE])

# The relative tolerance of the comparison of a scaled layer's midpoints with
# the metre layer's scaled: the two are computed from coordinates that differ
# by one rounding of the multiplication by `by`, about 1e-16 of the
# coordinates (5e5 to 4e6 m), so 1e-9 of them is far above that rounding and
# far below the factor of the unit the defect put on the distance.
.lmu_expect_scaled <- function(got, ref, by) {
  expect_equal(.lmu_xy(got), .lmu_xy(ref) * by, tolerance = 1e-9)
}

test_that("the scaled CRSs have a unit sf does not name, so sf::st_length() is in metres", {
  for (by in c(1e-6, 0.304801, 0.5, 2)) {
    s <- .lmu_scale(.lmu_lines(), by)
    expect_null(sf::st_crs(s)$ud_unit)
    expect_identical(as.character(units(sf::st_length(s))), "m")
  }
})

test_that("coerce_to_points() takes a line's midpoint at half its length in a unit sf does not name", {
  ref <- .lmu_lines()
  for (mode in c("auto", "line_midpoint")) {
    ref_pt <- coerce_to_points(ref, mode = mode)
    for (by in c(1e-6, 0.304801, 0.5, 2)) {
      s <- .lmu_scale(ref, by)
      for (tp in c(TRUE, FALSE)) {
        got <- coerce_to_points(s, mode = mode, tmp_project = tp)
        expect_identical(sf::st_crs(got), sf::st_crs(s))
        .lmu_expect_scaled(got, ref_pt, by)
      }
    }
  }
})

test_that("a straight line gets its middle, not its end or a quarter of the way along", {
  base <- sf::st_crs(32617)$proj4string
  ln <- sf::st_linestring(rbind(c(0, 0), c(10, 0), c(20, 0)))
  for (tm in c(2, 0.5)) {
    crs <- sf::st_crs(paste(trimws(gsub("\\+units=[^ ]+", "", base)), sprintf("+to_meter=%g", tm)))
    x <- sf::st_sf(id = 1, geometry = sf::st_sfc(ln, crs = crs))
    for (mode in c("auto", "line_midpoint")) {
      expect_equal(.lmu_xy(coerce_to_points(x, mode = mode)), cbind(10, 0))
      expect_equal(.lmu_xy(coerce_to_points(x, mode = mode, tmp_project = FALSE)), cbind(10, 0))
    }
  }
  # EPSG:2314 (Trinidad Grid) is in Clarke's feet, a unit sf has no name for.
  x <- sf::st_sf(id = 1, geometry = sf::st_sfc(ln + c(3e5, 3e5), crs = 2314))
  expect_null(sf::st_crs(x)$ud_unit)
  for (mode in c("auto", "line_midpoint"))
    expect_equal(.lmu_xy(coerce_to_points(x, mode = mode)), cbind(3e5 + 10, 3e5))
})

test_that("a MULTILINESTRING's longest part gets its midpoint in a unit sf does not name", {
  ref <- .lmu_mlines()
  ref_pt <- coerce_to_points(ref, mode = "auto")
  for (by in c(1e-6, 0.304801, 0.5, 2)) {
    s <- .lmu_scale(ref, by)
    .lmu_expect_scaled(coerce_to_points(s, mode = "auto"), ref_pt, by)
    .lmu_expect_scaled(coerce_to_points(s, mode = "auto", tmp_project = FALSE), ref_pt, by)
  }
})

test_that("prep_model_data() pointizes lines at their midpoints in a unit sf does not name", {
  ref <- .lmu_lines()
  ref$y <- c(1.5, 2.5, 4)
  ref_pd <- prep_model_data(ref, "y", character(0), pointize = "line_midpoint")
  for (by in c(1e-6, 0.304801)) {
    s <- .lmu_scale(ref, by)
    for (pz in c("auto", "line_midpoint")) {
      got <- prep_model_data(s, "y", character(0), pointize = pz)
      .lmu_expect_scaled(got, ref_pd, by)
    }
  }
})

test_that("lines in metres and in US survey feet get sf::st_line_sample()'s midpoints to the bit", {
  # The midpoints are still taken by sf::st_line_sample(), with the CRS set
  # aside: in metres and in a unit sf names, st_length() gives the same
  # number with the CRS or without it, so the points are those of before.
  for (crs in c(32617, 2264)) {
    x <- .lmu_lines(crs)
    before <- sf::st_cast(sf::st_line_sample(sf::st_geometry(x), sample = 0.5), "POINT")
    for (mode in c("auto", "line_midpoint")) {
      for (tp in c(TRUE, FALSE)) {
        got <- coerce_to_points(x, mode = mode, tmp_project = tp)
        expect_identical(.lmu_xy(got), .lmu_xy(before))
      }
    }
    m <- .lmu_mlines(crs)
    parts <- lapply(sf::st_geometry(m), function(g) sf::st_sfc(lapply(unclass(g), sf::st_linestring), crs = crs))
    before_m <- do.call(rbind, lapply(parts, function(p) {
      k <- which.max(as.numeric(sf::st_length(p)))
      .lmu_xy(sf::st_cast(sf::st_line_sample(p[k], sample = 0.5), "POINT"))
    }))
    expect_identical(.lmu_xy(coerce_to_points(m, mode = "auto")), before_m)
  }
})
