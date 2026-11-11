# tests/testthat/test-fix-antimeridian-box-centre.R
# ---------------------------------------------------------------------------
# BUG-269.  coerce_to_points(mode = "bbox_center"), and mode "centroid" with
# s2 off, read a lon/lat geometry written across longitude 180 on the plane:
# the box of a square from 179.5 to -179.5 spans the globe, and its centre,
# and GEOS's centroid in degrees, lay at longitude 0, about 17,600 km from
# it, with nothing said.  prep_model_data(pointize = "bbox_center") then
# projected a layer of such squares as one that circles the North Pole.
# Such a geometry is now read with the longitudes below 0 taken from 180 to
# 360, as the points on the surface are (BUG-222).
# ---------------------------------------------------------------------------

.am_sq <- function(x0, x1, y0 = 10, y1 = 11)
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))

# Rows 1 to 3 cross longitude 180 (a square, a line and a MULTIPOINT); rows
# 4 and 5 (a square and a line at 170) do not.
.am_layer <- function(crs = 4326)
  sf::st_sf(id = 1:5, geometry = sf::st_sfc(
    .am_sq(179.5, -179.5),
    sf::st_linestring(rbind(c(179.5, 10), c(-179.5, 10.5))),
    sf::st_multipoint(rbind(c(179.5, 10), c(-179.5, 10.5))),
    .am_sq(170, 171),
    sf::st_linestring(rbind(c(170, 10), c(171, 10.5))), crs = crs))

# The distance on the sphere, in km, from each point to longitude 180 at
# the point's own latitude.  The wrap adds and takes away 360, which moves
# a longitude by at most a few ulps of 360 (about 1e-13 degrees, a few
# nanometres): 1 km is far wider than that on every platform, and far
# narrower than the 17,600 km the planar reading put the points at.
.am_km_from_180 <- function(pts) {
  xy <- sf::st_coordinates(pts)
  at <- sf::st_sfc(lapply(xy[, 2L], function(y) sf::st_point(c(180, y))), crs = 4326)
  p  <- sf::st_set_crs(sf::st_geometry(pts), 4326)
  old <- suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old)))
  as.numeric(sf::st_distance(p, at, by_element = TRUE)) / 1000
}

.am_catch <- function(expr) {
  ws <- list()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws, texts = vapply(ws, conditionMessage, character(1)))
}

# The box centres and planar centroids of the rows that do not cross 180, as
# the planar reading takes them (the points before the change).
.am_planar_box <- function(g)
  t(vapply(g, function(s) {
    bb <- sf::st_bbox(s)
    c((bb[["xmin"]] + bb[["xmax"]]) / 2, (bb[["ymin"]] + bb[["ymax"]]) / 2)
  }, numeric(2)))

test_that("bbox_center puts a geometry across longitude 180 at the centre of its box read across 180", {
  x <- .am_layer()
  for (s2 in c(TRUE, FALSE)) {
    old <- suppressMessages(sf::sf_use_s2(s2))
    on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
    r <- .am_catch(coerce_to_points(x, "bbox_center"))
    expect_length(r$warnings, 0L)
    expect_true(all(.am_km_from_180(r$value[1:3, ]) < 1), info = paste("s2", s2))
    xy <- unname(sf::st_coordinates(r$value))
    expect_identical(xy[1:3, 2L], c(10.5, 10.25, 10.25))
    # The rows that do not cross 180 keep their planar box centre, to the bit.
    expect_identical(xy[4:5, ], unname(.am_planar_box(sf::st_geometry(x)[4:5])))
    expect_identical(xy[4:5, ], rbind(c(170.5, 10.5), c(170.5, 10.25)))
  }
})

test_that("centroid with s2 off takes GEOS's centroid read across 180; with s2 on s2's is kept", {
  x <- .am_layer()
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  r <- .am_catch(coerce_to_points(x, "centroid"))
  expect_length(r$warnings, 0L)
  expect_true(all(.am_km_from_180(r$value[1:3, ]) < 1))
  # The other rows keep GEOS's planar centroid.
  planar <- suppressWarnings(sf::st_centroid(sf::st_geometry(x)[4:5]))
  expect_identical(sf::st_geometry(r$value)[4:5], planar)
  # The centroid "line_midpoint" takes of a row that is not a LINESTRING.
  lm <- .am_catch(coerce_to_points(x[c(1, 3, 4), ], "line_midpoint"))
  expect_true(all(.am_km_from_180(lm$value[1:2, ]) < 1))
  # With s2 on, s2's centroids on the sphere are returned as before.
  suppressMessages(sf::sf_use_s2(TRUE))
  r_on <- .am_catch(coerce_to_points(x, "centroid"))
  expect_identical(sf::st_geometry(r_on$value), sf::st_centroid(sf::st_geometry(x)))
})

test_that("a layer without a CRS read as degrees gets the same points, with the plain warning", {
  x  <- .am_layer()
  xn <- .am_layer(NA_integer_)
  want <- sf::st_coordinates(coerce_to_points(x, "bbox_center"))
  r <- .am_catch(coerce_to_points(xn, "bbox_center"))
  expect_true(is.na(sf::st_crs(r$value)))
  expect_identical(sf::st_coordinates(r$value), want)
  expect_identical(r$texts, paste0(
    "coerce_to_points(): `x` has no CRS and its coordinates look like lon/lat, so rows 1, 2 ",
    "and 3 of `x`, which span more than 180 degrees of longitude, are read as EPSG:4326 and ",
    "given the centres of their boxes read across longitude 180. Set the CRS of `x` ",
    "explicitly to say how to read it."))
  rc <- .am_catch(coerce_to_points(xn[1, ], "centroid"))
  expect_true(.am_km_from_180(rc$value) < 1)
  expect_identical(rc$texts, paste0(
    "coerce_to_points(): `x` has no CRS and its coordinates look like lon/lat, so row 1 of ",
    "`x`, which spans more than 180 degrees of longitude, is read as EPSG:4326 and given its ",
    "centroid read across longitude 180. Set the CRS of `x` explicitly to say how to read it."))
  # A planar layer without a CRS keeps its planar points, with nothing said.
  pl <- sf::st_sf(id = 1, geometry = sf::st_sfc(.am_sq(-5e5, 5e5, 4e6, 4.1e6)))
  r <- .am_catch(coerce_to_points(pl, "bbox_center"))
  expect_length(r$warnings, 0L)
  expect_identical(unname(sf::st_coordinates(r$value)), cbind(0, 4.05e6))
})

test_that("a point read across 180 east or west of 180 is brought back into -180 to 180", {
  # A square from 179 to -177 is read from 179 to 183, centre 181, and comes
  # back at -179; one from 178 to -179 is read from 178 to 181, centre 179.5.
  # Each centre is a sum and a halving of whole and half degrees, and the
  # wrap takes away 360 from a whole number, so the box centres are exact.
  x <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(.am_sq(179, -177), .am_sq(178, -179), crs = 4326))
  want <- rbind(c(-179, 10.5), c(179.5, 10.5))
  for (s2 in c(TRUE, FALSE)) {
    old <- suppressMessages(sf::sf_use_s2(s2))
    on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
    r <- .am_catch(coerce_to_points(x, "bbox_center"))
    expect_length(r$warnings, 0L)
    expect_identical(unname(sf::st_coordinates(r$value)), want, info = paste("s2", s2))
  }
  # GEOS's centroid of the square read so is its centre up to rounding in
  # the area-weighted sum (a few ulps of 181); 1e-9 degrees, about 0.1 mm,
  # is far wider than that on every platform.
  suppressMessages(sf::sf_use_s2(FALSE))
  r <- .am_catch(coerce_to_points(x, "centroid"))
  expect_length(r$warnings, 0L)
  xy <- unname(sf::st_coordinates(r$value))
  expect_true(all(abs(xy - want) < 1e-9))
  expect_true(all(xy[, 1L] >= -180 & xy[, 1L] <= 180))
})

test_that("the reading across 180 takes a span of more than 180 as written and no more than 180 so read", {
  # A box from -90 to 90 spans exactly 180 degrees as written: it is not
  # read across 180, and keeps (0, y), not (180, y) half the globe away.
  b <- sf::st_sf(id = 1, geometry = sf::st_sfc(.am_sq(-90, 90), crs = 4326))
  # A line through 90, 180 and -90 spans 270 as written and exactly 180
  # read from 0 to 360 (90 to 270): it is read across 180, centre 180.
  l <- sf::st_sf(id = 1, geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(90, 10), c(180, 10.5), c(-90, 11))), crs = 4326))
  for (s2 in c(TRUE, FALSE)) {
    old <- suppressMessages(sf::sf_use_s2(s2))
    on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
    info <- paste("s2", s2)
    expect_identical(unname(sf::st_coordinates(coerce_to_points(b, "bbox_center"))),
                     cbind(0, 10.5), info = info)
    expect_identical(unname(sf::st_coordinates(coerce_to_points(l, "bbox_center"))),
                     cbind(180, 10.5), info = info)
  }
})

test_that("\"auto\" with s2 off takes the centroid of a MULTIPOINT across 180 read across 180", {
  # GEOS's centroid of (179.5, 10) and (-179.5, 10.5) on the plane was
  # (0, 10.25), about 17,700 km from either point.
  x <- .am_layer()
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  r <- .am_catch(coerce_to_points(x[c(3, 4), ], "auto"))
  expect_length(r$warnings, 0L)
  expect_true(.am_km_from_180(r$value[1, ]) < 1)
  expect_identical(sf::st_coordinates(r$value)[1, "Y"], c(Y = 10.25))
  # The same layer without a CRS, read as degrees, with the plain warning.
  rn <- .am_catch(coerce_to_points(.am_layer(NA_integer_)[3, ], "auto"))
  expect_true(.am_km_from_180(rn$value) < 1)
  expect_identical(rn$texts, paste0(
    "coerce_to_points(): `x` has no CRS and its coordinates look like lon/lat, so row 1 of ",
    "`x`, which spans more than 180 degrees of longitude, is read as EPSG:4326 and given its ",
    "centroid read across longitude 180. Set the CRS of `x` explicitly to say how to read it."))
  # A MULTIPOINT that does not cross 180 keeps GEOS's planar centroid.
  mp <- sf::st_sf(id = 1, geometry = sf::st_sfc(
    sf::st_multipoint(rbind(c(170, 10), c(171, 10.5))), crs = 4326))
  expect_identical(sf::st_geometry(coerce_to_points(mp, "auto")),
                   suppressWarnings(sf::st_centroid(sf::st_geometry(mp))))
})

test_that("a geometry wider than 180 degrees in both readings keeps its planar point", {
  # A ring round the North Pole at latitude 80.
  lon <- c(-170, -90, 0, 90, 170, -170)
  ring <- sf::st_sf(id = 1, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(cbind(lon, 80), c(-170, 80)))), crs = 4326))
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  r <- .am_catch(coerce_to_points(ring, "bbox_center"))
  expect_length(r$warnings, 0L)
  expect_identical(unname(sf::st_coordinates(r$value)),
                   unname(.am_planar_box(sf::st_geometry(ring))))
})

test_that("prep_model_data(pointize = \"bbox_center\") projects squares across 180 as \"auto\" does", {
  d <- sf::st_sf(resp = 1:3, x1 = c(1, 3, 2), geometry = sf::st_sfc(
    .am_sq(179.5, -179.5), .am_sq(178, 179), .am_sq(-178, -177), crs = 4326))
  pb <- suppressWarnings(prep_model_data(d, "resp", "x1", pointize = "bbox_center"))
  pa <- suppressWarnings(prep_model_data(d, "resp", "x1", pointize = "auto"))
  expect_false(grepl("+lat_0=90", sf::st_crs(pb)$proj4string, fixed = TRUE))
  expect_identical(sf::st_crs(pb), sf::st_crs(pa))
  # Every point lies within the layer's 4 degrees of longitude of the centre.
  expect_true(all(abs(sf::st_coordinates(pb)[, 1L]) < 3e5))
})
