# tests/testthat/test-fix-crs-lat-beyond-90.R
# ---------------------------------------------------------------------------
# BUG-205.  A row of a lon/lat layer with a vertex at a latitude beyond 90
# has no place on the Earth, and no CRS places it, but the CRS chooser read
# it: it shaped the box, the wrap test's longitudes, the centre and the
# Albers parallels of the projection chosen for the other rows.  39 stations
# over the conterminous United States get an Albers CRS (lat_1 32.9, lat_2
# 43.0); with a 40th row at (-96, 95) they got a Lambert azimuthal CRS
# centred on (-96.7, 40.7), and the Albers candidate's parallels moved to
# 41.2 and 84.2.
#
# Such a row now takes no part in the choice, the whole row, as a row with a
# coordinate beyond 1e150: the choice is the one made without it.  The layer
# returned is not changed: ensure_projected() moves the row, which comes
# back empty, and its spatialkit_rows_unplaced warning names it.
# ---------------------------------------------------------------------------

.clb_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
.clb_st39 <- function() {
  set.seed(1)
  sf::st_as_sf(data.frame(x = stats::runif(39, -120, -75), y = stats::runif(39, 30, 47)),
               coords = c("x", "y"), crs = 4326)
}
.clb_with <- function(st, geom) {
  g <- c(sf::st_geometry(st), sf::st_sfc(geom, crs = 4326))
  sf::st_sf(geometry = g)
}
# ensure_projected(), with the spatialkit_rows_unplaced warning kept.
.clb_project <- function(x, purpose) {
  w <- NULL
  e <- withCallingHandlers(suppressMessages(ensure_projected(x, purpose = purpose)),
    warning = function(cnd) {
      if (inherits(cnd, "spatialkit_rows_unplaced")) w <<- cnd
      invokeRestart("muffleWarning")
    })
  list(e = e, w = w)
}

test_that(".beyond_pole_features() finds the rows with a vertex at a latitude beyond 90", {
  g <- sf::st_sfc(
    sf::st_point(c(-96, 40)),
    sf::st_point(c(-96, 95)),
    sf::st_point(c(-96, -90)),
    sf::st_point(),
    sf::st_linestring(rbind(c(0, 0), c(1, 91))),
    sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0)))),
    sf::st_multipolygon(list(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0))),
                             list(rbind(c(5, 89), c(6, 89), c(6, 95), c(5, 89))))),
    sf::st_geometrycollection(list(sf::st_point(c(1, 1)), sf::st_point(c(2, -100)))),
    sf::st_point(c(Inf, 40)),
    sf::st_point(c(95, 40)),
    sf::st_point(c(0, 90.000001)),
    sf::st_point(c(0, -90.000001)),
    sf::st_point(c(0, 90)),
    crs = 4326)
  expect_identical(.beyond_pole_features(g),
                   c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE,
                     TRUE, TRUE, FALSE))
  # XYZ points read the latitude, not the height.
  z <- sf::st_sfc(sf::st_point(c(-96, 40, 500)), sf::st_point(c(-96, 95, 0)), crs = 4326)
  expect_identical(.beyond_pole_features(z), c(FALSE, TRUE))
  expect_identical(.beyond_pole_features(sf::st_sfc(crs = 4326)), logical(0))
})

test_that("a point at a latitude beyond 90 takes no part in the choice", {
  st39 <- .clb_st39()
  for (p in c("distance", "area")) {
    ref <- .clb_quiet(.pick_local_projected_crs(st39, purpose = p))
    expect_match(ref$crs$input, "+proj=aea", fixed = TRUE)
    ref_e <- .clb_project(st39, p)
    expect_null(ref_e$w)
    for (bad in list(c(-96, 95), c(-100, -95), c(40, -96))) {
      tag <- sprintf("(%g, %g), %s", bad[1L], bad[2L], p)
      x  <- .clb_with(st39, sf::st_point(bad))
      tr <- .clb_quiet(.pick_local_projected_crs(x, purpose = p))
      expect_identical(tr$crs, ref$crs, info = tag)
      expect_identical(tr$candidates, ref$candidates, info = tag)
      got <- .clb_project(x, p)
      expect_identical(sf::st_crs(got$e), sf::st_crs(ref_e$e), info = tag)
      expect_identical(attr(got$e, "crs_choice"), attr(ref_e$e, "crs_choice"), info = tag)
      expect_identical(sf::st_geometry(got$e)[1:39], sf::st_geometry(ref_e$e), info = tag)
      expect_identical(which(sf::st_is_empty(got$e)), 40L, info = tag)
      expect_identical(got$w$rows, 40L, info = tag)
    }
  }
})

test_that("a polygon with one vertex at a latitude beyond 90 takes no part in the choice", {
  st39 <- .clb_st39()
  sq   <- sf::st_polygon(list(rbind(c(-96, 40), c(-95, 40), c(-95, 95), c(-96, 41),
                                    c(-96, 40))))
  x     <- .clb_with(st39, sq)
  as_it <- .clb_with(st39, sf::st_polygon())
  for (p in c("distance", "area")) {
    ref <- .clb_quiet(.pick_local_projected_crs(st39, purpose = p))
    tr  <- .clb_quiet(.pick_local_projected_crs(x, purpose = p))
    expect_identical(tr, .clb_quiet(.pick_local_projected_crs(as_it, purpose = p)), info = p)
    expect_identical(tr$crs, ref$crs, info = p)
    expect_identical(tr$candidates, ref$candidates, info = p)
    # The polygon is moved as before: it comes back without the vertex the
    # CRS cannot place, and the warning names it.
    got <- .clb_project(x, p)
    expect_identical(sf::st_crs(got$e), ref$crs, info = p)
    expect_identical(sf::st_geometry(got$e)[1:39],
                     sf::st_geometry(.clb_project(st39, p)$e), info = p)
    expect_identical(nrow(sf::st_coordinates(sf::st_geometry(got$e)[40L])), 4L, info = p)
    expect_identical(got$w$rows, 40L, info = p)
  }
})

test_that("a lon/lat layer whose every row has a latitude beyond 90 falls back and says why", {
  x <- sf::st_as_sf(data.frame(x = c(-96, 10), y = c(95, -100)),
                    coords = c("x", "y"), crs = 4326)
  lines <- capture_spatialkit_log(tr <- .clb_quiet(.pick_local_projected_crs(x)))
  expect_identical(tr$crs, sf::st_crs(3857))
  expect_null(tr$candidates)
  expect_true(log_has(lines, "no row has a location.*a latitude beyond 90.*EPSG:3857"))
  lines <- capture_spatialkit_log(
    tr <- .clb_quiet(.pick_local_projected_crs(x, purpose = "area")))
  expect_match(tr$crs$input, "+proj=eqearth", fixed = TRUE)
  expect_true(log_has(lines, "no row has a location.*a latitude beyond 90.*Equal Earth"))
})

test_that("the latitude is read in degrees on the WGS84 copy, whatever the geographic CRS", {
  # NAD83 in degrees: the same choice as for the layer in EPSG:4326.
  st39 <- .clb_st39()
  x    <- .clb_with(st39, sf::st_point(c(-96, 95)))
  x83  <- sf::st_set_crs(x, NA)
  x83  <- sf::st_set_crs(x83, 4269)
  st83 <- sf::st_set_crs(sf::st_set_crs(st39, NA), 4269)
  expect_identical(.clb_quiet(.pick_local_projected_crs(x83)),
                   .clb_quiet(.pick_local_projected_crs(st83)))
  # EPSG:4807 is in grads.  A row at 95 grads lies at 85.5 degrees: it has a
  # location, and the choice reads it, as on the layer moved to WGS84.  A row
  # at 101 grads, beyond the pole, takes no part.
  st07 <- sf::st_set_crs(sf::st_set_crs(st39, NA), 4807)
  alone <- .clb_quiet(.pick_local_projected_crs(st07))
  w07 <- function(lat) {
    sf::st_sf(geometry = c(sf::st_geometry(st07),
                           sf::st_sfc(sf::st_point(c(-96, lat)), crs = 4807)))
  }
  x95  <- w07(95)
  tr95 <- .clb_quiet(.pick_local_projected_crs(x95))
  expect_identical(tr95, .clb_quiet(.pick_local_projected_crs(sf::st_transform(x95, 4326))))
  expect_false(identical(tr95$crs, alone$crs))
  expect_identical(.clb_quiet(.pick_local_projected_crs(w07(101))), alone)
})
