# tests/testthat/test-crs-selection.R
# ---------------------------------------------------------------------------
# Specification for .pick_local_projected_crs() / ensure_projected().
#
# The narrow-extent cases are regression guards: the common case must keep
# getting a UTM zone.  The wide-extent cases pin the equal-area branch that
# replaced forcing continental data into a single zone.
# ---------------------------------------------------------------------------

ll_points <- function(lon, lat, n = NULL, seed = 1) {
  if (is.null(n)) {
    df <- data.frame(lon = lon, lat = lat)
  } else {
    set.seed(seed)
    df <- data.frame(lon = runif(n, lon[1], lon[2]),
                     lat = runif(n, lat[1], lat[2]))
  }
  sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
}

proj4_of <- function(x) sf::st_crs(x)$proj4string


# ---- Narrow extents: UTM, unchanged ---------------------------------------

test_that("narrow northern-hemisphere extents get the right UTM zone", {
  # Wake County, NC.  Centroid ~ -78.6 lon -> zone 17 -> EPSG 32617.
  pts <- ll_points(c(-78.7, -78.5), c(35.7, 35.9))
  expect_equal(sf::st_crs(ensure_projected(pts))$epsg, 32617L)
})

test_that("narrow southern-hemisphere extents get a 327xx zone", {
  pts <- ll_points(c(151.1, 151.3), c(-33.9, -33.7))
  epsg <- sf::st_crs(ensure_projected(pts))$epsg
  expect_true(epsg >= 32701L && epsg <= 32760L)
})

test_that("narrow extents emit no wide-extent log line", {
  pts <- ll_points(c(-78.7, -78.5), c(35.7, 35.9))
  lines <- capture_spatialkit_log(ensure_projected(pts))
  expect_false(log_has(lines, "central meridian of UTM zone"))
})

test_that("already-projected input is left alone", {
  pts <- sf::st_as_sf(data.frame(x = c(1000, 2000), y = c(3000, 4000)),
                      coords = c("x", "y"), crs = 3857)
  expect_equal(sf::st_crs(ensure_projected(pts))$epsg, 3857L)
})


# ---- Wide extents: equal-area ---------------------------------------------

test_that("continental east-west extents get Albers equal-area", {
  # CONUS spans ~57 deg of longitude.  A single UTM zone is 6 deg wide, so the
  # extent edge would sit ~29 deg from the central meridian -- roughly 2470 km
  # at 40N, where k = k0 * (1 + x^2 / 2R^2) implies about +7.5% distance error.
  pts <- ll_points(c(-124, -67), c(25, 49))
  out  <- ensure_projected(pts)
  expect_match(proj4_of(out), "\\+proj=aea", fixed = FALSE)
  expect_false(grepl("utm", proj4_of(out), fixed = TRUE))
})

test_that("wide extents log the projection switch", {
  # .log_warn() routes through logger, which raises no R condition -- see
  # helper-logging.R.  expect_warning() here would pass vacuously.
  pts   <- ll_points(c(-124, -67), c(25, 49))
  lines <- capture_spatialkit_log(ensure_projected(pts))
  expect_true(log_has(lines, "central meridian of UTM zone"))
  expect_true(log_has(lines, "Albers equal-area"))
})

test_that("equatorial wide extents get LAEA, not Albers", {
  # abs(lat) <= 20 -- Albers standard parallels degenerate near the equator.
  pts <- ll_points(c(-20, 20), c(-8, 8))
  out <- ensure_projected(pts)
  expect_match(proj4_of(out), "\\+proj=laea", fixed = FALSE)
})

test_that("tall narrow north-south extents KEEP their UTM zone", {
  # Transverse Mercator error scales with distance from the central meridian,
  # and cos(lat) shrinks that distance -- a 4-deg-wide, 60-deg-tall strip peaks
  # around +0.02%.  This is UTM's design case.  Switching it to LAEA would be
  # roughly 50x worse (about -1.1% at 30 deg from the projection centre), so
  # latitude span must NOT trigger the equal-area branch.
  pts  <- ll_points(c(-72, -68), c(-40, 20))
  out  <- ensure_projected(pts)
  epsg <- sf::st_crs(out)$epsg
  expect_true(!is.na(epsg) && ((epsg >= 32601L && epsg <= 32660L) ||
                               (epsg >= 32701L && epsg <= 32760L)))

  lines <- capture_spatialkit_log(ensure_projected(pts))
  expect_false(log_has(lines, "central meridian"))
})

test_that("the switch is driven by offset from the zone central meridian", {
  # narrow: centroid -97.8 -> zone 14 (CM -99) -> max offset 3.5 deg -> UTM
  # wide:   centroid -96.0 -> zone 15 (CM -93) -> max offset 7.0 deg -> equal-area
  narrow <- ll_points(c(-100.0, -95.5), c(38, 42))
  wide   <- ll_points(c(-100.0, -92.0), c(38, 42))
  expect_true(!is.na(sf::st_crs(ensure_projected(narrow))$epsg))
  expect_match(proj4_of(ensure_projected(wide)),
               "\\+proj=(aea|laea)", fixed = FALSE)
})


# ---- Downstream propagation ------------------------------------------------

test_that("the projection choice changes downstream distance estimates", {
  # The reason this matters is the consumers, not the projection itself.
  # A continental extent forced into one UTM zone yields a materially
  # different autocorrelation range and block grid than an equal-area CRS.
  skip_if_not_installed("gstat")

  pts <- ll_points(c(-124, -67), c(25, 49), n = 400, seed = 3)
  set.seed(11)
  pts$z <- stats::rnorm(nrow(pts))

  auto <- ensure_projected(pts)                  # equal-area
  utm  <- ensure_projected(pts, target_crs = sf::st_crs(32615))    # forced zone

  bb_auto <- sf::st_bbox(auto)
  bb_utm  <- sf::st_bbox(utm)
  w_auto  <- as.numeric(bb_auto[["xmax"]] - bb_auto[["xmin"]])
  w_utm   <- as.numeric(bb_utm[["xmax"]]  - bb_utm[["xmin"]])
  h_auto  <- as.numeric(bb_auto[["ymax"]] - bb_auto[["ymin"]])
  h_utm   <- as.numeric(bb_utm[["ymax"]]  - bb_utm[["ymin"]])

  # Forcing one zone stretches east-west distances by several percent.
  expect_gt(abs(w_utm - w_auto) / w_auto, 0.02)
  expect_gt(abs(h_utm - h_auto) / h_auto, 0.02)

  # And that difference reaches block sizing, which works in CRS units.
  # `block_size` is a minimum block EDGE LENGTH, so make_folds() derives
  # nx = floor(w / block_size).  Sizing one block at exactly 1/50th of the
  # equal-area width therefore gives nx = 50 there by construction, and
  # floor(50 * w_utm / w_auto) in the forced zone -- which the >2% assertion
  # above guarantees is not 50.  Same argument for ny.
  #
  # The DEFAULT sizing would not show this: it is purely geometric (nx/ny come
  # from the aspect ratio and k), so it lands on the same grid in both CRSs
  # while each block covers a different distance on the ground.
  block_size <- w_auto / 50
  f_auto <- make_folds(auto, k = 4, method = "block_kfold", seed = 1,
                       block_size = block_size)
  f_utm  <- make_folds(utm,  k = 4, method = "block_kfold", seed = 1,
                       block_size = block_size)

  expect_equal(f_auto$params$grid_nx, 50L)
  expect_equal(f_utm$params$grid_nx, as.integer(floor(50 * w_utm / w_auto)))
  expect_false(f_utm$params$grid_nx == f_auto$params$grid_nx)
  expect_false(f_utm$params$grid_ny == f_auto$params$grid_ny)

  # The geometric default, by contrast, does NOT depend on the CRS -- worth
  # pinning so the assertions above are read as being about `block_size`.
  g_auto <- make_folds(auto, k = 4, method = "block_kfold", seed = 1)
  g_utm  <- make_folds(utm,  k = 4, method = "block_kfold", seed = 1)
  expect_equal(g_auto$params$grid_nx, g_utm$params$grid_nx)
  expect_equal(g_auto$params$grid_ny, g_utm$params$grid_ny)
})


test_that("antimeridian data are recognised, not flattened onto EPSG:3857", {
  # A bounding BOX cannot distinguish global coverage from data straddling the
  # antimeridian -- but the coordinates can: a wrapped layer leaves one very
  # large gap in the sorted longitudes.  This matters because EPSG:3857 SPLITS
  # such a layer: the two points below are 2 deg apart and Web Mercator puts
  # them on opposite sides of the world.
  pts <- ll_points(c(-179, 179), c(-10, 10))
  out <- ensure_projected(pts)

  lines <- capture_spatialkit_log(ensure_projected(pts))
  expect_true(log_has(lines, "straddle the antimeridian"))
  expect_false(isTRUE(sf::st_crs(out)$epsg == 3857L))

  # The test that matters: distance survives the projection.
  d_true <- as.numeric(sf::st_distance(pts)[1, 2])
  d_proj <- as.numeric(stats::dist(sf::st_coordinates(out)))
  expect_lt(abs(d_proj / d_true - 1), 0.05)

  # Genuinely global coverage still falls back, because no local projection
  # fits it -- and it says so.
  glob <- ll_points(seq(-170, 170, by = 20), rep(c(-40, 0, 40), length.out = 18))
  gout <- ensure_projected(glob)
  expect_equal(sf::st_crs(gout)$epsg, 3857L)
  glines <- capture_spatialkit_log(ensure_projected(glob))
  expect_true(log_has(glines, "global coverage"))
})
