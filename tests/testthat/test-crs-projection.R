# ===========================================================================
# CRS detection, projection, and harmonisation.
# Choosing a projection is upstream of everything the package measures:
# variogram ranges, block sizes, GWR bandwidths and GP length-scales are
# all in CRS units. See also test-crs-selection.R for the wide-extent rules.
# ===========================================================================

test_that(".is_longlat detects geographic CRS", {
  pts <- sf::st_sfc(sf::st_point(c(-122.4, 37.8)), crs = 4326)
  expect_true(spatialkit:::.is_longlat(pts))

  pts_proj <- sf::st_transform(pts, 32610)
  expect_false(spatialkit:::.is_longlat(pts_proj))
})


test_that("ensure_projected transforms geographic data", {
  pts <- sf::st_sf(
    id = 1:3,
    geometry = sf::st_sfc(
      sf::st_point(c(-122.4, 37.8)),
      sf::st_point(c(-122.3, 37.7)),
      sf::st_point(c(-122.5, 37.9)),
      crs = 4326
    )
  )
  result <- ensure_projected(pts)
  expect_false(sf::st_is_longlat(result))
})


test_that("ensure_projected does NOT assume 4326 for small integer-like coords", {
  # Simulates a local site survey in metres with coords outside the

  # geographic envelope (values > 180 / > 90) so the heuristic cannot
  # mistakenly interpret them as lon/lat.
  pts <- sf::st_sf(
    id = 1:3,
    geometry = sf::st_sfc(
      sf::st_point(c(500, 600)),
      sf::st_point(c(510, 620)),
      sf::st_point(c(505, 610))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # CRS should still be NA — the function must NOT guess 4326

  expect_true(is.na(sf::st_crs(result)))
})


test_that("ensure_projected assumes 4326 for CRS-less data with geographic precision", {
  # Realistic lon/lat coordinates without a CRS set
  pts <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(
      sf::st_point(c(-0.1278, 51.5074)),
      sf::st_point(c(-0.1385, 51.5013))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # Should have been assigned a projected CRS
  expect_false(is.na(sf::st_crs(result)))
  expect_false(sf::st_is_longlat(result))
})


test_that("ensure_projected assumes 4326 for CRS-less data with large extent", {
  pts <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(
      sf::st_point(c(10, 20)),
      sf::st_point(c(15, 25))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # Extent > 1 degree in both axes → should assume 4326
  expect_false(is.na(sf::st_crs(result)))
  expect_false(sf::st_is_longlat(result))
})


# ---------------------------------------------------------------------------
# harmonize_crs()
#
# The test point is (10 E, 50 N), NOT the origin.  (0, 0) is the origin of both
# EPSG:4326 and EPSG:3857, so its coordinates do not move under the transform
# -- an implementation that called st_set_crs() instead of st_transform() (the
# exact hazard R/crs-geometry.R warns about) would relabel the CRS and pass a
# test written on that point identically.
#
# The expected projected coordinates are computed from the Web Mercator
# forward equations rather than from sf, so the assertion does not check
# st_transform() against itself:
#     x = R * lambda,  y = R * ln(tan(pi/4 + phi/2)),  R = 6378137
# ---------------------------------------------------------------------------

.wm_R    <- 6378137
.wm_x    <- .wm_R * (10 * pi / 180)
.wm_y    <- .wm_R * log(tan(pi / 4 + (50 * pi / 180) / 2))
.ll_pt   <- function() sf::st_sf(id = 1L,
                                 geometry = sf::st_sfc(sf::st_point(c(10, 50)),
                                                       crs = 4326))
.wm_pt   <- function() sf::st_sf(id = 2L,
                                 geometry = sf::st_sfc(sf::st_point(c(.wm_x, .wm_y)),
                                                       crs = 3857))

test_that("harmonize_crs(prefer = 'b') reprojects a onto b's CRS", {
  a <- .ll_pt(); b <- .wm_pt()
  result <- harmonize_crs(a, b, prefer = "b")

  expect_equal(sf::st_crs(result$a), sf::st_crs(3857))
  expect_equal(sf::st_crs(result$b), sf::st_crs(3857))

  # `a` MOVED: metres, not degrees, matching the hand-computed projection.
  expect_equal(unname(sf::st_coordinates(result$a)[1, 1:2]),
               c(.wm_x, .wm_y), tolerance = 1e-6)
  # `b` was already in the target CRS and must be untouched.
  expect_equal(unname(sf::st_coordinates(result$b)[1, 1:2]),
               c(.wm_x, .wm_y), tolerance = 1e-9)
})

test_that("harmonize_crs(prefer = 'a') reprojects b onto a's CRS", {
  a <- .ll_pt(); b <- .wm_pt()
  result <- harmonize_crs(a, b, prefer = "a")

  expect_equal(sf::st_crs(result$a), sf::st_crs(4326))
  expect_equal(sf::st_crs(result$b), sf::st_crs(4326))
  expect_equal(unname(sf::st_coordinates(result$b)[1, 1:2]), c(10, 50),
               tolerance = 1e-7)
  expect_equal(unname(sf::st_coordinates(result$a)[1, 1:2]), c(10, 50),
               tolerance = 1e-9)
})

test_that("harmonize_crs(target_crs) moves BOTH inputs to a third CRS", {
  a <- .ll_pt(); b <- .wm_pt()
  result <- harmonize_crs(a, b, target_crs = 32632)   # UTM 32N

  expect_equal(sf::st_crs(result$a), sf::st_crs(32632))
  expect_equal(sf::st_crs(result$b), sf::st_crs(32632))
  # Both describe the same place, so they must land on the same coordinates.
  expect_equal(unname(sf::st_coordinates(result$a)[1, 1:2]),
               unname(sf::st_coordinates(result$b)[1, 1:2]), tolerance = 1e-3)
  # ... and that place is not either input's original coordinates.
  expect_gt(abs(sf::st_coordinates(result$a)[1, 1] - 10), 1000)
})

test_that("harmonize_crs stamps a CRS-less input and says so", {
  # st_set_crs() relabels without moving anything, so the assumption has to be
  # announced.  .log_warn() raises no R condition -- see helper-logging.R.
  b <- .wm_pt()
  a_na <- sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_point(c(10, 50))))

  lines <- capture_spatialkit_log(res <- harmonize_crs(a_na, b))
  expect_true(log_has(lines, "`a` has no CRS"))
  expect_true(log_has(lines, "WITHOUT reprojection"))
  expect_equal(sf::st_crs(res$a), sf::st_crs(3857))
  # Stamped, not transformed: the coordinates are exactly as supplied.
  expect_equal(unname(sf::st_coordinates(res$a)[1, 1:2]), c(10, 50))

  # The mirror branch, with the CRS-less object second.
  lines_b <- capture_spatialkit_log(res_b <- harmonize_crs(b, a_na))
  expect_true(log_has(lines_b, "`b` has no CRS"))
  expect_equal(sf::st_crs(res_b$b), sf::st_crs(3857))
  expect_equal(unname(sf::st_coordinates(res_b$b)[1, 1:2]), c(10, 50))

  # Neither input has a CRS: nothing to align, nothing to warn about.
  b_na <- sf::st_sf(id = 2L, geometry = sf::st_sfc(sf::st_point(c(1, 2))))
  lines_none <- capture_spatialkit_log(res_none <- harmonize_crs(a_na, b_na))
  expect_false(log_has(lines_none, "has no CRS"))
  expect_true(is.na(sf::st_crs(res_none$a)))
  expect_true(is.na(sf::st_crs(res_none$b)))
})

test_that("harmonize_crs on_transform_error controls a failed transform", {
  a <- .ll_pt(); b <- .wm_pt()
  # A syntactically parseable but geometrically impossible projection: PROJ
  # rejects |lat_0| > 90, so st_transform() fails rather than returning
  # something wrong.
  bad_crs <- paste("+proj=omerc +lat_0=91 +lonc=0 +alpha=0 +k=1",
                   "+x_0=0 +y_0=0 +datum=WGS84")

  # Default: surface the failure.
  expect_error(
    suppressWarnings(harmonize_crs(a, b, target_crs = bad_crs)),
    "st_transform\\(\\) failed"
  )

  # Opt-in override: no error, but the unsafe fallback is announced.
  lines <- capture_spatialkit_log(
    suppressWarnings(harmonize_crs(a, b, target_crs = bad_crs,
                                   on_transform_error = "set_crs"))
  )
  expect_true(log_has(lines, "coordinates are NOT reprojected"))
})

test_that("harmonize_crs rejects non-spatial input", {
  a <- .ll_pt()
  expect_error(harmonize_crs(data.frame(x = 1), a), "`a` must be sf or sfc")
  expect_error(harmonize_crs(a, data.frame(x = 1)), "`b` must be sf or sfc")
})
