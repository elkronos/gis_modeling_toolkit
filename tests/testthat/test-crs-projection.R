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


test_that("harmonize_crs aligns two objects", {
  a <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326))
  b <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 3857))
  result <- harmonize_crs(a, b, prefer = "b")
  expect_equal(sf::st_crs(result$a), sf::st_crs(b))
})
