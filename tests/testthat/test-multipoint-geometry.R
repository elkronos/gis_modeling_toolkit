# ===========================================================================
# MULTIPOINT handling across the package.
#
# The recurring hazard is that sf::st_coordinates() returns one row per
# VERTEX, not per feature, so any code that pairs it with attribute rows
# silently misaligns the moment a MULTIPOINT appears. These tests pin the
# feature-count contract at every entry point that accepts geometry.
# ===========================================================================

test_that(".dedup_points removes duplicate POINT features (fast path)", {
  g <- sf::st_sfc(
    sf::st_point(c(0, 0)),
    sf::st_point(c(1, 1)),
    sf::st_point(c(0, 0)),  # duplicate of feature 1
    crs = 32632
  )
  out <- spatialkit:::.dedup_points(g)
  expect_length(out, 2L)
})


test_that(".dedup_points handles MULTIPOINT without misaligning features", {
  mp  <- sf::st_multipoint(rbind(c(0, 0), c(1, 1), c(2, 2)))
  g <- sf::st_sfc(
    mp,                          # feature 1
    mp,                          # feature 2: duplicate of feature 1
    sf::st_point(c(5, 5)),       # feature 3
    sf::st_point(c(5, 5)),       # feature 4: duplicate of feature 3
    sf::st_point(c(9, 9)),       # feature 5: unique
    crs = 32632
  )
  # Before the fix, st_coordinates() returned one row per *vertex* (7 rows
  # for 5 features), so the dedup mask was misaligned with features.
  out <- spatialkit:::.dedup_points(g)
  expect_length(out, 3L)
  gtypes <- as.character(sf::st_geometry_type(out, by_geometry = TRUE))
  expect_equal(sort(gtypes), c("MULTIPOINT", "POINT", "POINT"))

  # sf-object variant must subset rows, not vertices
  sf_obj <- sf::st_sf(id = 1:5, geometry = g)
  out_sf <- spatialkit:::.dedup_points(sf_obj)
  expect_s3_class(out_sf, "sf")
  expect_equal(nrow(out_sf), 3L)
  expect_equal(out_sf$id, c(1L, 3L, 5L))
})


test_that("create_voronoi_polygons accepts MULTIPOINT input", {
  g <- sf::st_sfc(
    sf::st_multipoint(rbind(c(0, 0), c(10, 0))),
    sf::st_multipoint(rbind(c(0, 10), c(10, 10))),
    sf::st_point(c(5, 5)),
    crs = 32632
  )
  pts <- sf::st_sf(geometry = g)
  res <- suppressMessages(create_voronoi_polygons(pts, quiet = TRUE))
  expect_true(is.list(res))
  expect_s3_class(res$cells, "sf")
  expect_gt(nrow(res$cells), 0L)
})


test_that("prep_model_data coerces multi-vertex MULTIPOINT to POINT", {
  skip_if_not_installed("sf")

  g <- sf::st_sfc(
    sf::st_multipoint(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10))),
    sf::st_multipoint(rbind(c(100, 100), c(120, 100))),
    sf::st_point(c(50, 50)),
    crs = 32632
  )
  dat <- sf::st_sf(y = c(1, 2, 3), x1 = c(10, 20, 30), geometry = g)

  out <- prep_model_data(dat, "y", "x1")

  # One row per feature, all plain POINT (no per-vertex explosion)
  expect_equal(nrow(out), 3L)
  gtypes <- as.character(sf::st_geometry_type(out, by_geometry = TRUE))
  expect_true(all(gtypes == "POINT"))

  # MULTIPOINT features collapse to their centroid
  xy <- sf::st_coordinates(out)
  expect_equal(unname(xy[1, ]), c(5, 5))
  expect_equal(unname(xy[2, ]), c(110, 100))
  expect_equal(unname(xy[3, ]), c(50, 50))

  # st_coordinates() now aligns 1:1 with attribute rows
  expect_equal(nrow(xy), nrow(out))
})


test_that("fit_bayesian_spatial_model rejects MULTIPOINT when .already_prepped", {
  skip_if_not_installed("sf")
  skip_if_not_installed("brms")

  g <- sf::st_sfc(
    sf::st_multipoint(rbind(c(0, 0), c(10, 0))),
    sf::st_point(c(5, 5)),
    sf::st_point(c(2, 8)),
    crs = 32632
  )
  dat <- sf::st_sf(y = c(1, 2, 3), x1 = c(10, 20, 30), geometry = g)

  expect_error(
    fit_bayesian_spatial_model(dat, "y", "x1", .already_prepped = TRUE),
    "must be POINT"
  )
})


test_that("coerce_to_points() uses centroid for MULTIPOINT, not first sub-point", {
  skip_if_not_installed("sf")

  # Build a MULTIPOINT whose centroid differs from its first sub-point
  mp <- sf::st_multipoint(matrix(c(0, 0,
                                    10, 0,
                                    10, 10,
                                    0, 10), ncol = 2, byrow = TRUE))
  sf_obj <- sf::st_sf(geometry = sf::st_sfc(mp, crs = 4326))

  result <- coerce_to_points(sf_obj, mode = "auto")

  res_coords <- sf::st_coordinates(result)

  # Centroid of the four corners is (5, 5), NOT the first point (0, 0).
  # st_centroid on EPSG:4326 computes a geodesic centroid, so we use
  # unname() to drop the column-name attribute and a tolerance that
  # accommodates the slight spherical deviation from the arithmetic mean.
  expect_equal(unname(res_coords[1, "X"]), 5, tolerance = 0.1)
  expect_equal(unname(res_coords[1, "Y"]), 5, tolerance = 0.1)
})


test_that("coerce_to_points() returns POINT geometry for MULTIPOINT input", {
  skip_if_not_installed("sf")

  mp <- sf::st_multipoint(matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE))
  sf_obj <- sf::st_sf(geometry = sf::st_sfc(mp, crs = 4326))

  result <- coerce_to_points(sf_obj, mode = "auto")

  gtypes <- as.character(sf::st_geometry_type(result))
  expect_true(all(gtypes == "POINT"))
})
