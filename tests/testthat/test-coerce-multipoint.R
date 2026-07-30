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
