# ===========================================================================
# Small geometry utilities: input assertions, target clipping, grid cells.
# ===========================================================================

test_that(".assert_sf rejects non-sf objects", {
  expect_error(spatialkit:::.assert_sf(data.frame(x = 1)), "Expected an sf object")
})


test_that("clip_target_for works with points only", {
  pts <- sf::st_sf(
    id = 1:5,
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)), sf::st_point(c(1, 1)),
      sf::st_point(c(0.5, 0.5)),
      crs = 32632
    )
  )
  tgt <- clip_target_for(pts)
  expect_s3_class(tgt, "sf")
  expect_true(all(sf::st_geometry_type(tgt) %in% c("POLYGON", "MULTIPOLYGON")))
})


test_that("create_grid_polygons produces cells", {
  bnd <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))),
      crs = 32632
    )
  )
  grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  expect_s3_class(grid, "sf")
  expect_true(nrow(grid) >= 4)
})
