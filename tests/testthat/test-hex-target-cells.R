# ---------------------------------------------------------------------------
# Regression tests: hex-grid target_cells packing correction
#
# Guards against the inverted packing adjustment in create_grid_polygons(),
# which divided target_cells by sqrt(3)/2 instead of multiplying, so hex
# grids overshot the requested cell count by ~33%.
# ---------------------------------------------------------------------------

.square_boundary <- function(side = 100) {
  sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(side, 0), c(side, side), c(0, side), c(0, 0)
      ))),
      crs = 32632
    )
  )
}

test_that("hex grid cell count lands near target_cells", {
  bnd <- .square_boundary(100)
  target <- 400

  hex <- create_grid_polygons(bnd, target_cells = target, type = "hex")
  n_hex <- nrow(hex)

  # Approximate target: edge partials inflate the count somewhat, but the
  # inverted correction produced ~1.35x the target. Accept up to 1.25x.
  expect_gt(n_hex, 0.5 * target)
  expect_lt(n_hex, 1.25 * target)
})

test_that("hex and square grids give comparable counts for the same target", {
  bnd <- .square_boundary(100)
  target <- 400

  sq  <- create_grid_polygons(bnd, target_cells = target, type = "square")
  hex <- create_grid_polygons(bnd, target_cells = target, type = "hex")

  # With the inverted correction, n_hex / n_square was ~1.35-1.4.
  # After the fix the two grid types should land in the same ballpark.
  ratio <- nrow(hex) / nrow(sq)
  expect_gt(ratio, 0.6)
  expect_lt(ratio, 1.25)
})
