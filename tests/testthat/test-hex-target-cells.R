# ---------------------------------------------------------------------------
# Regression tests: hex-grid target_cells packing correction
#
# Guards against the inverted packing adjustment in create_grid_polygons(),
# which divided target_cells by sqrt(3)/2 instead of multiplying, so hex
# grids overshot the requested cell count by ~33%.
#
# The bounds below are two-sided on purpose.  A one-sided lower bound of
# 0.5 * target admitted the mirror-image defect -- applying sqrt(3)/2 twice
# in the *correct* direction, which undershoots to ~0.75x -- so a test that
# only caught the overshoot would have stayed green through it.
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

  hex <- create_grid_polygons(bnd, target_cells = target, type = "hex",
                              quiet = TRUE)
  n_hex <- nrow(hex)

  # Clipping keeps every hexagon that merely overhangs the boundary, so the
  # count runs slightly ABOVE target; nx/ny rounding moves it a few percent
  # either way.  Neither effect approaches 15%, so the band is symmetric.
  expect_gt(n_hex, 0.85 * target)
  expect_lt(n_hex, 1.25 * target)
})

test_that("a hexagon's area equals boundary area / target_cells", {
  # The derivation the packing correction exists for.  st_make_grid() reads
  # `cellsize` as the flat-to-flat distance, so a hexagon covers
  # (sqrt(3)/2) * cellsize^2 -- less than the square of the same cellsize.
  # Applying the factor exactly once therefore makes a hex cell cover the
  # same ground as a square cell built for the same target: A / target_cells.
  #
  # Applying it zero times (or inverting it) gives 1 / (sqrt(3)/2) = 1.155x
  # that area; applying it twice gives 0.866x.  Both fall outside the band
  # below, which the correct implementation clears at every target from 100
  # to 2000 (observed range 0.96-1.08; nx/ny rounding is the only slack).
  bnd  <- .square_boundary(100)
  area <- as.numeric(sf::st_area(bnd))

  for (target in c(100, 400, 900)) {
    # clip = FALSE so every cell is whole and edge partials cannot bias the
    # median.
    hex <- create_grid_polygons(bnd, target_cells = target, type = "hex",
                                clip = FALSE, quiet = TRUE)
    ideal <- area / target
    got   <- stats::median(as.numeric(sf::st_area(hex)))
    expect_gt(got / ideal, 0.93)
    expect_lt(got / ideal, 1.12)
  }
})

test_that("hex and square grids give comparable counts for the same target", {
  bnd <- .square_boundary(100)
  target <- 400

  sq  <- create_grid_polygons(bnd, target_cells = target, type = "square",
                              quiet = TRUE)
  hex <- create_grid_polygons(bnd, target_cells = target, type = "hex",
                              quiet = TRUE)

  # Both grid types aim at the same count over the same boundary, so the
  # ratio's expected value is 1 -- not sqrt(3)/2 and not its reciprocal.
  # Only hex edge overhang pushes it up, and never past 1.15.
  ratio <- nrow(hex) / nrow(sq)
  expect_gt(ratio, 0.9)
  expect_lt(ratio, 1.2)

  # The two grid types must also agree on cell size, which is the packing
  # correction stated without the edge effects.
  sq_u  <- create_grid_polygons(bnd, target_cells = target, type = "square",
                                clip = FALSE, quiet = TRUE)
  hex_u <- create_grid_polygons(bnd, target_cells = target, type = "hex",
                                clip = FALSE, quiet = TRUE)
  area_ratio <- stats::median(as.numeric(sf::st_area(hex_u))) /
    stats::median(as.numeric(sf::st_area(sq_u)))
  expect_gt(area_ratio, 0.92)
  expect_lt(area_ratio, 1.12)
})
