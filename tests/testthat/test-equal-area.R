# tests/testthat/test-equal-area.R
# ---------------------------------------------------------------------------
# ensure_projected(purpose = "area") and summarize_by_cell(area = TRUE): the
# CRS a density is computed in has to be equal-area, and that property is
# measured (planar against geodesic areas over probe polygons) rather than
# read off the projection's name.
# ---------------------------------------------------------------------------

ea_local <- function(n = 200, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(data.frame(lon = runif(n, 8, 10.5), lat = runif(n, 47, 49.5)),
               coords = c("lon", "lat"), crs = 4326)
}
ea_conus <- function(n = 200, seed = 2) {
  set.seed(seed)
  sf::st_as_sf(data.frame(lon = runif(n, -124, -67), lat = runif(n, 25, 49)),
               coords = c("lon", "lat"), crs = 4326)
}
area_err <- function(x, crs = NULL) spatialkit:::.crs_area_error(x, crs)


test_that("the area-distortion measurement separates equal-area from conformal projections", {
  loc <- ea_local(); conus <- ea_conus()
  laea_loc <- sf::st_crs("+proj=laea +lat_0=48 +lon_0=9 +datum=WGS84 +units=m +no_defs")
  laea_us  <- sf::st_crs("+proj=laea +lat_0=37 +lon_0=-96 +datum=WGS84 +units=m +no_defs")
  # Inside a zone the conformal error is small; forced over a continent it is
  # not; Web Mercator fails even locally; equal-area stays within the
  # sphere-versus-ellipsoid residual everywhere.
  expect_lt(area_err(loc, sf::st_crs(32632)), 0.003)
  expect_gt(area_err(conus, sf::st_crs(32614)), 0.05)
  expect_gt(area_err(loc, sf::st_crs(3857)), 0.01)
  expect_lt(area_err(loc, laea_loc), 0.005)
  expect_lt(area_err(conus, laea_us), 0.005)
  # Unmeasurable cases are NA, not errors.
  expect_true(is.na(area_err(sf::st_set_crs(loc, NA))))
})

test_that("purpose = 'area' picks an equal-area projection for lon/lat input", {
  loc <- ea_local()
  p_dist <- ensure_projected(loc)
  p_area <- ensure_projected(loc, purpose = "area")
  expect_equal(sf::st_crs(p_dist)$epsg, 32632)
  expect_match(sf::st_crs(p_area)$proj4string, "\\+proj=laea")
  expect_lt(area_err(p_area), 0.005)
  # A continental extent: an equal-area candidate, chosen by measured
  # distance error among the equal-area ones only.
  c_area <- ensure_projected(ea_conus(), purpose = "area")
  expect_match(sf::st_crs(c_area)$proj4string, "\\+proj=(laea|aea)")
  expect_lt(area_err(c_area), 0.01)
  # target_crs still wins.
  expect_equal(sf::st_crs(ensure_projected(loc, target_crs = 3035, purpose = "area"))$epsg, 3035)
  # Global coverage: Equal Earth, not Web Mercator.
  set.seed(3)
  glob <- sf::st_as_sf(data.frame(lon = runif(100, -170, 170), lat = runif(100, -60, 70)),
                       coords = c("lon", "lat"), crs = 4326)
  expect_match(sf::st_crs(ensure_projected(glob, purpose = "area"))$proj4string, "eqearth|moll")
})

test_that("projected input is returned untouched, with a logged warning when it cannot deliver areas", {
  loc_utm   <- sf::st_transform(ea_local(), 32632)
  conus_utm <- sf::st_transform(ea_conus(), 32614)
  lines_ok  <- capture_spatialkit_log(out <- ensure_projected(loc_utm, purpose = "area"))
  expect_identical(sf::st_crs(out), sf::st_crs(loc_utm))
  expect_false(log_has(lines_ok, "distorts areas"))
  lines_bad <- capture_spatialkit_log(out2 <- ensure_projected(conus_utm, purpose = "area"))
  expect_identical(sf::st_crs(out2), sf::st_crs(conus_utm))
  expect_true(log_has(lines_bad, "distorts areas across this extent by up to"))
  expect_true(log_has(lines_bad, "EPSG:32614"))
})

ea_cells <- function(pts_proj, target_cells = 16) {
  bnd  <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts_proj)))
  grid <- create_grid_polygons(bnd, target_cells = target_cells, type = "square")
  list(grid = grid, assigned = assign_features_to_polygons(pts_proj, grid))
}

test_that("summarize_by_cell(area = TRUE) adds the density in an equal-area CRS", {
  loc <- sf::st_transform(ea_local(), 32632)
  cc  <- ea_cells(loc)
  s   <- summarize_by_cell(cc$assigned, cells_sf = cc$grid, area = TRUE)
  expect_s3_class(s, "sf")
  expect_true(all(c("cell_area", "n_per_area") %in% names(s)))
  expect_equal(s$cell_area, as.numeric(sf::st_area(s)))
  expect_equal(s$n_per_area, s$n / s$cell_area)
  expect_true(is.finite(attr(s, "area_error")))
  expect_lt(attr(s, "area_error"), 0.01)
  # Without the flag nothing is added and nothing else changes.
  s0 <- summarize_by_cell(cc$assigned, cells_sf = cc$grid)
  expect_false(any(c("cell_area", "n_per_area") %in% names(s0)))
  expect_equal(sf::st_drop_geometry(s)[, names(sf::st_drop_geometry(s0))],
               sf::st_drop_geometry(s0))
  # A cell with no observations has NA, not 0, for the density.
  sub <- cc$assigned[sf::st_coordinates(cc$assigned)[, 1] < stats::median(sf::st_coordinates(cc$assigned)[, 1]), ]
  s2  <- summarize_by_cell(sub, cells_sf = cc$grid, area = TRUE)
  expect_true(any(is.na(s2$n)))
  expect_true(all(is.na(s2$n_per_area[is.na(s2$n)])))
})

test_that("summarize_by_cell(area = TRUE) refuses a CRS whose areas are not comparable", {
  conus <- sf::st_transform(ea_conus(), 3857)
  cc <- ea_cells(conus)
  expect_error(summarize_by_cell(cc$assigned, cells_sf = cc$grid, area = TRUE),
               "`area = TRUE` refused: the CRS of `cells_sf` \\(EPSG:3857\\) distorts areas")
  expect_error(summarize_by_cell(cc$assigned, cells_sf = cc$grid, area = TRUE),
               "ensure_projected\\(purpose = \"area\"\\)")
  conus_utm <- sf::st_transform(ea_conus(), 32614)
  cu <- ea_cells(conus_utm)
  expect_error(summarize_by_cell(cu$assigned, cells_sf = cu$grid, area = TRUE),
               "EPSG:32614")
  # And the rest of the contract: cells are required, with a CRS, and the
  # flag is logical.
  expect_error(summarize_by_cell(cc$assigned, area = TRUE), "needs `cells_sf`")
  expect_error(summarize_by_cell(cc$assigned, cells_sf = sf::st_set_crs(cc$grid, NA), area = TRUE),
               "carry a CRS")
  expect_error(summarize_by_cell(cc$assigned, cells_sf = cc$grid, area = NA), "`area` must be TRUE or FALSE")
})
