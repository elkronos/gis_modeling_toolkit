# ===========================================================================
# Regressions from the second adversarial review of tessellation, CRS
# selection and seeding.  Each section says what used to go wrong.
# ===========================================================================


.r2_sq <- function(x0, y0, s, crs = sf::NA_crs_) {
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0)))),
    crs = crs))
}

.r2_box_ll <- function(x0, x1, y0, y1) {
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)))), crs = 4326))
}

# 30 points over a kilometre at UTM-sized coordinates, with no CRS.
.r2_pts_na <- function(n = 30L, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = 5e5 + stats::runif(n, 0, 1000),
                          y = 5e6 + stats::runif(n, 0, 1000)),
               coords = c("x", "y"))
}

# Evaluate `code` with sf_use_s2() set to `on`, restoring the session's value.
.r2_with_s2 <- function(on, code) {
  was <- suppressMessages(sf::sf_use_s2(on))
  on.exit(suppressMessages(sf::sf_use_s2(was)), add = TRUE)
  force(code)
}


# ---------------------------------------------------------------------------
# Choosing the projection
# ---------------------------------------------------------------------------

test_that("ensure_projected() accepts a lon/lat polygon s2 rejects", {
  # A repeated vertex is valid for GEOS and a plain st_transform(), but s2's
  # st_centroid(st_union()) stopped with "Edge 1 is degenerate".
  rep_v <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(8, 47), c(9, 47), c(9, 47), c(9, 48), c(8, 48), c(8, 47)))), crs = 4326))
  expect_equal(sf::st_crs(ensure_projected(rep_v))$epsg, 32632L)
  expect_false(sf::st_is_longlat(ensure_projected(rep_v, purpose = "area")))
  g <- create_grid_polygons(rep_v, target_cells = 20, quiet = TRUE)
  expect_gt(nrow(g), 0L)
})


test_that("a single study-area polygon is scored, not left to the UTM fallback", {
  # One polygon reduced to one point: every candidate scored NA and CONUS
  # stayed in UTM zone 15 (13% worst-case distance error).
  conus <- .r2_box_ll(-124, -67, 25, 49)
  out <- ensure_projected(conus)
  ch  <- attr(out, "crs_choice")
  expect_true(all(is.finite(ch$distance_error)))
  expect_false(grepl("UTM", ch$name[ch$chosen]))
  expect_lt(ch$distance_error[ch$chosen], 0.05)
  expect_gt(ch$distance_error[grepl("UTM", ch$name)], 0.05)
})


test_that("the UTM zone does not depend on sf_use_s2()", {
  # Two clusters whose centroid on the sphere lies west of -78 (zone 17) and
  # whose planar centroid in degrees lies east of it (zone 18).
  two <- sf::st_as_sf(data.frame(x = c(-78.6, -78.59, -77.2, -77.19),
                                 y = c(0, 0.01, 60, 60.01)),
                      coords = c("x", "y"), crs = 4326)
  z_on  <- .r2_with_s2(TRUE,  sf::st_crs(ensure_projected(two))$epsg)
  z_off <- .r2_with_s2(FALSE, sf::st_crs(ensure_projected(two))$epsg)
  expect_identical(z_on, 32617L)
  expect_identical(z_off, z_on)
})


test_that("with s2 off, lon/lat input raises none of sf's centroid conditions", {
  pts_ll <- sf::st_as_sf(data.frame(x = c(9.1, 9.2, 9.15, 9.3),
                                    y = c(48.7, 48.8, 48.75, 48.72)),
                         coords = c("x", "y"), crs = 4326)
  .r2_with_s2(FALSE, {
    expect_no_warning(expect_no_message(clip_target_for(pts_ll, quiet = TRUE)))
    expect_no_warning(expect_no_message(ensure_projected(pts_ll)))
  })
})


test_that("grids over a near-global boundary are equal-area; local ones keep UTM", {
  glob <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    cbind(seq(-170, 170, 10), -60), cbind(seq(170, -170, -10), 70), c(-170, -60)))),
    crs = 4326))
  # Web Mercator, which the distance choice falls back to here, made the
  # full cells' true areas differ almost five-fold.
  g <- create_grid_polygons(glob, target_cells = 100, quiet = TRUE)
  expect_false(grepl("+proj=merc", sf::st_crs(g)$proj4string, fixed = TRUE))
  expect_lt(.crs_area_error(g), 0.01)
  # The cached builder lays the same grid in the same CRS.
  gc <- create_grid_polygons_cached(glob, target_cells = 100, cache_env = new.env())
  expect_equal(sf::st_crs(gc), sf::st_crs(g))

  # A local boundary is unchanged: its zone.
  loc <- create_grid_polygons(.r2_box_ll(9, 10, 47, 48), target_cells = 50, quiet = TRUE)
  expect_equal(sf::st_crs(loc)$epsg, 32632L)
})


# ---------------------------------------------------------------------------
# sf_use_s2(FALSE) without lwgeom
# ---------------------------------------------------------------------------

test_that("with s2 off the package needs no lwgeom", {
  set.seed(1)
  pts <- sf::st_as_sf(data.frame(x = 5e5 + stats::runif(15, 0, 100),
                                 y = 5e6 + stats::runif(15, 0, 100)),
                      coords = c("x", "y"), crs = 32632)
  bnd <- .r2_sq(0, 0, 100, crs = 32632)
  bnd_ll <- .r2_box_ll(8, 9, 47, 48)
  on <- .r2_with_s2(TRUE, list(
    ids = ensure_stable_poly_id(create_grid_polygons(bnd, target_cells = 9))))
  .r2_with_s2(FALSE, {
    expect_equal(nrow(create_voronoi_polygons(pts, quiet = TRUE)$cells), 15L)
    expect_equal(nrow(build_tessellation(pts, method = "voronoi", quiet = TRUE)$cells), 15L)
    ids <- ensure_stable_poly_id(create_grid_polygons(bnd, target_cells = 9))
    # The same IDs as with s2 on: the sort key is measured the same way.
    expect_identical(sf::st_as_text(sf::st_geometry(ids)),
                     sf::st_as_text(sf::st_geometry(on$ids)))
    expect_equal(nrow(create_grid_polygons_cached(bnd, target_cells = 9,
                                                  cache_env = new.env())), 9L)
    expect_true(is.finite(.crs_area_error(create_grid_polygons(bnd, target_cells = 9))))
    expect_equal(nrow(suppressWarnings(voronoi_seeds_random(bnd_ll, k = 5, set_seed = 1))), 5L)
    expect_equal(nrow(suppressWarnings(
      get_voronoi_seeds(bnd_ll, method = "kmeans", n = 3, set_seed = 1))), 3L)
    # And the session's setting is left as it was.
    expect_false(sf::sf_use_s2())
  })
})
