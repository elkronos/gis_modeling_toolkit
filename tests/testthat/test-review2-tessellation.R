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


# ---------------------------------------------------------------------------
# One of points and boundary without a CRS
# ---------------------------------------------------------------------------

test_that("CRS-less planar points tessellate with a CRS-less boundary", {
  # ensure_projected() marks planar CRS-less input crs_assumed = "none", and
  # build_tessellation() read that as a CRS: st_crs("none") is an error, so
  # every method failed, including the documented clip_target_for() pattern,
  # and CRS-less planar data could not be gridded at all.
  pts <- .r2_pts_na()
  bnd <- clip_target_for(pts, quiet = TRUE)
  for (m in c("voronoi", "triangles", "hex", "square")) {
    tess <- build_tessellation(pts, boundary = bnd, method = m, quiet = TRUE,
                               approx_n_cells = if (m %in% c("hex", "square")) 20)
    expect_gt(nrow(tess$cells), 0L)
    expect_true(is.na(sf::st_crs(tess$cells)), info = m)
    expect_false(anyNA(tess$index), info = m)
  }
})


test_that("a CRS-less boundary takes the points' CRS in every builder", {
  pts <- sf::st_set_crs(.r2_pts_na(), 32632)
  bnd <- .r2_sq(5e5 - 10, 5e6 - 10, 1020)            # same numbers, no CRS

  # create_voronoi_polygons() died on sf's "st_crs(x) == st_crs(y) is not TRUE".
  expect_warning(v <- create_voronoi_polygons(pts, boundary = bnd, quiet = TRUE),
                 "stamping")
  expect_equal(nrow(v$cells), 30L)
  expect_equal(sf::st_crs(v$cells), sf::st_crs(32632))

  # clip_target_for() returned a target with no CRS.
  expect_warning(tgt <- clip_target_for(pts, boundary = bnd, quiet = TRUE), "stamping")
  expect_equal(sf::st_crs(tgt), sf::st_crs(32632))

  # build_tessellation(method = "triangles") failed like create_voronoi_polygons().
  expect_warning(tri <- build_tessellation(pts, boundary = bnd, method = "triangles",
                                           quiet = TRUE), "stamping")
  expect_gt(nrow(tri$cells), 0L)

  # Lon/lat points with a CRS-less lon/lat boundary: the target came back in
  # degrees with no CRS, and expand = 20 buffered by 20 degrees.
  pts_ll <- sf::st_transform(pts, 4326)
  bnd_ll <- sf::st_set_crs(sf::st_transform(sf::st_set_crs(bnd, 32632), 4326), NA)
  expect_warning(tgt_ll <- clip_target_for(pts_ll, boundary = bnd_ll, expand = 20,
                                           quiet = TRUE), "look like lon/lat")
  expect_false(is.na(sf::st_crs(tgt_ll)))
  expect_false(sf::st_is_longlat(tgt_ll))
  bb <- sf::st_bbox(sf::st_transform(tgt_ll, 32632))
  expect_equal(as.numeric(bb["xmax"] - bb["xmin"]), 1020 + 2 * 20, tolerance = 1e-3)
})


test_that("CRS-less points take a projected boundary's CRS, and refuse a geographic one", {
  pts <- .r2_pts_na()                                  # UTM numbers from a CSV
  bnd <- .r2_sq(5e5 - 10, 5e6 - 10, 1020, crs = 32632)
  for (m in c("voronoi", "triangles", "hex", "square")) {
    expect_warning(
      tess <- build_tessellation(pts, boundary = bnd, method = m, quiet = TRUE,
                                 approx_n_cells = if (m %in% c("hex", "square")) 20),
      "stamping", info = m)
    expect_equal(sf::st_crs(tess$cells), sf::st_crs(32632), info = m)
    expect_false(anyNA(tess$index), info = m)
  }
  expect_warning(v <- create_voronoi_polygons(pts, boundary = bnd, quiet = TRUE),
                 "stamping")
  expect_equal(nrow(v$cells), 30L)

  # Stamping degrees on UTM numbers would be wrong: refused, with a reason.
  expect_error(
    suppressWarnings(build_tessellation(pts, boundary = sf::st_transform(bnd, 4326),
                                        method = "voronoi", quiet = TRUE)),
    "has no CRS and its coordinates do not look like lon/lat")
})


# ---------------------------------------------------------------------------
# A geographic `crs`
# ---------------------------------------------------------------------------

test_that("a geographic crs returns cells built in metres", {
  set.seed(3)
  pts <- sf::st_as_sf(data.frame(x = stats::runif(150, 10, 14),
                                 y = stats::runif(150, 53, 57)),
                      coords = c("x", "y"), crs = 4326)

  # Voronoi: every location in a cell is nearest (on the ground) to that
  # cell's point.  Built in degrees, about 20% were not.
  tv <- build_tessellation(pts, method = "voronoi", crs = 4326, quiet = TRUE)
  expect_equal(sf::st_crs(tv$cells), sf::st_crs(4326))
  cc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(tv$cells)))
  expect_true(all(abs(cc[, 1]) <= 180 & abs(cc[, 2]) <= 90))
  probe <- sf::st_as_sf(expand.grid(x = seq(10.2, 13.8, by = 0.1),
                                    y = seq(53.2, 56.8, by = 0.1)),
                        coords = c("x", "y"), crs = 4326)
  hit  <- sf::st_intersects(probe, tv$cells)
  one  <- lengths(hit) == 1L
  cell <- tv$cells$cell_id[unlist(hit[one])]
  near <- sf::st_nearest_feature(sf::st_transform(probe[one, ], 32632),
                                 sf::st_transform(pts, 32632))
  expect_lt(mean(tv$index[near] != cell), 0.01)

  # Grids sized by a count: no s2 "degenerate edge" error, and every point
  # inside the boundary has a cell.
  bnd <- .r2_box_ll(10, 14, 53, 57)
  for (m in c("hex", "square")) {
    tg <- build_tessellation(pts, boundary = bnd, method = m, approx_n_cells = 80,
                             crs = 4326, quiet = TRUE)
    expect_equal(sf::st_crs(tg$cells), sf::st_crs(4326), info = m)
    expect_false(anyNA(tg$index), info = m)
  }
  # create_grid_polygons(): the grid laid in the projected CRS it picks on its
  # own, returned in lon/lat, rather than a different grid laid in degrees.
  g   <- create_grid_polygons(bnd, target_cells = 50, crs = 4326, quiet = TRUE)
  ref <- create_grid_polygons(bnd, target_cells = 50, quiet = TRUE)
  expect_equal(sf::st_crs(g), sf::st_crs(4326))
  expect_equal(nrow(g), nrow(ref))
  expect_equal(sort(as.numeric(sf::st_area(sf::st_transform(g, sf::st_crs(ref))))),
               sort(as.numeric(sf::st_area(ref))), tolerance = 1e-4)
})


# ---------------------------------------------------------------------------
# Grid sizing, expand, collinear points
# ---------------------------------------------------------------------------

test_that("hex cell size does not depend on which way the boundary lies", {
  # st_make_grid() sizes hexagons by cellsize[1], which was w / nx: a tall
  # strip got hexagons as wide as the whole strip, 1734 of them against 89.
  rect <- function(w, h) sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(w, 0), c(w, h), c(0, h), c(0, 0)))), crs = 32632))
  hex_area <- function(g) max(as.numeric(sf::st_area(g)))   # a whole hexagon
  for (d in list(c(1000, 1, 9), c(3000, 1000, 100), c(5000, 700, 30))) {
    flat <- create_grid_polygons(rect(d[1], d[2]), target_cells = d[3], type = "hex",
                                 clip = FALSE)
    tall <- create_grid_polygons(rect(d[2], d[1]), target_cells = d[3], type = "hex",
                                 clip = FALSE)
    expect_equal(hex_area(tall), hex_area(flat), tolerance = 1e-9,
                 info = paste(d, collapse = " x "))
  }
  # The count follows (hexagons are not symmetric under a quarter turn, so
  # not exactly): 89 and 103 for the strip.
  flat <- nrow(create_grid_polygons(rect(1000, 1), target_cells = 9, type = "hex"))
  tall <- nrow(create_grid_polygons(rect(1, 1000), target_cells = 9, type = "hex"))
  expect_lt(max(flat, tall) / min(flat, tall), 1.3)

  # A boundary at least as wide as it is tall keeps the grid it had.
  g <- create_grid_polygons(rect(3000, 1000), target_cells = 100, type = "hex")
  eff <- 100 * sqrt(3) / 2
  expect_equal(hex_area(g), sqrt(3) / 2 * (3000 / round(sqrt(eff * 3)))^2,
               tolerance = 1e-9)
})


test_that("voronoi with expand > 0 returns the grown boundary its cells fill", {
  sq <- .r2_sq(0, 0, 1000, crs = 32632)
  set.seed(5)
  p <- sf::st_as_sf(rbind(data.frame(x = stats::runif(30, 0, 1000),
                                     y = stats::runif(30, 0, 1000)),
                          data.frame(x = c(-100, 1150), y = c(500, 500))),
                    coords = c("x", "y"), crs = 32632)
  tess <- build_tessellation(p, boundary = sq, expand = 200, quiet = TRUE)
  expect_false(anyNA(tess$index))      # the two outside points lie within 200 m
  expect_equal(as.numeric(sum(sf::st_area(tess$boundary))),
               as.numeric(sum(sf::st_area(tess$cells))), tolerance = 1e-6)
  expect_gt(as.numeric(sum(sf::st_area(tess$boundary))), 1.5e6)
  # expand = 0 is unchanged.
  t0 <- build_tessellation(p, boundary = sq, quiet = TRUE)
  expect_equal(as.numeric(sf::st_area(t0$boundary)), 1e6)
})


test_that("collinear points are refused for triangles, not returned empty", {
  col <- sf::st_as_sf(data.frame(x = 5e5 + 0:9 * 10, y = 5e6 + 0:9 * 5),
                      coords = c("x", "y"), crs = 32632)
  expect_error(build_tessellation(col, method = "triangles", quiet = TRUE),
               "collinear")
  # Voronoi handles a transect.
  v <- build_tessellation(col, method = "voronoi", quiet = TRUE)
  expect_equal(nrow(v$cells), 10L)
  # One point off the line is enough.
  col2 <- col; sf::st_geometry(col2)[[5]] <- sf::st_point(c(5e5 + 40, 5e6 + 21))
  expect_gt(nrow(build_tessellation(col2, method = "triangles", quiet = TRUE)$cells), 0L)
})


# ---------------------------------------------------------------------------
# Lattices from a profile's count
# ---------------------------------------------------------------------------

test_that("a lattice sized by a selection reports and warns about empty cells", {
  set.seed(21)
  ctrs <- matrix(stats::runif(12, 1000, 9000), ncol = 2)
  xy <- do.call(rbind, lapply(1:6, function(i)
    cbind(stats::rnorm(100, ctrs[i, 1], 300), stats::rnorm(100, ctrs[i, 2], 300))))
  clustered <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
                            coords = c("x", "y"), crs = 32632)
  sel <- structure(list(best = 56, criterion = "cp", edge = NA_character_),
                   class = "resolution_selection")
  bnd <- clip_target_for(clustered, quiet = TRUE)
  expect_warning(
    tess <- build_tessellation(clustered, boundary = bnd, method = "hex",
                               approx_n_cells = sel, quiet = TRUE),
    "of the .* cells hold a point")
  expect_lt(tess$params$cells_occupied, 0.75 * 56)
  expect_equal(tess$params$cells_occupied + tess$params$cells_empty, nrow(tess$cells))

  # Evenly spread points fill the lattice: no warning.
  set.seed(22)
  even <- sf::st_as_sf(data.frame(x = stats::runif(600, 0, 10000),
                                  y = stats::runif(600, 0, 10000)),
                       coords = c("x", "y"), crs = 32632)
  expect_no_warning(t2 <- build_tessellation(even, boundary = clip_target_for(even, quiet = TRUE),
                                             method = "hex", approx_n_cells = sel,
                                             quiet = TRUE))
  expect_gte(t2$params$cells_occupied, 0.75 * 56)
  # A plain number records occupancy but does not warn.
  expect_no_warning(t3 <- build_tessellation(clustered, boundary = bnd, method = "hex",
                                             approx_n_cells = 56, quiet = TRUE))
  expect_true(is.numeric(t3$params$cells_occupied))
})
