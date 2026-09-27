# tests/testthat/test-review3-S2-tessellation.R
# ---------------------------------------------------------------------------
# Regressions from the third review of tessellation, CRS selection and
# seeding.  Each test names the finding it closes and says what used to go
# wrong.
# ---------------------------------------------------------------------------

.r3s2_box <- function(x0, x1, y0, y1, crs = sf::NA_crs_) {
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)))), crs = crs))
}

.r3s2_pts <- function(n, x0, x1, y0, y1, crs = sf::NA_crs_, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = stats::runif(n, x0, x1), y = stats::runif(n, y0, y1)),
               coords = c("x", "y"), crs = crs)
}

.r3s2_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}


# ---------------------------------------------------------------------------
# S2-TESSELLATION-1: hex/square lattices over lon/lat data
# ---------------------------------------------------------------------------

test_that("build_tessellation() lays a near-global lattice in the CRS create_grid_polygons() picks", {
  # A near-global boundary: the distance CRS for the points is Web Mercator,
  # whose full hexagons differed 5.75-fold in true area, while
  # create_grid_polygons() on the same boundary used Equal Earth.  (Vertices
  # every 5 degrees, so the outline does not read as a band across 180.)
  lon <- seq(-170, 170, by = 5); lat <- seq(-60, 70, by = 5)
  ring <- rbind(cbind(lon, -60), cbind(170, lat[-1]), cbind(rev(lon)[-1], 70),
                cbind(-170, rev(lat)[-1]))
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326))
  pts <- .r3s2_pts(80, -165, 165, -55, 65, crs = 4326, seed = 3)
  ref <- create_grid_polygons(bnd, target_cells = 60, type = "hex", quiet = TRUE)
  expect_match(sf::st_crs(ref)$proj4string, "eqearth")
  for (m in c("hex", "square")) {
    tb <- build_tessellation(pts, bnd, method = m, approx_n_cells = 60, quiet = TRUE)
    expect_equal(sf::st_crs(tb$cells), sf::st_crs(ref), info = m)
    expect_equal(sf::st_crs(tb$boundary), sf::st_crs(ref), info = m)
    expect_false(anyNA(tb$index), info = m)
  }
  # The same grid as create_grid_polygons() lays.
  th <- build_tessellation(pts, bnd, method = "hex", approx_n_cells = 60, quiet = TRUE)
  expect_equal(nrow(th$cells), nrow(ref))
  # Returned in lon/lat on request, and still laid equal-area.
  tg <- build_tessellation(pts, bnd, method = "hex", approx_n_cells = 60,
                           crs = 4326, quiet = TRUE)
  expect_equal(sf::st_crs(tg$cells), sf::st_crs(4326))
  expect_equal(nrow(tg$cells), nrow(ref))
  expect_false(anyNA(tg$index))

  # A local extent keeps the points' UTM zone.
  loc <- build_tessellation(.r3s2_pts(40, 9.1, 9.9, 47.1, 47.9, crs = 4326),
                            .r3s2_box(9, 10, 47, 48, crs = 4326), method = "hex",
                            approx_n_cells = 20, quiet = TRUE)
  expect_equal(sf::st_crs(loc$cells)$epsg, 32632L)
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-2: an EMPTY point with method = "triangles"
# ---------------------------------------------------------------------------

test_that("an empty point is indexed NA by the triangulation, not an abort", {
  set.seed(7)
  g <- sf::st_sfc(c(lapply(1:20, function(i)
    sf::st_point(c(5e5 + stats::runif(1, 0, 1000), 5e6 + stats::runif(1, 0, 1000)))),
    list(sf::st_point())), crs = 32632)
  pts <- sf::st_sf(id = 1:21, geometry = g)
  bnd <- .r3s2_box(5e5, 5e5 + 1000, 5e6, 5e6 + 1000, crs = 32632)
  # Stopped with "NA/NaN/Inf in foreign function call (arg 1)".
  for (b in list(bnd, NULL)) {
    tri <- build_tessellation(pts, b, method = "triangles", quiet = TRUE)
    expect_gt(nrow(tri$cells), 0L)
    expect_length(tri$index, 21L)
    expect_true(is.na(tri$index[21]))
    expect_false(anyNA(tri$index[1:20]))
  }
  # Two real points and an empty one are still too few.
  expect_error(build_tessellation(pts[c(1, 2, 21), ], method = "triangles", quiet = TRUE),
               "at least 3 unique points")
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-3: near-degenerate extents and the max_cells message
# ---------------------------------------------------------------------------

test_that("clip_target_for() buffers a transect with sub-millimetre scatter", {
  set.seed(2)
  p <- sf::st_as_sf(data.frame(x = seq(0, 1000, length.out = 30),
                               y = stats::runif(30, 0, 1e-6)),
                    coords = c("x", "y"), crs = 32632)
  ct <- clip_target_for(p, quiet = TRUE)
  bb <- sf::st_bbox(ct)
  # Was a 1000 x 9e-7 sliver, over which a count of 25 built 166,536 squares.
  expect_gt(as.numeric(bb["ymax"] - bb["ymin"]), 1)
  sq <- build_tessellation(p, ct, method = "square", approx_n_cells = 25, quiet = TRUE)
  expect_lt(nrow(sq$cells), 100L)
})

test_that("the max_cells error names the argument that set the cell size", {
  sliver <- .r3s2_box(0, 1000, 0, 1.8e-8, crs = 32632)
  e <- tryCatch(create_grid_polygons(sliver, target_cells = 25, quiet = TRUE),
                error = conditionMessage)
  expect_match(e, "max_cells")
  expect_match(e, "`target_cells` = 25", fixed = TRUE)
  expect_no_match(e, "Check that `cellsize`", fixed = TRUE)
  sq <- .r3s2_box(0, 1000, 0, 1000, crs = 32632)
  expect_error(create_grid_polygons(sq, n = 2000, quiet = TRUE), "derived from `n` = 2000 x 2000")
  expect_error(create_grid_polygons(sq, cellsize = 0.5, quiet = TRUE),
               "Check that `cellsize` is in the boundary's CRS units")
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-4 / -14: a CRS-less boundary with lon/lat points
# ---------------------------------------------------------------------------

test_that("a CRS-less degree tile is read in the lon/lat points' CRS, not the projected one", {
  pts <- .r3s2_pts(50, 10.1, 10.9, 50.1, 50.9, crs = 4326, seed = 8)
  tile <- .r3s2_box(10, 11, 50, 51)        # integer corners: the heuristic declines it
  expect_false(.looks_like_lonlat(tile)$lonlat)
  # It was stamped with the UTM zone picked for the points: a one-metre square,
  # every point indexed NA.
  for (m in c("voronoi", "hex")) {
    w <- .r3s2_warnings(build_tessellation(pts, tile, method = m, quiet = TRUE,
                                           approx_n_cells = if (m == "hex") 20))
    expect_true(any(grepl("look like lon/lat", w$warnings)), info = m)
    expect_false(anyNA(w$value$index), info = m)
    expect_false(sf::st_is_longlat(w$value$cells), info = m)
  }
  v <- suppressWarnings(create_voronoi_polygons(pts, tile, quiet = TRUE))
  expect_false(anyNA(v$index))
  tgt <- suppressWarnings(clip_target_for(pts, tile, quiet = TRUE))
  expect_true(all(lengths(sf::st_intersects(sf::st_transform(pts, sf::st_crs(tgt)), tgt)) > 0))
})

test_that("a CRS-less boundary in metres with lon/lat points is refused, naming both", {
  pts <- .r3s2_pts(30, -1.9, -1.1, 52.1, 52.9, crs = 4326)
  bng <- sf::st_set_crs(sf::st_transform(.r3s2_box(-2, -1, 52, 53, crs = 4326), 27700), NA)
  for (m in c("voronoi", "hex"))
    expect_error(build_tessellation(pts, bng, method = m, quiet = TRUE,
                                    approx_n_cells = if (m == "hex") 20),
                 "cannot be placed in one space", info = m)
  expect_error(create_voronoi_polygons(pts, bng, quiet = TRUE), "cannot be placed in one space")
  expect_error(clip_target_for(pts, bng, quiet = TRUE), "cannot be placed in one space")
})

test_that("CRS-less lon/lat points with a CRS-less boundary in metres are refused (S2-14)", {
  pts <- .r3s2_pts(60, 12.2, 13.8, 54.2, 55.8, seed = 5)
  bnd <- sf::st_set_crs(sf::st_transform(.r3s2_box(12, 14, 54, 56, crs = 4326), 32633), NA)
  # Was stamped EPSG:4326, transformed to nothing, and refused as "not polygonal".
  for (m in c("voronoi", "hex"))
    expect_error(suppressWarnings(build_tessellation(pts, bnd, method = m, quiet = TRUE,
                                                     approx_n_cells = if (m == "hex") 20)),
                 "has no CRS either and was taken as lon/lat", info = m)
  # A CRS-less lon/lat boundary still takes the same assumption.
  ok <- suppressWarnings(build_tessellation(pts, .r3s2_box(12, 14, 54, 56),
                                            method = "voronoi", quiet = TRUE))
  expect_equal(sf::st_crs(ok$cells)$epsg, 32633L)
  expect_false(anyNA(ok$index))
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-5: create_voronoi_polygons() and CRS-less degrees
# ---------------------------------------------------------------------------

test_that("create_voronoi_polygons() applies the lon/lat heuristic to CRS-less points", {
  pts <- .r3s2_pts(60, 12.2, 13.8, 54.2, 55.8, seed = 5)
  # Built silently in degrees before, with CRS NA.
  expect_warning(v <- create_voronoi_polygons(pts, quiet = TRUE), "assuming EPSG:4326")
  bt <- suppressWarnings(build_tessellation(pts, method = "voronoi", quiet = TRUE))
  expect_equal(sf::st_crs(v$cells), sf::st_crs(bt$cells))
  expect_false(sf::st_is_longlat(v$cells))
  expect_identical(v$index, bt$index)
  # And with a CRS-less boundary, the same result as build_tessellation().
  bnd <- .r3s2_box(12, 14, 54, 56)
  v2 <- suppressWarnings(create_voronoi_polygons(pts, bnd, quiet = TRUE))
  b2 <- suppressWarnings(build_tessellation(pts, bnd, method = "voronoi", quiet = TRUE))
  expect_equal(sf::st_crs(v2$cells), sf::st_crs(b2$cells))
  expect_identical(v2$index, b2$index)
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-7: triangles returned in a geographic crs
# ---------------------------------------------------------------------------

test_that("every point lies in its indexed triangle on the returned lon/lat layer", {
  pts <- .r3s2_pts(80, 9, 11, 54, 56, crs = 4326, seed = 4)
  res <- build_tessellation(pts, method = "triangles", crs = 4326, quiet = TRUE)
  expect_equal(sf::st_crs(res$cells), sf::st_crs(4326))
  hit <- sf::st_intersects(pts, res$cells)
  expect_true(all(lengths(hit) > 0))
  expect_true(all(mapply(function(h, id) id %in% res$cells$cell_id[h], hit, res$index)))
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-8 / -9: scoring a projection on an outline
# ---------------------------------------------------------------------------

test_that("a detailed outline is scored on a bounded sample", {
  n  <- 1e5
  th <- seq(0, 2 * pi, length.out = n)
  ring <- cbind(10 + 8 * cos(th), 50 + 3 * sin(th)); ring[n, ] <- ring[1, ]
  poly <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326))
  e <- .crs_distance_error(poly, sf::st_crs(32632))
  expect_true(is.finite(e))
  expect_lt(e, 0.05)
})

test_that("an outline across the antimeridian is measured along its own edges", {
  fiji <- .r3s2_box(177, -178, -19, -15, crs = 4326)
  ch <- attr(ensure_projected(fiji), "crs_choice")
  # Reported 164%: the 177 -> -178 edge was filled with points near 0E.
  expect_lt(ch$distance_error[ch$chosen], 0.01)
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-11: the projection message names the CRS
# ---------------------------------------------------------------------------

test_that("build_tessellation() names the CRS it projects to", {
  pts <- .r3s2_pts(20, 9.1, 9.9, 47.1, 47.9, crs = 4326)
  msgs <- testthat::capture_messages(build_tessellation(pts, method = "voronoi"))
  expect_true(any(grepl("projecting points to EPSG:32632", msgs, fixed = TRUE)))
  expect_false(any(grepl("local UTM CRS", msgs, fixed = TRUE)))
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-13: an sf layer as `crs`
# ---------------------------------------------------------------------------

test_that("the builders take a layer's CRS as `crs`", {
  pts <- .r3s2_pts(30, 9, 10, 50, 51, crs = 4326, seed = 10)
  ref <- sf::st_transform(pts, 32632)
  # Stopped with "the condition has length > 1".
  a <- build_tessellation(pts, method = "voronoi", crs = ref, quiet = TRUE)
  b <- build_tessellation(pts, method = "voronoi", crs = 32632, quiet = TRUE)
  expect_identical(a$index, b$index)
  expect_equal(sf::st_crs(a$cells), sf::st_crs(32632))
  v <- create_voronoi_polygons(pts, crs = sf::st_geometry(ref), quiet = TRUE)
  expect_equal(sf::st_crs(v$cells), sf::st_crs(32632))
  box <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(ref)))
  g <- create_grid_polygons(box, target_cells = 9, crs = ref[1, ], quiet = TRUE)
  expect_equal(sf::st_crs(g), sf::st_crs(32632))
})


# ---------------------------------------------------------------------------
# S2-TESSELLATION-15: lon/lat seeding without lwgeom
# ---------------------------------------------------------------------------

test_that("random and k-means seeding on a lon/lat boundary raise no lwgeom warning", {
  # Looked up with system.file() rather than requireNamespace(), which R CMD
  # check counts as a use of a package DESCRIPTION would have to declare.
  skip_if(nzchar(system.file(package = "lwgeom")),
          "with lwgeom installed sf does not warn")
  bnd <- .r3s2_box(9, 11, 54, 56, crs = 4326)
  for (s2 in c(TRUE, FALSE)) {
    was <- suppressMessages(sf::sf_use_s2(s2))
    expect_no_warning(r <- voronoi_seeds_random(bnd, 5, set_seed = 1))
    expect_no_warning(k <- get_voronoi_seeds(bnd, method = "kmeans", n = 5, set_seed = 1))
    # With s2 off sf also printed its planar-union message.
    expect_no_message(get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 1))
    suppressMessages(sf::sf_use_s2(was))
    expect_equal(nrow(r), 5L)
    expect_equal(nrow(k), 5L)
  }
})


# ---------------------------------------------------------------------------
# S11-contracts-13 / -17 / -18
# ---------------------------------------------------------------------------

test_that("build_tessellation() says how to proceed with polygon features", {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)[1:5, ]
  e <- tryCatch(build_tessellation(nc, method = "voronoi", quiet = TRUE),
                error = conditionMessage)
  expect_match(e, "geometry must be one of: POINT, MULTIPOINT", fixed = TRUE)
  expect_match(e, "coerce_to_points(points_sf, \"auto\")", fixed = TRUE)
})

test_that("a whole build_tessellation() result is named where a polygon layer is expected", {
  pts <- .r3s2_pts(20, 5e5, 5e5 + 100, 5e6, 5e6 + 100, crs = 32632)
  tess <- build_tessellation(pts, method = "voronoi", quiet = TRUE)
  hint <- "this looks like a build_tessellation() result"
  # These stopped with sf's "no applicable method for 'st_crs<-' applied to
  # an object of class \"list\"" after a warning about stamping a CRS.
  for (f in list(function() build_tessellation(pts, tess, method = "voronoi", quiet = TRUE),
                 function() create_voronoi_polygons(pts, tess, quiet = TRUE),
                 function() clip_target_for(pts, tess, quiet = TRUE))) {
    w <- character(0)
    e <- tryCatch(withCallingHandlers(f(), warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }), error = conditionMessage)
    expect_match(e, hint, fixed = TRUE)
    expect_match(e, "`$boundary`", fixed = TRUE)
    expect_length(w, 0L)
  }
  expect_error(create_grid_polygons(tess, target_cells = 9, quiet = TRUE), hint, fixed = TRUE)
  expect_error(create_grid_polygons_cached(tess, target_cells = 9, cache_env = new.env()),
               hint, fixed = TRUE)
  expect_error(ensure_stable_poly_id(tess), "pass its `$cells`", fixed = TRUE)
})

test_that("triangles do not record the approx_n_cells they ignored", {
  skip_if_not_installed("geometry")
  pts <- .r3s2_pts(20, 0, 1, 0, 1, crs = 32632)
  w <- .r3s2_warnings(build_tessellation(pts, method = "triangles",
                                         approx_n_cells = 20, quiet = TRUE))
  expect_null(w$value$params$approx_n_cells)
  expect_null(w$value$params$approx_n_cells_from)
  expect_true(any(grepl("`params` does not record the request", w$warnings, fixed = TRUE)))
})
