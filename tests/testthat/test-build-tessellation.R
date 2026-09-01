# ===========================================================================
# build_tessellation(): the unified entry point over the four tessellation
# methods.  It had no tests, so several of its branches were only reachable
# in one direction.
#
# The properties every method must satisfy are asserted together in the first
# block, so a new method cannot be added without meeting them.
# ===========================================================================

.bt_sq <- function(x0, y0, x1, y1) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1),
                            c(x0, y0))))
}

.bt_points <- function(n = 30, seed = 2, crs = 32632) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100)),
               coords = c("x", "y"), crs = crs)
}

.bt_boundary <- function(crs = 32632) {
  sf::st_sf(geometry = sf::st_sfc(.bt_sq(0, 0, 100, 100), crs = crs))
}

.bt_all_methods <- c("voronoi", "triangles", "hex", "square")


test_that("every method returns the documented shape and indexes every point", {
  skip_if_not_installed("geometry")
  pts <- .bt_points()
  bnd <- .bt_boundary()

  for (m in .bt_all_methods) {
    res <- suppressMessages(build_tessellation(pts, boundary = bnd, method = m,
                                               approx_n_cells = 9, quiet = TRUE))

    expect_named(res, c("cells", "index", "boundary", "method", "params"),
                 info = m)
    expect_equal(res$method, m, info = m)
    expect_s3_class(res$cells, "sf")
    expect_gt(nrow(res$cells), 0L)
    expect_true(all(as.character(sf::st_geometry_type(res$cells,
                                                      by_geometry = TRUE)) %in%
                      c("POLYGON", "MULTIPOLYGON")), info = m)

    # `cell_id` is the id column EVERY method returns -- it used to be deleted
    # from hex/square grids after indexing, so plot_tessellation_map(fill_col =
    # "cell_id") drew nothing for a hex tessellation.
    expect_true("cell_id" %in% names(res$cells), info = m)
    expect_equal(sort(res$cells$cell_id), seq_len(nrow(res$cells)), info = m)

    # index: one entry per INPUT FEATURE, holding cell_id values (not row
    # indices), with every point placed.
    expect_length(res$index, nrow(pts))
    expect_false(anyNA(res$index), info = m)
    expect_true(all(res$index %in% res$cells$cell_id), info = m)
  }
})


test_that("hex and square keep poly_id alongside cell_id, and they agree", {
  pts <- .bt_points()
  bnd <- .bt_boundary()

  for (m in c("hex", "square")) {
    res <- suppressMessages(build_tessellation(pts, boundary = bnd, method = m,
                                               approx_n_cells = 9, quiet = TRUE))
    # poly_id is retained for backward compatibility and holds the same values.
    expect_true("poly_id" %in% names(res$cells), info = m)
    expect_equal(res$cells$cell_id, res$cells$poly_id, info = m)
  }

  # A square grid over a square boundary hits the requested count exactly.
  sqr <- suppressMessages(build_tessellation(pts, boundary = bnd,
                                             method = "square",
                                             approx_n_cells = 9, quiet = TRUE))
  expect_equal(nrow(sqr$cells), 9L)
  # The Voronoi and triangle methods carry no poly_id.
  vor <- suppressMessages(build_tessellation(pts, boundary = bnd,
                                             method = "voronoi", quiet = TRUE))
  expect_false("poly_id" %in% names(vor$cells))
})


test_that("voronoi gives one cell per unique point", {
  pts <- .bt_points(30)
  bnd <- .bt_boundary()
  res <- suppressMessages(build_tessellation(pts, boundary = bnd,
                                             method = "voronoi", quiet = TRUE))
  expect_equal(nrow(res$cells), 30L)

  # Duplicated positions collapse unless keep_duplicates asks otherwise.
  dup <- rbind(pts, pts[1:5, ])
  dedup <- suppressMessages(build_tessellation(dup, boundary = bnd,
                                               method = "voronoi", quiet = TRUE))
  expect_equal(nrow(dedup$cells), 30L)
  expect_length(dedup$index, 35L)          # still one index entry per feature
})


test_that("triangles is the Delaunay triangulation of the POINT SET", {
  skip_if_not_installed("geometry")
  pts <- .bt_points(30)
  res <- suppressMessages(build_tessellation(pts, method = "triangles",
                                             quiet = TRUE))

  # Euler's formula for a Delaunay triangulation of n points in general
  # position with h of them on the convex hull: 2n - h - 2 triangles.
  n <- nrow(pts)
  h <- nrow(sf::st_coordinates(
    sf::st_convex_hull(sf::st_union(sf::st_geometry(pts))))) - 1L
  expect_equal(nrow(res$cells), 2L * n - h - 2L)
  # Triangulating the convex-hull POLYGON instead -- which the no-qhull path
  # used to do -- discards every interior point and yields only h - 2.
  expect_gt(nrow(res$cells), h - 2L)
})


test_that("the no-qhull fallback triangulates the points, not the hull", {
  # The regression: the fallback used to call st_triangulate() on the convex
  # hull POLYGON, throwing away every interior point.  It now triangulates the
  # MULTIPOINT, so both paths agree on the triangle count.
  skip_if_not_installed("geometry")
  pts <- .bt_points(30)
  qhull <- suppressMessages(build_tessellation(pts, method = "triangles",
                                               quiet = TRUE))

  lines <- capture_spatialkit_log({
    local_mocked_bindings(delaunayn = function(...) stop("qhull unavailable"),
                          .package = "geometry")
    geos <- suppressMessages(build_tessellation(pts, method = "triangles",
                                                quiet = TRUE))
  })

  expect_true(log_has(lines, "Falling back to GEOS"))
  expect_equal(nrow(geos$cells), nrow(qhull$cells))
  expect_equal(sum(as.numeric(sf::st_area(geos$cells))),
               sum(as.numeric(sf::st_area(qhull$cells))), tolerance = 1e-6)
  expect_length(geos$index, nrow(pts))
  expect_false(anyNA(geos$index))
})


test_that("build_tessellation needs at least 3 unique points for triangles", {
  pts <- .bt_points(30)

  expect_error(build_tessellation(pts[1:2, ], method = "triangles",
                                  quiet = TRUE),
               "need at least 3 unique points")
  # Three COINCIDENT points are two too few once duplicates are dropped.
  coincident <- sf::st_as_sf(data.frame(x = rep(5, 3), y = rep(5, 3)),
                             coords = c("x", "y"), crs = 32632)
  expect_error(build_tessellation(coincident, method = "triangles",
                                  quiet = TRUE),
               "need at least 3 unique points")
  # ... unless duplicates are kept on purpose, which is what the flag is for.
  expect_error(build_tessellation(pts[1:2, ], method = "triangles",
                                  keep_duplicates = TRUE, quiet = TRUE),
               "need at least 3 unique points")

  # The guard is triangles-only: Voronoi copes with two points.
  two <- suppressMessages(build_tessellation(pts[1:2, ], method = "voronoi",
                                             quiet = TRUE))
  expect_equal(nrow(two$cells), 2L)
})


test_that("build_tessellation tessellates CRS-less input", {
  # A CRS-less layer yields NA_crs_, which is a LIST rather than NULL, so it
  # was not recognised as "no CRS supplied" and create_voronoi_polygons() tried
  # to st_transform() a CRS-less object.
  set.seed(4)
  nocrs <- sf::st_as_sf(data.frame(x = runif(25, 0, 100), y = runif(25, 0, 100)),
                        coords = c("x", "y"))
  expect_true(is.na(sf::st_crs(nocrs)))

  res <- suppressMessages(build_tessellation(nocrs, method = "voronoi",
                                             quiet = TRUE))
  expect_equal(nrow(res$cells), 25L)
  expect_length(res$index, 25L)
  expect_false(anyNA(res$index))
  expect_true(is.na(sf::st_crs(res$cells)))     # nothing invented

  # And the triangle path, which shares the same CRS normalisation.
  skip_if_not_installed("geometry")
  tri <- suppressMessages(build_tessellation(nocrs, method = "triangles",
                                             quiet = TRUE))
  expect_gt(nrow(tri$cells), 0L)
  expect_false(anyNA(tri$index))
})


test_that("an explicit crs is applied to grids without a mismatch error", {
  # `crs` used to be dropped on the hex/square path, so create_grid_polygons()
  # re-projected the boundary on its own and the st_intersects() index was then
  # computed between two different coordinate systems.
  pts <- .bt_points()
  bnd <- .bt_boundary()

  for (m in c("hex", "square")) {
    res <- suppressMessages(build_tessellation(pts, boundary = bnd, method = m,
                                               approx_n_cells = 9, crs = 3857,
                                               quiet = TRUE))
    expect_equal(sf::st_crs(res$cells), sf::st_crs(3857), info = m)
    expect_gt(nrow(res$cells), 0L)
    # Points and cells ended up in the same CRS, so every point is placed.
    expect_false(anyNA(res$index), info = m)
    expect_length(res$index, nrow(pts))
  }

  # lon/lat input is projected automatically, and the grid comes back in the
  # projected CRS rather than in degrees.
  ll <- sf::st_transform(pts, 4326)
  bnd_ll <- sf::st_transform(bnd, 4326)
  auto <- suppressMessages(build_tessellation(ll, boundary = bnd_ll,
                                              method = "square",
                                              approx_n_cells = 9, quiet = TRUE))
  expect_false(sf::st_is_longlat(auto$cells))
  expect_false(anyNA(auto$index))
})


test_that("params echoes the arguments the tessellation was built with", {
  pts <- .bt_points()
  bnd <- .bt_boundary()

  # `expand` used to be hard-coded to 0 in params for the grid and triangle
  # methods, so a caller could not tell what had actually been asked for.
  for (m in .bt_all_methods) {
    if (m == "triangles") skip_if_not_installed("geometry")
    res <- suppressMessages(build_tessellation(pts, boundary = bnd, method = m,
                                               approx_n_cells = 9, expand = 25,
                                               quiet = TRUE))
    expect_equal(res$params$expand, 25, info = m)
    expect_equal(res$params$clip, TRUE, info = m)
    expect_equal(res$params$keep_duplicates, FALSE, info = m)
  }

  # Only the Voronoi method acts on `expand` -- documented, and worth pinning
  # so the echoed value is not mistaken for an applied one.
  no_exp <- suppressMessages(build_tessellation(pts, boundary = bnd,
                                                method = "voronoi", quiet = TRUE))
  with_exp <- suppressMessages(build_tessellation(pts, boundary = bnd,
                                                  method = "voronoi",
                                                  expand = 25, quiet = TRUE))
  expect_gt(sum(as.numeric(sf::st_area(with_exp$cells))),
            sum(as.numeric(sf::st_area(no_exp$cells))))

  sq_no  <- suppressMessages(build_tessellation(pts, boundary = bnd,
                                                method = "square",
                                                approx_n_cells = 9, quiet = TRUE))
  sq_exp <- suppressMessages(build_tessellation(pts, boundary = bnd,
                                                method = "square",
                                                approx_n_cells = 9, expand = 25,
                                                quiet = TRUE))
  expect_equal(sum(as.numeric(sf::st_area(sq_exp$cells))),
               sum(as.numeric(sf::st_area(sq_no$cells))), tolerance = 1e-6)
})


test_that("build_tessellation validates its inputs", {
  pts <- .bt_points()
  bnd <- .bt_boundary()

  expect_error(build_tessellation(pts, method = "hex", approx_n_cells = 9,
                                  quiet = TRUE),
               "`boundary` is required for hex/square grids")
  expect_error(build_tessellation(pts, method = "square", approx_n_cells = 9,
                                  quiet = TRUE),
               "`boundary` is required for hex/square grids")
  expect_error(build_tessellation(pts, boundary = pts, method = "voronoi",
                                  quiet = TRUE),
               "`boundary` must be polygonal")
  expect_error(build_tessellation(bnd, method = "voronoi", quiet = TRUE),
               "geometry must be one of: POINT, MULTIPOINT")
  expect_error(build_tessellation(sf::st_drop_geometry(pts), method = "voronoi",
                                  quiet = TRUE),
               "Expected an sf object")
  expect_error(build_tessellation(pts, method = "nonsense", quiet = TRUE),
               "'arg' should be one of")
})


test_that("clip = FALSE leaves cells extending past the boundary", {
  pts <- .bt_points()
  # A boundary covering only the lower-left quadrant, so clipping visibly bites.
  small <- sf::st_sf(geometry = sf::st_sfc(.bt_sq(0, 0, 50, 50), crs = 32632))

  clipped <- suppressMessages(build_tessellation(pts, boundary = small,
                                                 method = "square",
                                                 approx_n_cells = 16,
                                                 quiet = TRUE))
  whole <- suppressMessages(build_tessellation(pts, boundary = small,
                                               method = "square",
                                               approx_n_cells = 16, clip = FALSE,
                                               quiet = TRUE))

  expect_equal(sum(as.numeric(sf::st_area(clipped$cells))),
               as.numeric(sf::st_area(small)), tolerance = 1e-6)
  expect_gte(sum(as.numeric(sf::st_area(whole$cells))),
             sum(as.numeric(sf::st_area(clipped$cells))))
  expect_equal(clipped$params$clip, TRUE)
  expect_equal(whole$params$clip, FALSE)
})


# ---------------------------------------------------------------------------
# create_grid_polygons(): a caller-supplied `cellsize` wins over `n`
# ---------------------------------------------------------------------------

test_that("create_grid_polygons ignores `n` when the caller fixes `cellsize`", {
  # sf::st_make_grid() does NOT ignore `n` when `cellsize` is given: it uses
  # cellsize for the cell dimensions and nx/ny for the COUNTS, anchored at the
  # bounding-box corner.  create_grid_polygons(cellsize = 25, n = 2) on a
  # 100x100 boundary therefore returned 4 cells covering the bbox 0,0,50,50 --
  # a quarter of the study area -- and because clip = TRUE discards nothing
  # there, it looked like a perfectly ordinary grid.
  bnd  <- .bt_boundary()
  side <- 100
  expect_equal(as.numeric(sf::st_area(bnd)), side^2)

  lines <- capture_spatialkit_log(g <- create_grid_polygons(bnd, cellsize = 25,
                                                            n = 2))
  # `n` is dropped, and loudly: silently ignoring it would be its own trap.
  expect_true(log_has(lines, "both `cellsize` and `n` were supplied"))
  expect_true(log_has(lines, "`cellsize` wins and `n` \\(2 x 2\\) is ignored"))

  # 100 / 25 = 4 columns and 4 rows ...
  expect_equal(nrow(g), 16L)
  # ... covering the boundary completely, not a corner of it.
  expect_equal(sum(as.numeric(sf::st_area(g))), side^2, tolerance = 1e-8)
  expect_equal(unname(sf::st_bbox(g)[c("xmin", "ymin", "xmax", "ymax")]),
               c(0, 0, side, side), tolerance = 1e-8)
  # Every cell is the size that was asked for.
  expect_equal(unique(round(as.numeric(sf::st_area(g)), 8)), 25^2)

  # Identical to passing `cellsize` alone -- which is the point: the extra `n`
  # changed nothing.
  alone <- create_grid_polygons(bnd, cellsize = 25)
  expect_equal(nrow(alone), nrow(g))
  expect_equal(sf::st_bbox(alone), sf::st_bbox(g))

  # And no warning when only one of them is given.
  quiet <- capture_spatialkit_log(create_grid_polygons(bnd, cellsize = 25))
  expect_false(log_has(quiet, "both `cellsize` and `n`"))
  quiet_n <- capture_spatialkit_log(create_grid_polygons(bnd, n = 4))
  expect_false(log_has(quiet_n, "both `cellsize` and `n`"))
})

test_that("create_grid_polygons still honours a package-derived cell size", {
  # The `n` argument is forwarded to st_make_grid() when the PACKAGE derived
  # cellsize from it, because ceiling(w / (w / n)) floating-point-rounds one
  # cell too far (100 / (100/9) = 9.0000...4 -> 10 columns).  Dropping `n`
  # unconditionally would have reintroduced that off-by-one.
  bnd  <- .bt_boundary()
  side <- 100

  for (k in c(3L, 4L, 7L, 9L)) {
    g <- create_grid_polygons(bnd, n = k)
    expect_equal(nrow(g), k^2, info = paste("n =", k))
    expect_equal(sum(as.numeric(sf::st_area(g))), side^2, tolerance = 1e-8,
                 info = paste("n =", k))
    expect_equal(unique(round(as.numeric(sf::st_area(g)), 6)), (side / k)^2,
                 info = paste("n =", k))
  }

  # target_cells derives both, and covers the boundary in full.
  for (tc in c(16L, 25L, 100L)) {
    g <- create_grid_polygons(bnd, target_cells = tc)
    expect_equal(nrow(g), tc, info = paste("target_cells =", tc))
    expect_equal(sum(as.numeric(sf::st_area(g))), side^2, tolerance = 1e-8,
                 info = paste("target_cells =", tc))
  }
})
