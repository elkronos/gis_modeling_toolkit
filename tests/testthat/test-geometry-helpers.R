# ===========================================================================
# Small geometry utilities: input assertions, target clipping, grid cells.
# ===========================================================================

test_that(".assert_sf rejects non-sf objects", {
  expect_error(spatialkit:::.assert_sf(data.frame(x = 1)), "Expected an sf object")
})


test_that(".assert_sf requires EVERY geometry to match, not just one", {
  # `any()` let a mixed layer through a POINT-only check, contradicting the
  # error text and handing non-point geometry to code that assumes points.
  mixed <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)),
      sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))),
      crs = 32632
    )
  )
  expect_error(spatialkit:::.assert_sf(mixed, "POINT", "mixed"),
               "geometry must be one of: POINT")
  # The message names what it actually found, both types.
  expect_error(spatialkit:::.assert_sf(mixed, "POINT", "mixed"), "POLYGON")

  # A homogeneous layer of an accepted type still passes ...
  pts <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(
    sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = 32632))
  expect_silent(spatialkit:::.assert_sf(pts, "POINT", "pts"))

  # ... and a zero-row layer keeps failing: all() is vacuously TRUE on the
  # empty set, so the explicit length check is what holds this line.
  expect_error(spatialkit:::.assert_sf(pts[0, , drop = FALSE], "POINT", "empty"),
               "geometry must be one of")
})


# ---------------------------------------------------------------------------
# clip_target_for()
# ---------------------------------------------------------------------------

.gh_points <- function(crs = 32632) {
  sf::st_sf(
    id = 1:5,
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)), sf::st_point(c(1, 1)),
      sf::st_point(c(0.5, 0.5)),
      crs = crs
    )
  )
}

test_that("clip_target_for contains the points it was built from", {
  pts <- .gh_points()
  tgt <- clip_target_for(pts, quiet = TRUE)

  expect_s3_class(tgt, "sf")
  expect_true(all(sf::st_geometry_type(tgt) %in% c("POLYGON", "MULTIPOLYGON")))
  expect_equal(sf::st_crs(tgt), sf::st_crs(pts))

  # The point that matters: the polygon actually covers the input.  Without
  # this the test passed for any polygon anywhere on the planet.
  expect_true(all(lengths(sf::st_covered_by(pts, tgt)) > 0L))

  # With expand = 0 it is the bounding box exactly -- a unit square here.
  expect_equal(as.numeric(sf::st_area(tgt)), 1, tolerance = 1e-9)
})

test_that("clip_target_for expands by an absolute distance and by a fraction", {
  pts <- .gh_points()

  # Absolute: >= 1 is taken as a distance in CRS units, so the bbox grows by
  # exactly that much on every side.
  abs_tgt <- clip_target_for(pts, expand = 10, quiet = TRUE)
  bb <- sf::st_bbox(abs_tgt)
  expect_equal(unname(as.numeric(bb[c("xmin", "ymin", "xmax", "ymax")])),
               c(-10, -10, 11, 11), tolerance = 1e-6)

  # Fractional: 0 < expand < 1 is a fraction of the larger bbox side (1 here),
  # so expand = 0.5 buffers by 0.5.
  frac_tgt <- clip_target_for(pts, expand = 0.5, quiet = TRUE)
  bb_f <- sf::st_bbox(frac_tgt)
  expect_equal(unname(as.numeric(bb_f[c("xmin", "ymin", "xmax", "ymax")])),
               c(-0.5, -0.5, 1.5, 1.5), tolerance = 1e-6)

  # Expanding can only grow the target, never shrink it.
  expect_gt(as.numeric(sf::st_area(frac_tgt)), as.numeric(sf::st_area(
    clip_target_for(pts, quiet = TRUE))))
})

test_that("clip_target_for uses a supplied boundary in place of the bbox", {
  pts <- .gh_points()
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-5, -5), c(5, -5), c(5, 5), c(-5, 5), c(-5, -5)))), crs = 32632))

  tgt <- clip_target_for(pts, boundary = bnd, quiet = TRUE)
  expect_equal(as.numeric(sf::st_area(tgt)), 100, tolerance = 1e-9)
  expect_true(all(lengths(sf::st_covered_by(pts, tgt)) > 0L))

  # `expand` is then a fraction of the BOUNDARY's extent (10), not the points'.
  tgt_e <- clip_target_for(pts, boundary = bnd, expand = 0.1, quiet = TRUE)
  bb <- sf::st_bbox(tgt_e)
  expect_equal(unname(as.numeric(bb[c("xmin", "xmax")])), c(-6, 6),
               tolerance = 1e-6)

  expect_error(clip_target_for(pts, boundary = pts, quiet = TRUE),
               "must be polygonal")
})

test_that("clip_target_for projects lon/lat before applying `expand`", {
  # Regression: the fraction-of-extent form of `expand` is derived from the
  # bounding box, which for lon/lat is measured in DEGREES, while st_buffer()
  # reads `dist` as METRES -- five orders of magnitude apart.  A 0.2-degree
  # extent used to buffer by 0.02 METRES instead of ~2.2 km, i.e. not at all.
  ll <- sf::st_sf(
    id = 1:4,
    geometry = sf::st_sfc(
      sf::st_point(c(10, 50)),   sf::st_point(c(10.2, 50)),
      sf::st_point(c(10, 50.2)), sf::st_point(c(10.2, 50.2)),
      crs = 4326
    )
  )

  plain <- clip_target_for(ll, quiet = TRUE)
  grown <- clip_target_for(ll, expand = 0.1, quiet = TRUE)

  # The documented return: a local projected CRS, not the input's lon/lat.
  expect_false(sf::st_is_longlat(plain))
  expect_equal(sf::st_crs(grown), sf::st_crs(plain))

  bb_p <- sf::st_bbox(plain); bb_g <- sf::st_bbox(grown)
  side_p <- max(as.numeric(bb_p[["xmax"]] - bb_p[["xmin"]]),
                as.numeric(bb_p[["ymax"]] - bb_p[["ymin"]]))

  # expand = 0.1 buffers by 0.1 * (longer side), applied on both sides.
  expect_equal(as.numeric(bb_g[["xmax"]] - bb_g[["xmin"]]),
               as.numeric(bb_p[["xmax"]] - bb_p[["xmin"]]) + 0.2 * side_p,
               tolerance = 0.02 * side_p)
  # Degrees-as-metres would have grown this by 0.04 m out of ~15 km.
  expect_gt(as.numeric(bb_g[["xmax"]] - bb_g[["xmin"]]) -
              as.numeric(bb_p[["xmax"]] - bb_p[["xmin"]]), 1000)

  # It still covers the input, once the input is put in the same CRS.
  expect_true(all(lengths(sf::st_covered_by(
    sf::st_transform(ll, sf::st_crs(grown)), grown)) > 0L))

  # The CRS change is announced unless silenced.
  expect_message(clip_target_for(ll, expand = 0.1), "projecting to a local")
})

test_that("clip_target_for buffers a degenerate (collinear) point set", {
  # A zero-width bbox is not a polygon; the function must fall back to a
  # buffer rather than returning something empty.
  collinear <- sf::st_sf(
    id = 1:3,
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(0, 1)),
                          sf::st_point(c(0, 2)), crs = 32632)
  )
  tgt <- suppressMessages(clip_target_for(collinear))
  expect_gt(as.numeric(sf::st_area(tgt)), 0)
  expect_true(all(lengths(sf::st_covered_by(collinear, tgt)) > 0L))
})


# ---------------------------------------------------------------------------
# create_grid_polygons()
# ---------------------------------------------------------------------------

test_that("create_grid_polygons hits target_cells exactly on a square boundary", {
  bnd <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))),
      crs = 32632
    )
  )
  grid <- create_grid_polygons(bnd, target_cells = 9, type = "square",
                               quiet = TRUE)

  expect_s3_class(grid, "sf")
  # 9 over a square is a 3 x 3 grid with no partial cells, so the count is
  # exact -- `>= 4` also passed for 4, 900 and 90,000.
  expect_equal(nrow(grid), 9L)
  expect_true(all(as.character(sf::st_geometry_type(grid, by_geometry = TRUE)) ==
                    "POLYGON"))
  # poly_id is a dense 1..n key; downstream joins rely on it being unique.
  expect_equal(sort(grid$poly_id), seq_len(nrow(grid)))
  expect_false(anyDuplicated(grid$poly_id) > 0L)

  # The cells tile the boundary: same total area, equal-sized cells.
  expect_equal(sum(as.numeric(sf::st_area(grid))),
               as.numeric(sf::st_area(bnd)), tolerance = 1e-6)
  expect_equal(as.numeric(sf::st_area(grid)), rep(100 * 100 / 9, 9L),
               tolerance = 1e-6)
  expect_equal(sf::st_crs(grid), sf::st_crs(bnd))

  # Other perfect squares land exactly too.
  for (tc in c(4L, 16L, 25L))
    expect_equal(nrow(create_grid_polygons(bnd, target_cells = tc,
                                           type = "square", quiet = TRUE)), tc)
})

test_that("create_grid_polygons clips to an irregular boundary", {
  # A triangle covering half the square: about half the cells survive the clip,
  # and the surviving area equals the boundary's.
  tri <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(100, 0), c(0, 100), c(0, 0)))), crs = 32632))

  clipped <- create_grid_polygons(tri, target_cells = 100, type = "square",
                                  quiet = TRUE)
  unclipped <- create_grid_polygons(tri, target_cells = 100, type = "square",
                                    clip = FALSE, quiet = TRUE)

  expect_lt(nrow(clipped), nrow(unclipped))
  expect_gt(nrow(clipped), 0.4 * nrow(unclipped))     # ~half, plus the edge row
  expect_equal(sum(as.numeric(sf::st_area(clipped))),
               as.numeric(sf::st_area(tri)), tolerance = 1e-6)
  expect_equal(sort(clipped$poly_id), seq_len(nrow(clipped)))
})

test_that("create_grid_polygons validates its sizing arguments", {
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)))), crs = 32632))

  expect_error(create_grid_polygons(bnd, quiet = TRUE),
               "supply either 'cellsize', 'n', or a positive 'target_cells'")
  expect_error(create_grid_polygons(bnd, n = 0, quiet = TRUE),
               "'n' must be integer >= 1")
  # `n` is validated the same way whether or not cellsize is also given.
  expect_error(create_grid_polygons(bnd, cellsize = 1, n = -2, quiet = TRUE),
               "'n' must be integer >= 1")
  expect_error(create_grid_polygons(bnd, cellsize = 0, quiet = TRUE),
               "'cellsize' must be positive numeric")
  expect_error(create_grid_polygons(.gh_points(), target_cells = 4, quiet = TRUE),
               "must be polygonal")

  # n = c(2, 5) means 2 columns and 5 rows -- 10 cells, not 4 or 25.
  expect_equal(nrow(create_grid_polygons(bnd, n = c(2, 5), quiet = TRUE)), 10L)
  # cellsize divides the 10 x 10 extent into 5 x 5 = 25 cells.
  expect_equal(nrow(create_grid_polygons(bnd, cellsize = 2, quiet = TRUE)), 25L)
})
