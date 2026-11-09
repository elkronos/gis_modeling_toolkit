# tests/testthat/test-fix-boundary-unlocated.R
# ---------------------------------------------------------------------------
# A `boundary` with no polygon that has a location.
#
# A row of `boundary` that is empty, or has a vertex with a coordinate that
# is not finite or lies beyond 1e150, covers no ground.  The readers of
# `boundary` never asked: clip_target_for() returned a target of area 0 or
# 5e204, create_voronoi_polygons() returned 0 cells and blamed the seeds,
# the seeders stopped with R's "no applicable method for 'st_geometry'" or
# put their seeds at y = 1e199, and build_tessellation(hex, square) refused
# a square beside an empty row as "not polygonal".  Such rows are now left
# out, read as given and again after each move, and a boundary with no row
# left is refused by name.
#
# The repair of a lon/lat boundary of rings with no height leaves no
# polygon, and the seeders stopped with R's "no applicable method".  They
# now refuse it by name.  A lon/lat ring of three vertices on one parallel is
# valid under s2 but its bounding box has no height, and the seeders drew in
# that box without end; they now refuse it at once.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
sq <- function(x0, y0, d)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d), c(x0, y0))))
with_s2_setting <- function(on, code) {
  old <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  suppressMessages(sf::sf_use_s2(on))
  code
}
err_of <- function(f) tryCatch({ quiet(f()); NULL }, error = function(e) e)

# 120 points in a 1e5 m square in EPSG:32617, and the square.
set.seed(5)
P <- sf::st_as_sf(data.frame(x = 5e5 + stats::runif(120, 0, 1e5),
                             y = 5e6 + stats::runif(120, 0, 1e5)),
                  coords = c("x", "y"), crs = 32617)
seeds <- P[1:10, ]
square <- sf::st_sfc(sq(5e5, 5e6, 1e5), crs = 32617)
# Rows without a location: empty, a vertex at Inf, at (x, 1e200) and at
# (-1e200, y); in the points' CRS and in EPSG:32616.
unlocated_rows <- function(crs) {
  x0 <- if (crs == 32617) 5e5 else 3.2e5
  list(
    empty  = sf::st_polygon(),
    inf    = sf::st_polygon(list(rbind(c(x0, 5e6), c(x0 + 1e4, 5e6), c(x0 + 1e4, Inf),
                                       c(x0, 5e6)))),
    far_y  = sf::st_polygon(list(rbind(c(x0, 5e6), c(x0 + 1e4, 5e6), c(x0 + 1e4, 1e200),
                                       c(x0, 5e6)))),
    far_x  = sf::st_polygon(list(rbind(c(x0, 5e6), c(-1e200, 5e6), c(x0, 5.01e6),
                                       c(x0, 5e6)))))
}
unlocated_boundaries <- function() {
  out <- list(multi_empty = sf::st_sf(geometry = sf::st_sfc(sf::st_multipolygon(),
                                                            crs = 32617)))
  for (crs in c(32617, 32616)) {
    rows <- unlocated_rows(crs)
    for (nm in names(rows))
      out[[sprintf("%s_%d", nm, crs)]] <- sf::st_sf(geometry = sf::st_sfc(rows[[nm]], crs = crs))
  }
  out
}
as_sf <- function(b) if (inherits(b, "sfc")) sf::st_sf(geometry = b) else b
calls_on <- function(b) list(
  clip      = function() clip_target_for(P, boundary = b, quiet = TRUE),
  clip_x    = function() clip_target_for(P, boundary = b, expand = 0.1, quiet = TRUE),
  voronoi   = function() create_voronoi_polygons(seeds, boundary = b, quiet = TRUE),
  bt_vor    = function() build_tessellation(P, boundary = b, method = "voronoi", quiet = TRUE),
  bt_hex    = function() build_tessellation(P, boundary = b, method = "hex",
                                            approx_n_cells = 20, quiet = TRUE),
  bt_sq     = function() build_tessellation(P, boundary = b, method = "square",
                                            approx_n_cells = 20, quiet = TRUE),
  bt_sq_crs = function() build_tessellation(P, boundary = b, method = "square",
                                            approx_n_cells = 20, crs = 32617, quiet = TRUE),
  grid_size = function() create_grid_polygons(b, cellsize = 2e4, quiet = TRUE),
  grid_n    = function() create_grid_polygons(b, target_cells = 20, quiet = TRUE),
  cached    = function() create_grid_polygons_cached(b, target_cells = 20, quiet = TRUE),
  # get_voronoi_seeds() takes `boundary` as an sf layer only.
  random    = function() get_voronoi_seeds(as_sf(b), n = 5, method = "random", set_seed = 1),
  kmeans    = function() get_voronoi_seeds(as_sf(b), n = 5, method = "kmeans", set_seed = 1),
  vsr       = function() voronoi_seeds_random(b, 5, set_seed = 1),
  folds     = function() make_folds(P, method = "block_kfold", k = 2, boundary = b,
                                    seed = 1))

test_that("a boundary with no polygon that has a location is refused by every reader", {
  for (nm in names(unlocated_boundaries())) {
    for (as_sfc in c(FALSE, TRUE)) {
      b <- unlocated_boundaries()[[nm]]
      if (as_sfc) b <- sf::st_geometry(b)
      calls <- calls_on(b)
      for (f in names(calls)) {
        lab <- sprintf("%s (%s) through %s", nm, if (as_sfc) "sfc" else "sf", f)
        e <- err_of(calls[[f]])
        expect_s3_class(e, "error")
        expect_match(conditionMessage(e),
                     "^[a-z_]+\\(\\): `boundary` has no polygon with a location \\(",
                     label = lab)
        expect_null(e$call, label = lab)
      }
    }
  }
})

test_that("a square beside a row without a location gives what the square alone gives", {
  for (crs in c(32617, 32616)) {
    rows <- unlocated_rows(crs)
    for (nm in names(rows)) {
      lab <- sprintf("%s in EPSG:%d", nm, crs)
      # The square in the points' CRS, the other row in its own CRS, moved
      # into the square's first where needed.
      other <- sf::st_sfc(rows[[nm]], crs = crs)
      if (crs != 32617) other <- spatialkit:::.transform_or_stamp(other, sf::st_crs(32617))
      with  <- sf::st_sf(geometry = c(square, other))
      alone <- sf::st_sf(geometry = square)
      a <- calls_on(alone)
      w <- calls_on(with)
      for (f in names(a))
        expect_identical(quiet(w[[f]]()), quiet(a[[f]]()), label = paste(lab, f))
    }
  }
})

test_that("a boundary whose only row the move empties is refused, naming the move", {
  # In Web Mercator, a vertex at x = 1e20 has a location as given (it is
  # finite and within 1e150), and the move into degrees cannot place it, so
  # the row comes back empty.
  far_x <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(1e6, 5e6), c(1e20, 5e6), c(1e6, 5.1e6), c(1e6, 5e6)))), crs = 3857))
  e <- err_of(function() create_grid_polygons(far_x, cellsize = 0.1, crs = 4326, quiet = TRUE))
  expect_identical(conditionMessage(e), paste0(
    "create_grid_polygons(): `boundary` has no polygon with a location (every row is empty ",
    "after its move to the CRS the call works in, which cannot place them), so no cell can ",
    "be built in it."))
})

# The flat lon/lat boundary: every vertex at latitude 60, so the ring has no
# height and s2 refuses it; its repair leaves an empty collection.
flat <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(10, 60), c(11, 60), c(11, 60), c(10, 60), c(10, 60)))), crs = 4326))
square_ll <- sf::st_sfc(sq(20, 50, 0.5), crs = 4326)

test_that("a lon/lat boundary its repair empties is refused by the seeders by name", {
  for (on in c(TRUE, FALSE)) with_s2_setting(on, {
    for (f in list(
      function() voronoi_seeds_random(flat, 5, set_seed = 1),
      function() get_voronoi_seeds(flat, n = 5, method = "random", set_seed = 1),
      function() get_voronoi_seeds(flat, n = 5, method = "kmeans", set_seed = 1))) {
      e <- err_of(f)
      expect_s3_class(e, "error")
      expect_match(conditionMessage(e), paste0(
        "`boundary` has no polygon with a location (the repair of the rows s2 refuses ",
        "leaves no polygon with an area)"), fixed = TRUE)
      expect_false(grepl("no applicable method", conditionMessage(e), fixed = TRUE))
    }
    beside <- sf::st_sf(geometry = c(sf::st_geometry(flat), square_ll))
    alone  <- sf::st_sf(geometry = square_ll)
    expect_identical(sf::st_geometry(quiet(voronoi_seeds_random(beside, 5, set_seed = 1))),
                     sf::st_geometry(quiet(voronoi_seeds_random(alone, 5, set_seed = 1))))
    expect_identical(
      sf::st_geometry(quiet(get_voronoi_seeds(beside, n = 5, method = "random", set_seed = 1))),
      sf::st_geometry(quiet(get_voronoi_seeds(alone, n = 5, method = "random", set_seed = 1))))
    expect_identical(
      sf::st_geometry(quiet(get_voronoi_seeds(beside, n = 5, method = "kmeans", set_seed = 1))),
      sf::st_geometry(quiet(get_voronoi_seeds(alone, n = 5, method = "kmeans", set_seed = 1))))
  })
})

# A ring of three vertices on one parallel: valid under s2 (its edges are
# great circles) and its bounding box has no height.
sliver <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(10, 60), c(11, 60), c(12, 60), c(10, 60)))), crs = 4326))

test_that("the box test refuses a region whose bounding box has no height or width", {
  # Checked on its own first: a regression fails here, before any call
  # below that could draw without end.
  e <- err_of(function() spatialkit:::.sample_box_check(sf::st_geometry(sliver), "probe",
                                                      "no seed can be drawn in it"))
  expect_identical(conditionMessage(e), paste0(
    "probe(): `boundary` cannot be sampled (its bounding box, read from its coordinates, ",
    "has no height), so no seed can be drawn in it."))
  tall <- sf::st_sfc(sf::st_polygon(list(rbind(c(10, 60), c(10, 61), c(10, 62), c(10, 60)))),
                     crs = 4326)
  expect_match(conditionMessage(err_of(function()
    spatialkit:::.sample_box_check(tall, "probe", "x"))), "has no width", fixed = TRUE)
  expect_null(spatialkit:::.sample_box_check(square_ll, "probe", "x"))
})

test_that("the seeders refuse a lon/lat sliver whose box has no height at once", {
  expect_true(spatialkit:::.with_s2(sf::st_is_valid(sliver)))
  for (f in list(
    function() voronoi_seeds_random(sliver, 5, set_seed = 1),
    function() get_voronoi_seeds(sliver, n = 5, method = "random", set_seed = 1),
    function() get_voronoi_seeds(sliver, n = 5, method = "kmeans", set_seed = 1))) {
    e <- err_of(f)
    expect_match(conditionMessage(e), paste0(
      "\\(\\): `boundary` cannot be sampled \\(its bounding box, read from its coordinates, ",
      "has no height\\), so no (seed|point of the k-means cloud) can be drawn in it\\.$"))
  }
  # A valid square beside a valid square samples as before.
  two <- sf::st_sf(geometry = c(square_ll, sf::st_sfc(sq(21, 50, 0.5), crs = 4326)))
  expect_identical(nrow(quiet(voronoi_seeds_random(two, 5, set_seed = 1))), 5L)
})

# A flat row: every vertex on one line (projected) or one parallel (in
# degrees), so the ring has no height.  GEOS and s2 refuse it, and their
# repair leaves a line or an empty collection.  The builders repaired
# `boundary` with no test of what the repair left: beside a square the hex
# and square grids refused it as "not polygonal", and alone the clip target
# was a line and the Voronoi cells held none of the seeds.
flat_p <- sf::st_sfc(sf::st_polygon(list(rbind(
  c(5e5, 5.05e6), c(5.5e5, 5.05e6), c(5.5e5, 5.05e6), c(5e5, 5.05e6), c(5e5, 5.05e6)))),
  crs = 32617)
# The message of an error, or a word that the call returned.
msg_of <- function(e) if (inherits(e, "condition")) conditionMessage(e) else "<returned>"
repair_refusal <- paste0(
  "^[a-z_]+\\(\\): `boundary` has no polygon with a location \\(the repair of the rows ",
  "(GEOS|s2) refuses leaves no polygon with an area\\), so ")

test_that("a projected row its repair leaves without a polygon is left out by every reader", {
  calls <- calls_on(sf::st_sf(geometry = flat_p))
  for (f in names(calls)) {
    e <- err_of(calls[[f]])
    expect_match(msg_of(e), repair_refusal, label = f)
  }
  a <- calls_on(sf::st_sf(geometry = square))
  w <- calls_on(sf::st_sf(geometry = c(square, flat_p)))
  for (f in names(a))
    expect_identical(quiet(w[[f]]()), quiet(a[[f]]()), label = f)
})

# The same in degrees, with points inside a 0.5 degree square.
set.seed(5)
P_ll <- sf::st_as_sf(data.frame(x = 20 + stats::runif(60, 0, 0.5),
                                y = 50 + stats::runif(60, 0, 0.5)),
                     coords = c("x", "y"), crs = 4326)
calls_ll <- function(b) list(
  clip      = function() clip_target_for(P_ll, boundary = b, quiet = TRUE),
  voronoi   = function() create_voronoi_polygons(P_ll[1:10, ], boundary = b, quiet = TRUE),
  bt_vor    = function() build_tessellation(P_ll, boundary = b, method = "voronoi", quiet = TRUE),
  bt_hex    = function() build_tessellation(P_ll, boundary = b, method = "hex",
                                            approx_n_cells = 20, quiet = TRUE),
  bt_sq     = function() build_tessellation(P_ll, boundary = b, method = "square",
                                            approx_n_cells = 20, quiet = TRUE),
  bt_sq_crs = function() build_tessellation(P_ll, boundary = b, method = "square",
                                            approx_n_cells = 20, crs = 4326, quiet = TRUE),
  grid_deg  = function() create_grid_polygons(b, cellsize = 0.1, crs = 4326, quiet = TRUE),
  grid_n    = function() create_grid_polygons(b, target_cells = 20, quiet = TRUE),
  cached    = function() create_grid_polygons_cached(b, target_cells = 20, quiet = TRUE),
  random    = function() get_voronoi_seeds(as_sf(b), n = 5, method = "random", set_seed = 1),
  kmeans    = function() get_voronoi_seeds(as_sf(b), n = 5, method = "kmeans", set_seed = 1),
  vsr       = function() voronoi_seeds_random(b, 5, set_seed = 1),
  folds     = function() make_folds(P_ll, method = "block_kfold", k = 2, boundary = b,
                                    seed = 1))

test_that("a lon/lat row its repair leaves without a polygon is left out by every reader", {
  for (on in c(TRUE, FALSE)) with_s2_setting(on, {
    calls <- calls_ll(flat)
    for (f in names(calls)) {
      lab <- sprintf("%s, s2 %s", f, on)
      expect_match(msg_of(err_of(calls[[f]])), repair_refusal, label = lab)
    }
    a <- calls_ll(sf::st_sf(geometry = square_ll))
    w <- calls_ll(sf::st_sf(geometry = c(square_ll, sf::st_geometry(flat))))
    for (f in names(a))
      expect_identical(quiet(w[[f]]()), quiet(a[[f]]()), label = sprintf("%s, s2 %s", f, on))
  })
})

test_that("build_tessellation() leaves out a row its own move empties", {
  # Points in degrees with `crs = 4326`: the cells are built in the CRS
  # picked for the points, and a boundary in Web Mercator is moved there.  A
  # vertex at x = 1e20 has a location as given, and the move cannot place
  # it, so the row comes back empty after the test of `boundary` as given.
  set.seed(5)
  pts <- sf::st_as_sf(data.frame(x = 10 + stats::runif(60, 0, 1),
                                 y = 50 + stats::runif(60, 0, 1)),
                      coords = c("x", "y"), crs = 4326)
  far <- sf::st_sfc(sf::st_polygon(list(rbind(c(1e6, 5e6), c(1e20, 5e6), c(1e6, 5.1e6),
                                              c(1e6, 5e6)))), crs = 3857)
  sqm <- sf::st_transform(sf::st_sfc(sq(10, 50, 1), crs = 4326), 3857)
  calls <- function(b) list(
    square  = function() build_tessellation(pts, boundary = b, method = "square",
                                            approx_n_cells = 20, crs = 4326, quiet = TRUE),
    voronoi = function() build_tessellation(pts, boundary = b, method = "voronoi",
                                            crs = 4326, quiet = TRUE))
  only <- calls(sf::st_sf(geometry = far))
  for (f in names(only))
    expect_identical(msg_of(err_of(only[[f]])), paste0(
      "build_tessellation(): `boundary` has no polygon with a location (every row is ",
      "empty after its move to the CRS the call works in, which cannot place them), so ",
      "no cell can be built in it."), label = f)
  a <- calls(sf::st_sf(geometry = sqm))
  w <- calls(sf::st_sf(geometry = c(sqm, far)))
  for (f in names(a))
    expect_identical(quiet(w[[f]]()), quiet(a[[f]]()), label = f)
})

test_that("a row without a location that also reaches beyond latitude 90 is left out", {
  # Its vertex at Inf leaves it without a location, so it is left out as
  # every such row is, and not refused for its latitude.
  both <- sf::st_polygon(list(rbind(c(20, 50), c(20.2, 50), c(20.2, Inf), c(20, 95),
                                    c(20, 50))))
  beside <- sf::st_sf(geometry = c(square_ll, sf::st_sfc(both, crs = 4326)))
  alone  <- sf::st_sf(geometry = square_ll)
  expect_identical(quiet(voronoi_seeds_random(beside, 5, set_seed = 1)),
                   quiet(voronoi_seeds_random(alone, 5, set_seed = 1)))
})

test_that("make_folds() refuses a boundary that is not an sf or sfc layer by name", {
  bb <- sf::st_bbox(square)
  for (b in list(bb, data.frame(x = 1:3, y = 1:3))) {
    e <- err_of(function() make_folds(P, method = "block_kfold", k = 2, boundary = b,
                                      seed = 1))
    expect_match(msg_of(e),
                 "^make_folds\\(\\): `boundary` must be an sf or sfc polygon layer")
  }
})

test_that("the builders' repair leaves out a row it empties, and keeps the rest", {
  # Rows that only the builders' own repair empties (after a move, or with
  # s2 off in degrees) meet .boundary_repaired() and not the test of
  # `boundary` as given.
  rep <- spatialkit:::.boundary_repaired
  alone <- sf::st_sf(geometry = square)
  expect_identical(rep(sf::st_sf(geometry = c(square, flat_p)), "probe", "x"),
                   spatialkit:::.safe_make_valid(alone))
  expect_identical(rep(c(square, flat_p), "probe", "x"),
                   spatialkit:::.safe_make_valid(square))
  expect_identical(rep(alone, "probe", "x"), spatialkit:::.safe_make_valid(alone))
  expect_identical(msg_of(err_of(function() rep(sf::st_sf(geometry = flat_p), "probe",
                                                "nothing can be done"))), paste0(
    "probe(): `boundary` has no polygon with a location (the repair of the rows GEOS ",
    "refuses leaves no polygon with an area), so nothing can be done."))
  # A row that was not a polygon as given is left to the callers' type tests.
  line <- sf::st_sfc(sf::st_linestring(rbind(c(5e5, 5e6), c(6e5, 5e6))), crs = 32617)
  expect_identical(nrow(rep(sf::st_sf(geometry = c(square, line)), "probe", "x")), 2L)
  # With s2 off a grid in degrees is clipped by GEOS, for which a ring of
  # three vertices on one parallel has no area.
  with_s2_setting(FALSE, {
    e <- err_of(function() create_grid_polygons(sliver, cellsize = 0.1, crs = 4326, quiet = TRUE))
    expect_match(msg_of(e), paste0(
      "create_grid_polygons(): `boundary` has no polygon with a location (the repair of the ",
      "rows GEOS refuses leaves no polygon with an area)"), fixed = TRUE)
  })
})
