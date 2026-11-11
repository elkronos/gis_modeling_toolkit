# tests/testthat/test-fix-movers-unplaced.R
# ---------------------------------------------------------------------------
# BUG-203.  GDAL returns a point it cannot move as an empty point and drops
# such a vertex from a line or a ring.  ensure_projected() said so (class
# spatialkit_rows_unplaced) for its moves out of lon/lat only; the other
# movers of a caller's layer returned the row emptied or short of a vertex
# without a word: harmonize_crs(), .transform_or_stamp(), .align_crs(),
# .cover_follow(), the moves of assign_features_to_polygons(), the free
# cells of cell_values_at(), and ensure_projected(target_crs =) on a
# projected layer.  Each now compares the rows before and after the move
# and warns in the words of the function called and the argument passed.
# ---------------------------------------------------------------------------

# Row 2 runs from longitude 10 to 170 along the equator: UTM zone 46
# (central meridian 93) places 90 and 170 but not 10.
.mvu_lines <- function(crs = 4326) {
  x <- sf::st_sf(v = 1:2, geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(80, 0), c(100, 0))),
    sf::st_linestring(rbind(c(10, 0), c(90, 0), c(170, 0))), crs = 4326))
  if (identical(crs, 4326)) x else sf::st_transform(x, crs)
}
.mvu_points <- function()
  sf::st_sf(v = 1:2, geometry = sf::st_sfc(sf::st_point(c(91, 0.1)), sf::st_point(c(10, 0.1)),
                                           crs = 4326))
.mvu_utm46 <- function()
  sf::st_sf(v = 1, geometry = sf::st_sfc(sf::st_point(c(5e5, 0)), crs = 32646))

# The warnings of class spatialkit_rows_unplaced an expression raises, and
# its value; every other warning and message is muffled.
.mvu_catch <- function(expr) {
  ws <- list()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    if (inherits(w, "spatialkit_rows_unplaced")) ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}

test_that("harmonize_crs() warns about a row the target CRS cannot place in full", {
  a <- .mvu_lines(); b <- .mvu_utm46()
  for (how in c("prefer", "target")) {
    r <- .mvu_catch(if (how == "prefer") harmonize_crs(a, b, prefer = "b")
                    else harmonize_crs(a, b, target_crs = 32646))
    expect_length(r$warnings, 1L)
    w <- r$warnings[[1L]]
    expect_identical(w$rows, 2L, info = how)
    expect_identical(w$n, 2L, info = how)
    expect_match(conditionMessage(w), "^harmonize_crs\\(\\): the move of `a` to EPSG:32646",
                 info = how)
    expect_match(conditionMessage(w), "row 2 of its 2, which has a location in longitude and latitude",
                 fixed = TRUE, info = how)
    # The geometry is the one the move returns, as before: row 2 without
    # its vertex at longitude 10.
    expect_identical(sf::st_geometry(r$value$a), sf::st_geometry(sf::st_transform(a, 32646)),
                     info = how)
    expect_identical(nrow(sf::st_coordinates(r$value$a[2, ])), 2L, info = how)
    expect_identical(r$value$b, b, info = how)
  }
  # The point the zone cannot place comes back empty, and is named.
  p <- .mvu_points()
  r <- .mvu_catch(harmonize_crs(p, b, prefer = "b"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  expect_identical(which(sf::st_is_empty(r$value$a)), 2L)
  # `b` moved into `a`'s CRS is named as `b`.
  r <- .mvu_catch(harmonize_crs(b, p, prefer = "a"))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), "the move of `b` to EPSG:32646", fixed = TRUE)
})

test_that("a ring with a vertex at latitude 95 is named by harmonize_crs()", {
  r95 <- sf::st_sf(v = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(90, 0), c(91, 0), c(91, 95), c(90, 1), c(90, 0)))), crs = 4326))
  r <- .mvu_catch(harmonize_crs(r95, .mvu_utm46(), prefer = "b"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 1L)
  expect_identical(nrow(sf::st_coordinates(r$value$a)), 4L)
})

test_that("from a projected CRS the warning names that CRS, and ensure_projected(target_crs =) gives it", {
  a <- .mvu_lines(3857)
  r <- .mvu_catch(harmonize_crs(a, .mvu_utm46(), prefer = "b"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  m <- conditionMessage(r$warnings[[1L]])
  expect_match(m, "which has a location in EPSG:3857", fixed = TRUE)
  expect_no_match(m, "longitude and latitude", fixed = TRUE)
  e <- .mvu_catch(ensure_projected(a, target_crs = 32646))
  expect_length(e$warnings, 1L)
  expect_identical(e$warnings[[1L]]$rows, 2L)
  m <- conditionMessage(e$warnings[[1L]])
  expect_match(m, "^ensure_projected\\(\\): the move of `x` to EPSG:32646")
  expect_match(m, "which has a location in EPSG:3857", fixed = TRUE)
  expect_identical(nrow(sf::st_coordinates(e$value[2, ])), 2L)
})

test_that(".transform_or_stamp() and .align_crs() warn in the words of their caller", {
  a <- .mvu_lines()
  r <- .mvu_catch(.transform_or_stamp(a, sf::st_crs(32646), what = "boundary",
                                      caller = "predict_surface"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^predict_surface\\(\\): the move of `boundary` to EPSG:32646")
  expect_match(conditionMessage(r$warnings[[1L]]), "Correct or drop those rows.$")
  # The same layer without a CRS, read as lon/lat and moved.
  r <- .mvu_catch(.transform_or_stamp(sf::st_set_crs(a, NA), sf::st_crs(32646),
                                      what = "points_sf", caller = "build_tessellation"))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^build_tessellation\\(\\): the move of `points_sf` to EPSG:32646")
  expect_match(conditionMessage(r$warnings[[1L]]), "a location in longitude and latitude",
               fixed = TRUE)
  r <- .mvu_catch(.align_crs(a, .mvu_utm46(), caller = "kriging_adequacy", arg = "cells_sf"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^kriging_adequacy\\(\\): the move of `cells_sf` to EPSG:32646")
  expect_identical(sf::st_geometry(r$value), sf::st_geometry(sf::st_transform(a, 32646)))
})

test_that("the builders' cover is named when the working CRS cannot place it", {
  cov <- sf::st_sfc(sf::st_point(c(90, 0)), sf::st_point(c(10, 0.1)), crs = 4326)
  pts <- .mvu_utm46()
  r <- .mvu_catch(.cover_follow(cov, pts, pts, caller = "build_tessellation"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^build_tessellation\\(\\): the move of `cover` to EPSG:32646")
  expect_identical(which(sf::st_is_empty(r$value)), 2L)
})

test_that("a public reader names its boundary, its seeds and its cover once", {
  set.seed(2)
  pts <- sf::st_transform(sf::st_as_sf(data.frame(x = stats::runif(30, 85, 95),
                                                  y = stats::runif(30, -4, 4)),
                                       coords = c("x", "y"), crs = 4326), 32646)
  # A ring with one vertex at longitude 10, which UTM zone 46 cannot place.
  bnd <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(80, -5), c(100, -5), c(100, 5), c(80, 5), c(10, 0), c(80, -5)))), crs = 4326))
  for (f in c("clip_target_for", "create_voronoi_polygons", "build_tessellation")) {
    r <- .mvu_catch(get(f)(pts, bnd))
    expect_length(r$warnings, 1L)
    expect_match(conditionMessage(r$warnings[[1L]]),
                 sprintf("^%s\\(\\): the move of `boundary` to EPSG:32646", f), info = f)
    expect_identical(r$warnings[[1L]]$rows, 1L, info = f)
  }
  seeds <- rbind(sf::st_transform(pts[1:5, ], 4326),
                 sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(10, 0.1)), crs = 4326)))
  r <- .mvu_catch(get_voronoi_seeds(sf::st_as_sf(sf::st_buffer(pts[1, ], 1e5)),
                                    method = "provided", seeds = seeds))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 6L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^get_voronoi_seeds\\(\\): the move of `seeds` to EPSG:32646")
  cov <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(90, 0)), sf::st_point(c(10, 0.1)),
                                         crs = 4326))
  r <- .mvu_catch(build_tessellation(pts, cover = cov))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^build_tessellation\\(\\): the move of `cover` to EPSG:32646")
})

test_that("assign_features_to_polygons() names its features and polygons, not harmonize_crs()'s `a`", {
  cells <- sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(3e5, -1e5), c(5e5, -1e5), c(5e5, 1e5), c(3e5, 1e5), c(3e5, -1e5)))),
    sf::st_polygon(list(rbind(c(5e5, -1e5), c(7e5, -1e5), c(7e5, 1e5), c(5e5, 1e5), c(5e5, -1e5)))),
    crs = 32646))
  f <- sf::st_sf(v = 1:3, geometry = sf::st_sfc(sf::st_point(c(92, 0.1)), sf::st_point(c(10, 0.1)),
                                                sf::st_point(c(94, 0.1)), crs = 4326))
  r <- .mvu_catch(assign_features_to_polygons(f, cells, keep_unassigned = TRUE, quiet = TRUE))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  m <- conditionMessage(r$warnings[[1L]])
  expect_match(m, "^assign_features_to_polygons\\(\\): the move of `features_sf` to EPSG:32646")
  expect_no_match(m, "harmonize_crs", fixed = TRUE)
  expect_identical(r$value$poly_id, c(1L, NA, 2L))
})

test_that("the scored join moves its features with the same warning", {
  # Cells in lon/lat whose record names a projected CRS are joined there; the
  # record is replaced by one that names UTM zone 46.
  local_mocked_bindings(.scored_crs = function(rule) sf::st_crs(32646))
  cells <- sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(90, -1), c(93, -1), c(93, 1), c(90, 1), c(90, -1)))),
    sf::st_polygon(list(rbind(c(93, -1), c(96, -1), c(96, 1), c(93, 1), c(93, -1)))),
    crs = 4326))
  f <- sf::st_sf(v = 1:3, geometry = sf::st_sfc(sf::st_point(c(92, 0.1)), sf::st_point(c(10, 0.1)),
                                                sf::st_point(c(94, 0.1)), crs = 4326))
  r <- .mvu_catch(assign_features_to_polygons(f, cells, keep_unassigned = TRUE, quiet = TRUE))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^assign_features_to_polygons\\(\\): the move of `features_sf` to EPSG:32646")
  expect_identical(r$value$poly_id, c(1L, NA, 2L))
})

test_that("cell_values_at() has the assignment's warning name `newdata` and `values`", {
  # cell_values_at() hands assign_features_to_polygons() the rows of
  # `newdata` that have a location and the polygons of `values` that carry an
  # ID; the warning names those arguments and their rows.
  f <- sf::st_sf(v = 1:3, geometry = sf::st_sfc(sf::st_point(c(92, 0.1)), sf::st_point(c(10, 0.1)),
                                                sf::st_point(c(94, 0.1)), crs = 4326))
  cells <- sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(3e5, -1e5), c(5e5, -1e5), c(5e5, 1e5), c(3e5, 1e5), c(3e5, -1e5)))),
    sf::st_polygon(list(rbind(c(5e5, -1e5), c(7e5, -1e5), c(7e5, 1e5), c(5e5, 1e5), c(5e5, -1e5)))),
    crs = 32646))
  close <- .unplaced_in_words_of("cell_values_at", list(
    features_sf = list(arg = "newdata", at = c(1L, 4L, 6L), n = 6L),
    polygons_sf = list(arg = "values", at = 2:3, n = 3L)))
  r <- tryCatch(.mvu_catch(assign_features_to_polygons(f, cells, keep_unassigned = TRUE,
                                                       quiet = TRUE)),
                finally = close())
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 4L)
  expect_identical(r$warnings[[1L]]$n, 6L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^cell_values_at\\(\\): the move of `newdata` to EPSG:32646 cannot place every vertex of row 4 of its 6")
  # The polygons, moved into the CRS of projected features.
  cells_ll <- sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(90, -1), c(93, -1), c(93, 1), c(90, 1), c(90, -1)))),
    sf::st_polygon(list(rbind(c(93, -1), c(96, -1), c(96, 1), c(10, 0), c(93, -1)))),
    crs = 4326))
  f46 <- sf::st_transform(f[-2, ], 32646)
  close <- .unplaced_in_words_of("cell_values_at", list(
    features_sf = list(arg = "newdata", at = 1:2, n = 2L),
    polygons_sf = list(arg = "values", at = 2:3, n = 3L)))
  r <- tryCatch(.mvu_catch(assign_features_to_polygons(f46, cells_ll, keep_unassigned = TRUE,
                                                       quiet = TRUE)),
                finally = close())
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 3L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^cell_values_at\\(\\): the move of `values` to EPSG:32646 cannot place every vertex of row 3 of its 3")
  # Outside the span the function speaks for itself again.
  r <- .mvu_catch(assign_features_to_polygons(f46, cells_ll, keep_unassigned = TRUE, quiet = TRUE))
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^assign_features_to_polygons\\(\\): the move of `polygons_sf`.* row 2 of its 2")
})

test_that("cell_values_at()'s free cells name the row of `values` the move cannot place", {
  values <- sf::st_sf(poly_id = c("a", "b"), geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(90, -1), c(92, -1), c(92, 1), c(90, 1), c(90, -1)))),
    sf::st_polygon(list(rbind(c(93, -1), c(96, -1), c(96, 1), c(10, 0), c(93, -1)))),
    crs = 4326))
  seed_xy <- matrix(c(-5e6, 0), 1L)
  r <- .mvu_catch(.values_at_free_cells(NA_character_, seed_xy, values, "poly_id",
                                        sf::st_crs(32646), sf::NA_crs_))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  expect_identical(r$warnings[[1L]]$n, 2L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^cell_values_at\\(\\): the move of `values` to EPSG:32646")
  # A row of `values` without an ID is no cell and is not moved; the row the
  # move cannot place is still named by its own row number, of all rows.
  values3 <- sf::st_sf(poly_id = c(NA, "a", "b"), geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(80, -1), c(82, -1), c(82, 1), c(80, 1), c(80, -1)))),
    sf::st_geometry(values)[[1L]], sf::st_geometry(values)[[2L]], crs = 4326))
  r <- .mvu_catch(.values_at_free_cells(NA_character_, seed_xy, values3, "poly_id",
                                        sf::st_crs(32646), sf::NA_crs_))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 3L)
  expect_identical(r$warnings[[1L]]$n, 3L)
  expect_match(conditionMessage(r$warnings[[1L]]), "row 3 of its 3", fixed = TRUE)
})

test_that("cell_values_at() names the rows of `values` and `newdata` it was given, end to end", {
  # The rows the assignment moves are some of the caller's rows: the
  # polygons of `values` that carry an ID, and the rows of `newdata` with a
  # location.  The warning names each by its row in the caller's layer.
  set.seed(3)
  pts <- sf::st_as_sf(data.frame(x = stats::runif(40, 4.8e5, 5.2e5),
                                 y = stats::runif(40, 1e4, 5e4)),
                      coords = c("x", "y"), crs = 32617)
  pts$y <- stats::rnorm(40)
  cells <- sf::st_sf(poly_id = c("a", "b"), geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(4.7e5, 0), c(5e5, 0), c(5e5, 6e4), c(4.7e5, 6e4), c(4.7e5, 0)))),
    sf::st_polygon(list(rbind(c(5e5, 0), c(5.3e5, 0), c(5.3e5, 6e4), c(5e5, 6e4), c(5e5, 0)))),
    crs = 32617))
  asg  <- suppressMessages(assign_features_to_polygons(pts, cells, quiet = TRUE))
  vals <- suppressWarnings(suppressMessages(
    summarize_by_cell(asg, "y", cells_sf = cells, estimate = "mean")))
  # `values` in lon/lat with a row without an ID first and a vertex at
  # latitude 95 in its last row: the polygons the assignment moves are rows
  # 2 and 3, and the one it cannot place in full is row 3.
  v3 <- rbind(vals[1, ], vals)
  v3$poly_id[1L] <- NA
  vll <- sf::st_transform(v3, 4326)
  g <- sf::st_geometry(vll)
  ring <- g[[3L]][[1L]]
  g[[3L]][[1L]] <- rbind(ring[1:2, ], c(ring[2L, 1L], 95), ring[-(1:2), ])
  vll <- sf::st_set_geometry(vll, g)
  r <- .mvu_catch(cell_values_at(pts[1:4, ], vll, asg, "y", id_col = "poly_id"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 3L)
  expect_identical(r$warnings[[1L]]$n, 3L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^cell_values_at\\(\\): the move of `values` to EPSG:32617 cannot place every vertex of row 3 of its 3")
  # `newdata` with a row without a location first, so the rows the
  # assignment moves are rows 2 and 3, and `values` in UTM zone 16, which
  # does not place the point (-2.5, 0.1) that zone 17 places: row 3.
  far <- sf::st_transform(sf::st_sfc(sf::st_point(c(-2.5, 0.1)), crs = 4326), 32617)
  nd <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(), sf::st_geometry(pts)[[1L]], far[[1L]],
                                        crs = 32617))
  r <- .mvu_catch(cell_values_at(nd, sf::st_transform(vals, 32616), asg, "y",
                                 id_col = "poly_id"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 3L)
  expect_identical(r$warnings[[1L]]$n, 3L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "^cell_values_at\\(\\): the move of `newdata` to EPSG:32616 cannot place every vertex of row 3 of its 3")
  # The row comes back as before: filled from the nearest scored row.
  expect_identical(r$value$filled, c(NA, FALSE, TRUE))
})

test_that("cell_values_at() reports seeds without a CRS in its own words, without the warning", {
  pts <- .mvu_utm46()
  seeds <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(sf::st_point(c(91.5, 0.1)),
                                                     sf::st_point(c(10.5, 0.1))))
  r <- .mvu_catch(.values_at_place(seeds, "seeds", "seeds", sf::st_crs(32646),
                                   sf::st_crs(32646), FALSE))
  expect_length(r$warnings, 0L)
  expect_identical(which(sf::st_is_empty(r$value)), 2L)
})

test_that("layers the target places in full raise no new warning and move as before", {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  set.seed(1)
  st39 <- sf::st_as_sf(data.frame(x = stats::runif(39, -120, -75), y = stats::runif(39, 30, 47)),
                       coords = c("x", "y"), crs = 4326)
  for (case in list(list(x = nc, to = 32119), list(x = st39, to = 5070),
                    list(x = sf::st_transform(st39, 3857), to = 5070))) {
    to <- sf::st_crs(case$to)
    ref <- sf::st_transform(case$x, to)
    tgt <- sf::st_sf(v = 1, geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = to))
    r <- .mvu_catch(harmonize_crs(case$x, tgt, prefer = "b"))
    expect_length(r$warnings, 0L)
    expect_identical(r$value$a, ref)
    r <- .mvu_catch(harmonize_crs(case$x, tgt, target_crs = to))
    expect_length(r$warnings, 0L)
    expect_identical(r$value$a, ref)
    r <- .mvu_catch(.transform_or_stamp(case$x, to, what = "x", caller = "f"))
    expect_length(r$warnings, 0L)
    expect_identical(r$value, ref)
    r <- .mvu_catch(.align_crs(case$x, tgt, caller = "f", arg = "x"))
    expect_length(r$warnings, 0L)
    expect_identical(r$value, ref)
    r <- .mvu_catch(ensure_projected(case$x, target_crs = to))
    expect_length(r$warnings, 0L)
    expect_identical(r$value, ref)
  }
  # An XYZ point among XY points is moved in XY, placed, without its Z: the
  # rows are compared on their vertices and their x and y.
  g <- sf::st_sfc(sf::st_point(c(-80, 40)), sf::st_point(c(-80.1, 40.1, 5)),
                  sf::st_point(c(-80.2, 40.2)), crs = 4326)
  r <- .mvu_catch(.transform_or_stamp(sf::st_sf(v = 1:3, geometry = g), sf::st_crs(32617),
                                      what = "x", caller = "f"))
  expect_length(r$warnings, 0L)
  expect_true(all(is.finite(sf::st_coordinates(r$value))))
  # So is an XYZ line among XY rows: its vertices are compared in x and y,
  # not with their Z.
  g <- sf::st_sfc(sf::st_point(c(-80, 40)),
                  sf::st_linestring(rbind(c(-80, 40, 1), c(-80.1, 40.1, 2))), crs = 4326)
  x <- sf::st_sf(v = 1:2, geometry = g)
  r <- .mvu_catch(.transform_or_stamp(x, sf::st_crs(32617), what = "x", caller = "f"))
  expect_length(r$warnings, 0L)
  expect_identical(r$value, sf::st_transform(x, 32617))
})

test_that("a layer without a CRS read as lon/lat is named in longitude and latitude", {
  # harmonize_crs() takes such a layer for EPSG:4326 and moves it; the
  # warning of assign_features_to_polygons() names the same CRS.
  cells <- sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(3e5, -1e5), c(5e5, -1e5), c(5e5, 1e5), c(3e5, 1e5), c(3e5, -1e5)))),
    sf::st_polygon(list(rbind(c(5e5, -1e5), c(7e5, -1e5), c(7e5, 1e5), c(5e5, 1e5), c(5e5, -1e5)))),
    crs = 32646))
  f <- sf::st_sf(v = 1:3, geometry = sf::st_sfc(sf::st_point(c(92, 0.1)), sf::st_point(c(10, 0.1)),
                                                sf::st_point(c(94, 0.1))))
  r <- .mvu_catch(assign_features_to_polygons(f, cells, keep_unassigned = TRUE, quiet = TRUE))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  m <- conditionMessage(r$warnings[[1L]])
  expect_match(m, "^assign_features_to_polygons\\(\\): the move of `features_sf` to EPSG:32646")
  expect_match(m, "row 2 of its 3, which has a location in longitude and latitude", fixed = TRUE)
  expect_identical(r$value$poly_id, c(1L, NA, 2L))
})

test_that("a collection the move empties is named, and a placed one is not", {
  gc <- function(x) sf::st_geometrycollection(list(sf::st_point(c(x, 0.1)),
                                                   sf::st_linestring(rbind(c(90, 0), c(91, 1)))))
  x <- sf::st_sfc(gc(92), gc(10), crs = 4326)
  r <- .mvu_catch(.align_crs(x, .mvu_utm46(), caller = "f", arg = "x"))
  expect_length(r$warnings, 1L)
  expect_identical(r$warnings[[1L]]$rows, 2L)
  expect_identical(.rows_unplaced(x[1], sf::st_transform(x[1], 32646)), integer(0))
})
