# tests/testthat/test-fix-s2-boundaries.R
# ---------------------------------------------------------------------------
# Boundaries and grids in degrees.
#
# A grid laid in degrees (`cellsize` with a geographic `crs`) is clipped with
# s2, which writes some clipped cells with a corner listed twice in a row; s2
# then refuses those cells, so build_tessellation() stopped with s2's bare
# message on a plain square.  The cells are now returned without the repeat,
# and a cell s2 refuses even then stops the call with the classed error.
#
# voronoi_seeds_random(), the key of create_grid_polygons_cached() and
# make_folds(block_kfold) unioned `boundary` before repairing it, so a
# boundary the other builders repair stopped them with GEOS's or s2's bare
# message.  The rows GEOS or s2 refuses are now repaired before the union.
#
# In degrees sf::st_make_valid() leaves crossing edges as they are, so a
# lon/lat bowtie stopped get_voronoi_seeds() and the grid in degrees.  The
# rows s2 still refuses are now rebuilt with crossing edges split, and a row
# s2 refuses even then stops the call with the classed error.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
sq <- function(x0, y0, d)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d), c(x0, y0))))
sq_rep <- function(x0, y0, d)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0), c(x0 + d, y0 + d),
                            c(x0, y0 + d), c(x0, y0))))
bowtie <- function(x0, y0, d)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0 + d), c(x0 + d, y0), c(x0, y0 + d), c(x0, y0))))
fig8 <- function(x0, y0, d) {
  m <- d / 2
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + m, y0 + m), c(x0 + d, y0), c(x0 + d, y0 + d),
                            c(x0 + m, y0 + m), c(x0, y0 + d), c(x0, y0))))
}
layer <- function(..., crs = 4326) sf::st_sf(geometry = sf::st_sfc(list(...), crs = crs))
box <- function(x0, y0, x1, y1)
  layer(sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)))))
# Positions s2 refuses, each polygon on its own.
s2_refused <- function(x) {
  g <- sf::st_geometry(x)
  spatialkit:::.with_s2(which(vapply(seq_along(g), function(i) inherits(tryCatch(
    quiet(sf::st_as_s2(g[i])), error = function(e) e), "error"), logical(1))))
}
# A ring rebuilt from its coordinates without the vertices that repeat the
# one before them in x and y (independent of the package's helper).
drop_repeats <- function(g1) {
  stopifnot(inherits(g1[[1]], "POLYGON"))
  xy <- sf::st_coordinates(g1)
  parts <- lapply(split(seq_len(nrow(xy)), xy[, "L1"]), function(i) {
    r <- xy[i, c("X", "Y"), drop = FALSE]
    keep <- c(TRUE, r[-1, 1] != r[-nrow(r), 1] | r[-1, 2] != r[-nrow(r), 2])
    unname(r[keep, , drop = FALSE])
  })
  sf::st_sfc(sf::st_polygon(unname(parts)), crs = sf::st_crs(g1))
}
planar_area <- function(x)
  as.numeric(sf::st_area(sf::st_set_crs(sf::st_geometry(x), NA)))
# Is every point of `pts` in `region`, on the sphere?
all_in <- function(pts, region)
  all(lengths(spatialkit:::.with_s2(sf::st_intersects(pts, region))) > 0L)
local_s2 <- function(on, env = parent.frame()) {
  old <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}

# The 0.5 degree square of the scan, and points inside it.
b_sq <- box(-94, 39, -93.5, 39.5)
pts_sq <- local({
  set.seed(1)
  sf::st_sf(v = rnorm(200), geometry = sf::st_sfc(lapply(1:200, function(i)
    sf::st_point(c(-94 + runif(1, 0.05, 0.45), 39 + runif(1, 0.05, 0.45)))), crs = 4326))
})
# The grid as the clip in degrees writes it: create_grid_polygons()'s steps,
# without the drop.
clip_in_degrees <- function(b, cellsize, type) {
  b    <- quiet(sf::st_make_valid(b))
  env  <- sf::st_as_sfc(sf::st_bbox(b))
  grid <- sf::st_make_grid(env, cellsize = c(cellsize, cellsize), what = "polygons",
                           square = identical(type, "square"))
  grid <- sf::st_sf(poly_id = seq_along(grid), geometry = sf::st_sfc(grid, crs = sf::st_crs(b)))
  out <- quiet(sf::st_intersection(sf::st_make_valid(grid), sf::st_union(sf::st_geometry(b))))
  out <- spatialkit:::.keep_polygonal(out)
  out$poly_id <- seq_len(nrow(out))
  out
}

# --- BUG-127 ---------------------------------------------------------------

test_that("a grid clipped in degrees comes back without repeated corners, which s2 accepts", {
  local_s2(TRUE)
  for (ty in c("square", "hex")) {
    raw <- clip_in_degrees(b_sq, 0.1, ty)
    bad <- s2_refused(raw)
    # The scan's grids: one square cell of 25 and three hexagons of 39 on
    # macOS arm64 (s2 1.1.11), two square cells on Linux x86-64 (s2 1.1.13).
    # How many cells s2's clip leaves with a repeated corner follows the last
    # bits of the corners it computes, so the test asks for at least one,
    # which gives the clean-up below something to do, and checks each of
    # them.
    expect_gte(length(bad), 1L)
    g <- quiet(create_grid_polygons(b_sq, cellsize = 0.1, crs = 4326, type = ty, quiet = TRUE))
    expect_identical(nrow(g), nrow(raw))
    expect_length(s2_refused(g), 0L)
    expect_identical(g$poly_id, seq_len(nrow(g)))
    # Cells s2 accepted are returned as they were; the others without the
    # repeat, with their points and planar areas.
    ok <- setdiff(seq_len(nrow(raw)), bad)
    expect_identical(unclass(sf::st_geometry(g))[ok], unclass(sf::st_geometry(raw))[ok])
    for (i in bad) {
      expect_identical(sf::st_geometry(g)[[i]], drop_repeats(sf::st_geometry(raw)[i])[[1]])
      expect_identical(unique(sf::st_coordinates(sf::st_geometry(g)[i])[, 1:2]),
                       unique(sf::st_coordinates(sf::st_geometry(raw)[i])[, 1:2]))
    }
    expect_equal(planar_area(g), planar_area(raw), tolerance = 1e-12)
    # Measured and joined on the sphere.
    a <- spatialkit:::.with_s2(as.numeric(sf::st_area(g)))
    expect_true(all(is.finite(a) & a > 0))
    one <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(-93.75, 39.25)), crs = 4326))
    expect_identical(nrow(quiet(assign_features_to_polygons(one, g, quiet = TRUE))), 1L)
  }
})

test_that("a grid clipped in degrees with corners that differ by rounding is returned whole", {
  local_s2(TRUE)
  # The clip also writes two corners that differ by rounding (7e-15 degrees),
  # which s2 reads as one point or as edges that cross: on the 0.5 degree
  # square, grids at 0.01 to 0.05 degrees stopped with the classed error, and
  # build_tessellation() with them.  The cells as the clip writes them: the
  # repair in degrees turned off.
  raw_grid <- function(cs, ty) {
    testthat::local_mocked_bindings(
      .s2_mend = function(g, split = TRUE) list(g = g, refused = integer(0), mended = integer(0)),
      .package = "spatialkit")
    quiet(create_grid_polygons(b_sq, cellsize = cs, crs = 4326, type = ty, quiet = TRUE))
  }
  near <- 0L
  for (ty in c("square", "hex")) for (cs in c(seq(0.01, 0.13, by = 0.01), 0.025)) {
    raw <- raw_grid(cs, ty)
    bad <- s2_refused(raw)
    g <- quiet(create_grid_polygons(b_sq, cellsize = cs, crs = 4326, type = ty, quiet = TRUE))
    expect_identical(nrow(g), nrow(raw))
    expect_length(s2_refused(g), 0L)
    ok <- setdiff(seq_len(nrow(raw)), bad)
    expect_identical(unclass(sf::st_geometry(g))[ok], unclass(sf::st_geometry(raw))[ok])
    for (i in bad) {
      xy0 <- sf::st_coordinates(sf::st_geometry(raw)[i])[, 1:2]
      xy1 <- sf::st_coordinates(sf::st_geometry(g)[i])[, 1:2]
      # No corner moves: each corner kept is one of the cell's, and each
      # corner dropped lies within 1e-12 degrees of one kept.
      key <- function(xy) paste(sprintf("%.17g", xy[, 1]), sprintf("%.17g", xy[, 2]))
      expect_true(all(key(xy1) %in% key(xy0)))
      gone <- xy0[!key(xy0) %in% key(xy1), , drop = FALSE]
      if (nrow(gone)) {
        near <- near + 1L
        d <- apply(gone, 1, function(p) min(pmax(abs(xy1[, 1] - p[1]), abs(xy1[, 2] - p[2]))))
        expect_true(all(d <= 1e-12))
      }
      a0 <- planar_area(sf::st_geometry(raw)[i])
      expect_lt(abs(planar_area(sf::st_geometry(g)[i]) - a0) / a0, 1e-10)
    }
  }
  # Some cells had a corner dropped that was not an exact repeat.
  expect_gt(near, 0L)
  # The grids that stopped run through build_tessellation().
  for (cfg in list(list(0.01, "square"), list(0.01, "hex"), list(0.02, "hex"),
                   list(0.03, "square"), list(0.04, "hex"), list(0.05, "hex"))) {
    tess <- quiet(build_tessellation(pts_sq, boundary = b_sq, method = cfg[[2]],
                                     cellsize = cfg[[1]], crs = 4326, quiet = TRUE))
    expect_length(s2_refused(tess$cells), 0L)
    expect_false(anyNA(tess$index))
  }
})

test_that("corners within the tolerance are read as one, and no corner moves", {
  drv <- spatialkit:::.drop_repeated_vertices
  ring <- function(m) sf::st_sfc(sf::st_polygon(list(m)), crs = 4326)
  xy <- function(g) unname(sf::st_coordinates(g)[, 1:2])
  e <- 1e-13
  # A corner just after the first, one in the middle, and one just before
  # the closing corner: each dropped, the closing corner kept.
  r <- rbind(c(0, 0), c(e, 0), c(1, 0), c(1, 1), c(1, 1 + e), c(0, 1), c(e, e), c(0, 0))
  expect_identical(xy(drv(ring(r), tol = 1e-12)), rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))
  # Exact repeats only by default.
  expect_identical(xy(drv(ring(r))), r)
  # A ring that would keep fewer than three corners is kept as it was.
  thin <- rbind(c(0, 0), c(1, 0), c(1, e), c(e, 0), c(0, 0))
  expect_identical(xy(drv(ring(thin), tol = 1e-12)), thin)
  # Farther apart than the tolerance: kept.
  wide <- rbind(c(0, 0), c(1e-11, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
  expect_identical(xy(drv(ring(wide), tol = 1e-12)), wide)
})

test_that("build_tessellation() indexes the points in a grid laid in degrees", {
  local_s2(TRUE)
  for (ty in c("square", "hex")) {
    tess <- quiet(build_tessellation(pts_sq, boundary = b_sq, method = ty, cellsize = 0.1,
                                     crs = 4326, quiet = TRUE))
    g <- quiet(create_grid_polygons(b_sq, cellsize = 0.1, crs = 4326, type = ty, quiet = TRUE))
    expect_identical(unclass(sf::st_geometry(tess$cells)), unclass(sf::st_geometry(g)))
    expect_length(tess$index, nrow(pts_sq))
    expect_false(anyNA(tess$index))
    # Each point's cell holds it.
    hit <- spatialkit:::.with_s2(sf::st_intersects(pts_sq, tess$cells))
    expect_true(all(mapply(function(h, i) i %in% h, hit, tess$index)))
  }
})

test_that("grids s2 accepts, and grids not laid in degrees, are those they were", {
  local_s2(TRUE)
  # A clipped grid in degrees s2 accepts whole, and the unclipped lattice:
  # nothing is mended.  Projected and target_cells grids do not reach the
  # repair at all.
  tripped <- function(...) stop("the repair in degrees was reached")
  plain <- function(expr) {
    testthat::local_mocked_bindings(.s2_mend = tripped, .boundary_mend = tripped,
                                    .package = "spatialkit")
    expr
  }
  tri <- layer(sf::st_polygon(list(rbind(c(-94, 39), c(-93.5, 39.1), c(-93.8, 39.5), c(-94, 39)))))
  for (ty in c("square", "hex")) {
    raw <- clip_in_degrees(tri, 0.1, ty)
    expect_length(s2_refused(raw), 0L)
    g1 <- quiet(create_grid_polygons(tri, cellsize = 0.1, crs = 4326, type = ty, quiet = TRUE))
    expect_identical(unclass(sf::st_geometry(g1)), unclass(sf::st_geometry(raw)))
    expect_identical(g1$poly_id, raw$poly_id)
  }
  g2 <- quiet(create_grid_polygons(b_sq, cellsize = 0.1, crs = 4326, type = "hex",
                                   clip = FALSE, quiet = TRUE))
  expect_length(s2_refused(g2), 0L)
  expect_identical(nrow(g2), length(sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(b_sq)),
                                                     cellsize = c(0.1, 0.1), square = FALSE)))
  bp <- sf::st_transform(b_sq, 32615)
  for (ty in c("square", "hex")) {
    expect_identical(
      plain(quiet(create_grid_polygons(bp, cellsize = 10000, type = ty, quiet = TRUE))),
      quiet(create_grid_polygons(bp, cellsize = 10000, type = ty, quiet = TRUE)))
    expect_no_error(plain(quiet(create_grid_polygons(b_sq, target_cells = 25, crs = 4326,
                                                     type = ty, quiet = TRUE))))
    expect_no_error(plain(quiet(build_tessellation(pts_sq, boundary = b_sq, method = ty,
                                                   approx_n_cells = 25, crs = 4326,
                                                   quiet = TRUE))))
    expect_no_error(plain(quiet(build_tessellation(sf::st_transform(pts_sq, 32615),
                                                   boundary = bp, method = ty,
                                                   cellsize = 10000, quiet = TRUE))))
  }
})

test_that("a cell s2 refuses even repaired stops the grid with the classed error", {
  local_s2(TRUE)
  # Rows of a lattice laid in degrees past latitude 90: s2 refuses them.
  pole <- box(-94, 89.45, -93.5, 89.99)
  e <- tryCatch(quiet(create_grid_polygons(pole, cellsize = 0.1, crs = 4326, clip = FALSE,
                                           quiet = TRUE)), error = function(e) e)
  expect_s3_class(e, "spatialkit_cells_degenerate")
  expect_identical(e$site, "grid")
  expect_identical(e$rows, 26:30)
  expect_identical(e$n, 30L)
  expect_true(all(grepl("crosses", e$reasons)))
  expect_null(conditionCall(e))
  msg <- conditionMessage(e)
  expect_match(msg, "^create_grid_polygons\\(\\): s2")
  # Laid with clip = FALSE, the cells were not clipped, and the error does
  # not say they were.
  expect_match(msg, "refuses 5 of the 30 cells of the grid laid in degrees (rows 26, 27, 28, 29 and 30;",
               fixed = TRUE)
  expect_no_match(msg, "clipped", fixed = TRUE)
  expect_match(msg, "Pass a projected `crs` and a `cellsize` in its units, or size the grid by `target_cells` or `n`",
               fixed = TRUE)
  expect_no_match(msg, "four seeds", fixed = TRUE)
  # With the clip, s2's clip meets those rows first.
  e2 <- tryCatch(quiet(create_grid_polygons(pole, cellsize = 0.1, crs = 4326, quiet = TRUE)),
                 error = function(e) e)
  expect_s3_class(e2, "spatialkit_cells_degenerate")
  expect_identical(e2$site, "lattice")
  expect_identical(e2$rows, 26:30)
  expect_match(conditionMessage(e2), "before its clip", fixed = TRUE)
  expect_match(conditionMessage(e2), "so the clip of the grid to `boundary` cannot run in degrees",
               fixed = TRUE)
  # build_tessellation() names its own argument.
  set.seed(2)
  pp <- sf::st_sf(geometry = sf::st_sfc(lapply(1:30, function(i)
    sf::st_point(c(-94 + runif(1, 0.05, 0.45), 89.5 + runif(1, 0, 0.4)))), crs = 4326))
  e3 <- tryCatch(quiet(build_tessellation(pp, boundary = pole, method = "square", cellsize = 0.1,
                                          crs = 4326, quiet = TRUE)), error = function(e) e)
  expect_s3_class(e3, "spatialkit_cells_degenerate")
  expect_identical(e3$site, "lattice")
  expect_match(conditionMessage(e3), "^build_tessellation\\(\\): ")
  expect_match(conditionMessage(e3), "size the grid by `approx_n_cells`", fixed = TRUE)
  expect_match(conditionMessage(e3), "call build_tessellation() again", fixed = TRUE)
  expect_no_match(conditionMessage(e3), "target_cells", fixed = TRUE)
  # The remedies run.
  g <- quiet(create_grid_polygons(pole, target_cells = 20, crs = 4326, quiet = TRUE))
  expect_gt(nrow(g), 0L)
  expect_length(s2_refused(g), 0L)
  expect_gt(nrow(quiet(create_grid_polygons(pole, cellsize = 1000, crs = 3413, quiet = TRUE))), 0L)
})

test_that("a cell with two corners s2 reads as one point names the clip in degrees as the cause", {
  local_s2(TRUE)
  # The box is 20 degrees tall, 4 rows of cells of 5 degrees, but the
  # boundary's height after its repair in degrees is 20.000000000000004, so
  # sf counted a hair over 4 rows and laid a fifth, whose clip is a sliver
  # with corners that differ by rounding, which s2 refuses.  A side that
  # exceeds a whole number of cells by no more than the margin of an edge
  # now gets that number (BUG-112): the lattice is 2 columns by 4 rows, and
  # the clip keeps one part of each of the 4 cells of the west column (s2
  # reads the box's 340 degree edges the short way, across 180 degrees).
  g <- quiet(create_grid_polygons(box(-170, -10, 170, 10), cellsize = c(190, 5),
                                  crs = 4326, quiet = TRUE))
  expect_identical(nrow(g), 4L)
  expect_length(s2_refused(g), 0L)
  # The side counted as sf counts it, as before, still lays that row, and
  # the error names the clip.
  testthat::local_mocked_bindings(.grid_side_counts = function(span, cellsize, bb)
    as.integer(pmax(1, ceiling(span / cellsize))))
  e <- tryCatch(quiet(create_grid_polygons(box(-170, -10, 170, 10), cellsize = c(190, 5),
                                           crs = 4326, quiet = TRUE)), error = function(e) e)
  expect_s3_class(e, "spatialkit_cells_degenerate")
  expect_identical(e$site, "grid")
  expect_match(e$reasons, "duplicate vertex", fixed = TRUE)
  msg <- conditionMessage(e)
  expect_match(msg, "two corners that s2 reads as one point", fixed = TRUE)
  expect_match(msg, "A grid laid in degrees and clipped to `boundary` can hold such corners",
               fixed = TRUE)
  expect_no_match(msg, "four seeds", fixed = TRUE)
  expect_no_match(msg, "returned in longitude and latitude (by", fixed = TRUE)
})

test_that("build_tessellation() gives the classed error when its index meets a cell s2 refuses", {
  local_s2(TRUE)
  raw <- clip_in_degrees(b_sq, 0.1, "square")
  bad <- s2_refused(raw)
  testthat::local_mocked_bindings(create_grid_polygons = function(...) raw)
  e <- tryCatch(quiet(build_tessellation(pts_sq, boundary = b_sq, method = "square",
                                         cellsize = 0.1, crs = 4326, quiet = TRUE)),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_cells_degenerate")
  expect_identical(e$site, "grid")
  expect_identical(e$rows, bad)
  expect_match(conditionMessage(e), "^build_tessellation\\(\\): ")
  expect_match(conditionMessage(e), "the index of the points in the cells cannot run in degrees",
               fixed = TRUE)
  expect_match(conditionMessage(e), "`approx_n_cells`", fixed = TRUE)
})

# --- BUG-128 ---------------------------------------------------------------

test_that("voronoi_seeds_random() repairs a lon/lat boundary before the union, with s2 on and off", {
  for (on in c(TRUE, FALSE)) {
    local_s2(on)
    for (shape in list(sq_rep(-94.2, 39.5, 0.4), fig8(-94.2, 39.5, 0.4))) {
      b <- layer(shape)
      s <- quiet(voronoi_seeds_random(b, 20, set_seed = 1))
      expect_identical(nrow(s), 20L)
      expect_true(all_in(s, spatialkit:::.with_s2(spatialkit:::.s2_mend(sf::st_geometry(b))$g)))
      expect_identical(sf::st_crs(s), sf::st_crs(b))
      expect_identical(nrow(quiet(get_voronoi_seeds(b, n = 20, method = "random", set_seed = 1))), 20L)
    }
    expect_identical(sf::sf_use_s2(), on)
  }
})

test_that("voronoi_seeds_random() gives the rows s2 refuses, and only those, s2's rebuild first", {
  # As get_voronoi_seeds() rebuilds every row with s2 on; valid rows are
  # left as they arrived.
  for (on in c(TRUE, FALSE)) {
    local_s2(on)
    seen <- list()
    testthat::local_mocked_bindings(.safe_make_valid = function(x) {
      seen[[length(seen) + 1L]] <<- list(n = length(sf::st_geometry(x)), s2 = sf::sf_use_s2())
      sf::st_make_valid(x)
    }, .package = "spatialkit")
    b <- layer(sq(-94.8, 39.5, 0.4), sq_rep(-94.2, 39.5, 0.4))
    s <- quiet(voronoi_seeds_random(b, 10, set_seed = 1))
    expect_identical(nrow(s), 10L)
    # Every rebuild is of the one row s2 refuses, with s2 on.  The test of
    # `boundary` as given asks what the rebuild leaves of that row (a ring
    # with no height leaves no polygon) before the rebuild that is kept, so
    # the row is rebuilt twice; the valid row never is.
    expect_gte(length(seen), 1L)
    for (x in seen) expect_identical(x, list(n = 1L, s2 = TRUE))
    seen <- list()
    quiet(voronoi_seeds_random(layer(sq(-94.8, 39.5, 0.4)), 10, set_seed = 1))
    expect_length(seen, 0L)
  }
})

test_that("a projected boundary with an invalid row gives seeds, a cached grid and block folds", {
  local_s2(TRUE)
  s  <- sq(400000, 4300000, 10000)
  bw <- bowtie(400000, 4300000, 10000)
  B2 <- layer(s, bw, crs = 32615)
  seeds <- quiet(voronoi_seeds_random(B2, 10, set_seed = 1))
  expect_identical(nrow(seeds), 10L)
  expect_true(all(lengths(sf::st_intersects(seeds, sf::st_union(sf::st_make_valid(B2)))) > 0L))
  # The cached grid is create_grid_polygons()'s with stable IDs.
  cached <- quiet(create_grid_polygons_cached(B2, target_cells = 16, quiet = TRUE))
  plain  <- quiet(ensure_stable_poly_id(create_grid_polygons(B2, target_cells = 16, quiet = TRUE)))
  expect_identical(nrow(cached), 16L)
  expect_identical(cached, plain)
  # The second call is the cache's.
  expect_identical(quiet(create_grid_polygons_cached(B2, target_cells = 16, quiet = TRUE)), cached)
  set.seed(1)
  pts <- sf::st_sf(v = rnorm(200), geometry = sf::st_sfc(lapply(1:200, function(i)
    sf::st_point(c(400000 + runif(1, 100, 9900), 4300000 + runif(1, 100, 9900)))), crs = 32615))
  f <- quiet(make_folds(pts, k = 4, method = "block_kfold", boundary = B2, block_size = 2000,
                        seed = 1))
  expect_identical(nrow(f$assignment), 200L)
  expect_identical(sort(unlist(lapply(f$folds, `[[`, "test"))), 1:200)
  # The bowtie lies inside the square, so the region is the square's: the
  # same folds as the square alone.
  f1 <- quiet(make_folds(pts, k = 4, method = "block_kfold", boundary = layer(s, crs = 32615),
                         block_size = 2000, seed = 1))
  expect_identical(f$folds, f1$folds)
  expect_identical(f$assignment, f1$assignment)
})

test_that("a boundary of valid rows gives its seeds, cache key and folds as before", {
  # The union as these sites took it before the repair: the rows as they came.
  old_union <- function(b, ...) {
    if (spatialkit:::.is_longlat(b)) spatialkit:::.with_s2(sf::st_union(b)) else sf::st_union(b)
  }
  old <- function(expr) {
    testthat::local_mocked_bindings(.boundary_union = old_union, .package = "spatialkit")
    expr
  }
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  layers <- list(ll = layer(sq(-94.2, 39.5, 0.4), sq(-93.7, 39.5, 0.3)),
                 nc = nc[1:6, ],
                 pr = layer(sq(400000, 4300000, 10000), sq(410000, 4300000, 10000), crs = 32615))
  for (on in c(TRUE, FALSE)) {
    local_s2(on)
    for (b in layers) {
      expect_identical(quiet(voronoi_seeds_random(b, 15, set_seed = 3)),
                       old(quiet(voronoi_seeds_random(b, 15, set_seed = 3))))
      bp <- spatialkit:::.project_for_grid(b, "x")
      expect_identical(spatialkit:::.cache_key(bp, "square", 16),
                       old(spatialkit:::.cache_key(bp, "square", 16)))
    }
  }
  local_s2(TRUE)
  bp <- sf::st_transform(layers$ll, 32615)
  set.seed(4)
  pts <- sf::st_sf(v = rnorm(150), geometry = sf::st_sample(bp, 150))
  expect_identical(quiet(make_folds(pts, k = 4, method = "block_kfold", boundary = bp,
                                    block_size = 5000, seed = 1)),
                   old(quiet(make_folds(pts, k = 4, method = "block_kfold", boundary = bp,
                                        block_size = 5000, seed = 1))))
})

test_that("the cache key hashes the rows when the union cannot be taken", {
  local_s2(TRUE)
  # A square and a bowtie row: GEOS cannot union them as they come.
  B2 <- layer(sq(400000, 4300000, 10000), bowtie(400000, 4300000, 10000), crs = 32615)
  expect_s3_class(try(sf::st_union(sf::st_geometry(B2)), silent = TRUE), "try-error")
  testthat::local_mocked_bindings(.boundary_union = function(...) stop("no union"))
  k <- spatialkit:::.cache_key(B2, "square", 16)
  expect_type(k, "character")
  expect_length(k, 1L)
  expect_identical(spatialkit:::.cache_key(B2, "square", 16), k)
  expect_false(identical(k, spatialkit:::.cache_key(layer(sq(400000, 4300000, 9000), crs = 32615),
                                                    "square", 16)))
})

# --- BUG-129 ---------------------------------------------------------------

test_that("a lon/lat bowtie is split into its lobes for the seeds and the grid in degrees", {
  bt <- layer(bowtie(-94.2, 39.5, 0.4))
  lobes <- spatialkit:::.with_s2(spatialkit:::.split_crossing_edges(sf::st_geometry(bt)))
  expect_false(is.null(lobes))
  for (on in c(TRUE, FALSE)) {
    local_s2(on)
    for (m in c("random", "kmeans")) {
      s <- quiet(get_voronoi_seeds(bt, n = 20, method = m, set_seed = 1))
      expect_identical(nrow(s), 20L)
      expect_true(all_in(s, lobes))
    }
    s <- quiet(voronoi_seeds_random(bt, 20, set_seed = 1))
    expect_identical(nrow(s), 20L)
    expect_true(all_in(s, lobes))
  }
  local_s2(TRUE)
  g <- quiet(create_grid_polygons(bt, cellsize = 0.05, crs = 4326, quiet = TRUE))
  expect_gt(nrow(g), 0L)
  expect_length(s2_refused(g), 0L)
  a_cells <- spatialkit:::.with_s2(sum(as.numeric(sf::st_area(g))))
  a_lobes <- spatialkit:::.with_s2(as.numeric(sf::st_area(lobes)))
  expect_equal(a_cells, a_lobes, tolerance = 1e-6)
  # A bowtie row beside a square: only that row is rebuilt.
  two <- layer(sq(-94.8, 39.5, 0.4), bowtie(-94.2, 39.5, 0.4))
  s <- quiet(voronoi_seeds_random(two, 30, set_seed = 1))
  expect_identical(nrow(s), 30L)
  expect_true(all_in(s, c(sf::st_geometry(two)[1], lobes)))
})

test_that("a valid boundary gives its seeds and its grid in degrees as before", {
  # get_voronoi_seeds() repaired every row before the union; that is kept.
  # In degrees the repair runs on the sphere whatever sf_use_s2() says: with
  # s2 off it ran in the plane, which turned a ring s2 accepts (three
  # vertices on one parallel) into a line, and the seeds of a valid boundary
  # differed from those with s2 on in the last bits.
  old_union <- function(b, ...) {
    if (spatialkit:::.is_longlat(b))
      spatialkit:::.with_s2(sf::st_union(spatialkit:::.safe_make_valid(b)))
    else sf::st_union(spatialkit:::.safe_make_valid(b))
  }
  old <- function(expr) {
    testthat::local_mocked_bindings(.boundary_union = old_union, .package = "spatialkit")
    expr
  }
  b <- layer(sq(-94.2, 39.5, 0.4), sq(-93.7, 39.5, 0.3))
  seeds_with <- list()
  for (on in c(TRUE, FALSE)) {
    local_s2(on)
    for (m in c("random", "kmeans")) {
      s <- quiet(get_voronoi_seeds(b, n = 12, method = m, set_seed = 3))
      expect_identical(s, old(quiet(get_voronoi_seeds(b, n = 12, method = m, set_seed = 3))))
      seeds_with[[paste(m, on)]] <- s
    }
  }
  # The same seeds with s2 on and off.
  for (m in c("random", "kmeans"))
    expect_identical(seeds_with[[paste(m, FALSE)]], seeds_with[[paste(m, TRUE)]], label = m)
  local_s2(TRUE)
  # A valid boundary is clipped as it was: the repair of its rows changes
  # nothing.
  for (bb in list(b_sq, b)) {
    expect_identical(
      quiet(create_grid_polygons(bb, cellsize = 0.05, crs = 4326, quiet = TRUE)),
      local({
        testthat::local_mocked_bindings(.boundary_mend = function(b, ...) b, .package = "spatialkit")
        quiet(create_grid_polygons(bb, cellsize = 0.05, crs = 4326, quiet = TRUE))
      }))
  }
})

test_that("a boundary row s2 refuses even split stops the call with the classed error", {
  local_s2(TRUE)
  bt  <- layer(bowtie(-94.2, 39.5, 0.4))
  two <- layer(sq(-94.8, 39.5, 0.4), bowtie(-94.2, 39.5, 0.4))
  # The split rebuild that does nothing: s2 refuses the bowtie even then.
  testthat::local_mocked_bindings(.split_crossing_edges = function(g1) NULL)
  check <- function(e, caller, rows = 1L, n = 1L) {
    expect_s3_class(e, "spatialkit_cells_degenerate")
    expect_identical(e$site, "boundary")
    expect_identical(e$rows, rows)
    expect_identical(e$n, n)
    expect_match(e$reasons, "crosses", fixed = TRUE)
    expect_s3_class(e$parent, "error")
    expect_null(conditionCall(e))
    msg <- conditionMessage(e)
    expect_match(msg, paste0("^", caller, "\\(\\): s2"))
    expect_match(msg, sprintf("of `boundary` (row %d", rows), fixed = TRUE)
    expect_match(msg, "Edge 0 crosses edge 2", fixed = TRUE)
    expect_no_match(msg, "two corners", fixed = TRUE)
    expect_no_match(msg, "four seeds", fixed = TRUE)
    msg
  }
  for (f in list(
    list("get_voronoi_seeds", function(b) get_voronoi_seeds(b, n = 10, method = "random", set_seed = 1)),
    list("get_voronoi_seeds", function(b) get_voronoi_seeds(b, n = 10, method = "kmeans", set_seed = 1)),
    list("voronoi_seeds_random", function(b) voronoi_seeds_random(b, 10, set_seed = 1)))) {
    msg <- check(tryCatch(quiet(f[[2]](bt)), error = function(e) e), f[[1]])
    expect_match(msg, "Transform `boundary` to a projected CRS", fixed = TRUE)
    expect_match(msg, sprintf("call %s() again", f[[1]]), fixed = TRUE)
    check(tryCatch(quiet(f[[2]](two)), error = function(e) e), f[[1]], rows = 2L, n = 2L)
    # The remedy runs.
    expect_identical(nrow(quiet(f[[2]](sf::st_transform(bt, 32615)))), 10L)
  }
  msg <- check(tryCatch(quiet(create_grid_polygons(bt, cellsize = 0.05, crs = 4326, quiet = TRUE)),
                        error = function(e) e), "create_grid_polygons")
  expect_match(msg, "size the grid by `target_cells` or `n`", fixed = TRUE)
  # With clip = FALSE only the bounding box of `boundary` is read: the grid
  # is laid as it was, whatever s2 says of the rows.
  g <- quiet(create_grid_polygons(bt, cellsize = 0.05, crs = 4326, clip = FALSE, quiet = TRUE))
  expect_identical(g, quiet(create_grid_polygons(box(-94.2, 39.5, -93.8, 39.9), cellsize = 0.05,
                                                 crs = 4326, clip = FALSE, quiet = TRUE)))
  expect_length(s2_refused(g), 0L)
  expect_gt(nrow(quiet(create_grid_polygons(bt, target_cells = 40, crs = 4326, quiet = TRUE))), 0L)
  set.seed(3)
  pp <- sf::st_sf(geometry = sf::st_sfc(lapply(1:20, function(i)
    sf::st_point(c(-94.15 + runif(1, 0, 0.05), 39.6 + runif(1, 0, 0.2)))), crs = 4326))
  msg <- check(tryCatch(quiet(build_tessellation(pp, boundary = bt, method = "hex", cellsize = 0.05,
                                                 crs = 4326, quiet = TRUE)),
                        error = function(e) e), "build_tessellation")
  expect_match(msg, "size the grid by `approx_n_cells`", fixed = TRUE)
  # With s2 off the error says that the package used s2 here.
  local_s2(FALSE)
  e <- tryCatch(quiet(voronoi_seeds_random(bt, 10, set_seed = 1)), error = function(e) e)
  expect_true(isTRUE(e$s2_forced))
})
