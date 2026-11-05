# tests/testthat/test-fix-join-cells-degenerate.R
# ---------------------------------------------------------------------------
# A join in longitude and latitude runs on the sphere, where s2 refuses a
# polygon with two corners it reads as one point.  Cells returned in lon/lat
# hold such corners on gridded locations, and the joins that do not run in
# the record's CRS (cells without the record, cells whose record says
# scored_cells = FALSE, cells without a CRS that keep a record) stopped with
# s2's bare "Loop 0 is not valid: Edge 20 is degenerate (duplicate vertex)".
# They stop with an error of class "spatialkit_cells_degenerate" that names
# the polygons s2 refuses, the cause and the transform to a projected CRS;
# so do the other refusals of a polygon by s2 (a hole that crosses its outer
# ring, two holes that share an edge), on the largest-overlap join too when
# the repair leaves the polygon refused, and cell_values_at() gives the
# error in its own words.  Every other error keeps its own message.
# ---------------------------------------------------------------------------

.jd_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# s2 on for the length of the calling test, as it is by default.
.jd_s2 <- function(on = TRUE, env = parent.frame()) {
  old <- suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}

# 30 squares of 0.1 degree in EPSG:4326, 6 across and 5 up from (-94.3,
# 39.6), numbered by row from the lower left.  `repeat_in` lists the squares
# whose second corner is listed twice: GEOS accepts such a ring, and s2
# refuses it ("duplicate vertex").  `bowtie_in` lists the squares drawn as a
# ring whose edges cross.
.jd_square <- function(x0, y0, d = 0.1, repeated = FALSE, bowtie = FALSE) {
  m <- if (bowtie) rbind(c(x0, y0), c(x0 + d, y0 + d), c(x0 + d, y0), c(x0, y0 + d), c(x0, y0))
       else rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d), c(x0, y0))
  if (repeated) m <- rbind(m[1:2, ], m[2:5, ])
  sf::st_polygon(list(m))
}
.jd_squares <- function(repeat_in = integer(0), bowtie_in = integer(0), d = 0.1) {
  k <- seq_len(30L)
  x0 <- -94.3 + d * ((k - 1L) %% 6L)
  y0 <- 39.6 + d * ((k - 1L) %/% 6L)
  g <- lapply(k, function(i) .jd_square(x0[i], y0[i], d, i %in% repeat_in, i %in% bowtie_in))
  sf::st_sf(cell_id = k, geometry = sf::st_sfc(g, crs = 4326))
}
# A point a little off the middle of each square, numbered as the squares.
.jd_points <- function(d = 0.1) {
  k <- seq_len(30L)
  sf::st_as_sf(data.frame(x = -94.3 + d * ((k - 1L) %% 6L) + 0.4 * d,
                          y = 39.6 + d * ((k - 1L) %/% 6L) + 0.45 * d, row = k),
               coords = c("x", "y"), crs = 4326)
}
# Small squares around those points, the second corner of the first listed
# twice when `repeated`.
.jd_patches <- function(repeated = FALSE) {
  xy <- sf::st_coordinates(.jd_points())
  g <- lapply(seq_len(nrow(xy)), function(i)
    .jd_square(xy[i, 1L] - 0.01, xy[i, 2L] - 0.01, 0.02, repeated = repeated && i == 1L))
  sf::st_sf(row = seq_len(nrow(xy)), geometry = sf::st_sfc(g, crs = 4326))
}

# The squares with square `i` given holes s2 refuses: "crossing", a hole
# that runs out across the square's right edge ("Loop 0 edge 1 crosses loop
# 1 edge 0"); "twin", two identical holes ("Loop 1 edge 0 has duplicate near
# loop 2 edge 1"); "same_way", with `oriented`, a hole that runs the way its
# outer ring runs ("Inconsistent loop orientations detected").  The holes
# keep clear of the square's point.
.jd_holed <- function(i, kind = c("crossing", "twin", "same_way"), oriented = FALSE) {
  kind <- match.arg(kind)
  cells <- .jd_squares()
  b <- sf::st_bbox(cells[i, ])
  x0 <- b[["xmin"]]; y0 <- b[["ymin"]]
  ring <- function(x, y, w, h) rbind(c(x, y), c(x + w, y), c(x + w, y + h), c(x, y + h), c(x, y))
  shell <- ring(x0, y0, 0.1, 0.1)
  g <- switch(kind,
    crossing = list(shell, ring(x0 + 0.05, y0 + 0.02, 0.1, 0.06)),
    twin     = { h <- ring(x0 + 0.06, y0 + 0.06, 0.02, 0.02); list(shell, h, h) },
    same_way = list(shell, ring(x0 + 0.06, y0 + 0.06, 0.02, 0.02)))
  geom <- sf::st_geometry(cells)
  geom[[i]] <- sf::st_polygon(g)
  if (oriented) attr(geom, "oriented") <- TRUE
  sf::st_set_geometry(cells, geom)
}

.jd_assign <- function(x, cells, ...)
  assign_features_to_polygons(x, cells, polygon_id_col = "cell_id", quiet = TRUE, ...)
.jd_catch <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)

# A lattice of `nx` x `ny` locations at spacing `sp` in EPSG:32615 with a
# smooth field plus noise drawn under seed 9001, as in the gridded-ties tests.
.jd_lattice <- function(nx, ny, sp) {
  g <- expand.grid(x = 4e5 + (seq_len(nx) - 1) * sp, y = 4.4e6 + (seq_len(ny) - 1) * sp)
  z <- local({
    .jd_seed(9001)
    sin(g$x / (sp * 6)) + cos(g$y / (sp * 5)) + stats::rnorm(nrow(g), 0, 0.4)
  })
  sf::st_as_sf(data.frame(x = g$x, y = g$y, z = z, row = seq_len(nrow(g))),
               coords = c("x", "y"), crs = 32615)
}

test_that("rows in degrees on lon/lat polygons s2 refuses stop with an error that names them, the cause and the transform", {
  .jd_s2(TRUE)
  pts <- .jd_points()
  cells <- .jd_squares(repeat_in = c(7L, 20L))
  expect_identical(which(!sf::st_is_valid(cells)), c(7L, 20L))
  expect_true(all(sf::st_is_valid(sf::st_transform(cells, 32615))))
  err <- .jd_catch(.jd_assign(pts, cells))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, c(7L, 20L))
  expect_identical(err$n, 30L)
  expect_length(err$reasons, 2L)
  expect_true(all(grepl("duplicate vertex", err$reasons, fixed = TRUE)))
  msg <- conditionMessage(err)
  expect_match(msg, paste0("assign_features_to_polygons(): s2, the engine sf uses for geometry ",
                           "in longitude and latitude, refuses 2 of the 30 polygons of ",
                           "`polygons_sf` (rows 7 and 20; \"Loop 0 is not valid: "), fixed = TRUE)
  expect_match(msg, "duplicate vertex)\"), so the join cannot run in degrees.", fixed = TRUE)
  expect_match(msg, "Each has two corners that s2 reads as one point.", fixed = TRUE)
  expect_match(msg, paste0("two cells meet there at two corners that differ only by rounding, ",
                           "which the transform to degrees puts at one point."), fixed = TRUE)
  expect_match(msg, paste0("Transform `polygons_sf` to a projected CRS, such as the one its cells ",
                           "were built in (sf::st_transform()) or the one ensure_projected() ",
                           "picks for it, and assign again; the join then runs in that CRS."),
               fixed = TRUE)
  expect_null(conditionCall(err))

  # One square, and more than six.
  err1 <- .jd_catch(.jd_assign(pts, .jd_squares(repeat_in = 7L)))
  expect_s3_class(err1, "spatialkit_cells_degenerate")
  expect_identical(err1$rows, 7L)
  expect_match(conditionMessage(err1), "refuses 1 of the 30 polygons of `polygons_sf` (row 7; ",
               fixed = TRUE)
  expect_match(conditionMessage(err1), "It has two corners that s2 reads as one point.",
               fixed = TRUE)
  err8 <- .jd_catch(.jd_assign(pts, .jd_squares(repeat_in = c(2L, 4L, 6L, 8L, 10L, 12L, 14L, 16L))))
  expect_s3_class(err8, "spatialkit_cells_degenerate")
  expect_identical(err8$rows, c(2L, 4L, 6L, 8L, 10L, 12L, 14L, 16L))
  expect_match(conditionMessage(err8), "(rows 2, 4, 6, 8, 10, 12 and 2 more; ", fixed = TRUE)
  # Six are all named, and seven name six and one more.
  err6 <- .jd_catch(.jd_assign(pts, .jd_squares(repeat_in = c(1L, 3L, 5L, 7L, 9L, 11L))))
  expect_s3_class(err6, "spatialkit_cells_degenerate")
  expect_match(conditionMessage(err6),
               "refuses 6 of the 30 polygons of `polygons_sf` (rows 1, 3, 5, 7, 9 and 11; ",
               fixed = TRUE)
  err7 <- .jd_catch(.jd_assign(pts, .jd_squares(repeat_in = c(1L, 3L, 5L, 7L, 9L, 11L, 13L))))
  expect_match(conditionMessage(err7), "(rows 1, 3, 5, 7, 9, 11 and 1 more; ", fixed = TRUE)

  # Polygon features joined by `predicate` meet the same refusal.
  errp <- .jd_catch(.jd_assign(.jd_patches(), cells, largest = FALSE))
  expect_s3_class(errp, "spatialkit_cells_degenerate")
  expect_identical(errp$rows, c(7L, 20L))

  # The transform the error names lets the join run, with every point in its
  # square, as on the squares without the repeated corners.
  want <- seq_len(30L)
  got <- .jd_assign(pts, sf::st_transform(cells, 32615))
  expect_identical(got$cell_id[match(want, got$row)], want)
  got <- suppressWarnings(.jd_assign(pts, ensure_projected(cells)))
  expect_identical(got$cell_id[match(want, got$row)], want)
  got <- .jd_assign(pts, .jd_squares())
  expect_identical(got$cell_id[match(want, got$row)], want)
  got <- .jd_assign(pts, cells[-c(7L, 20L), ], keep_unassigned = TRUE)
  expect_identical(got$cell_id[match(want, got$row)], replace(want, c(7L, 20L), NA_integer_))
  # Rows in a projected CRS are joined there, and GEOS accepts the squares.
  got <- .jd_assign(sf::st_transform(pts, 32615), cells)
  expect_identical(got$cell_id[match(want, got$row)], want)
})

test_that("with s2 off the join in degrees runs on the polygons s2 refuses", {
  .jd_s2(FALSE)
  pts <- .jd_points()
  got <- suppressMessages(.jd_assign(pts, .jd_squares(repeat_in = c(7L, 20L))))
  expect_identical(got$cell_id[match(seq_len(30L), got$row)], seq_len(30L))
})

test_that("a ring s2 refuses for another reason is named by s2's reason, without the corners' cause", {
  .jd_s2(TRUE)
  pts <- .jd_points()
  err <- .jd_catch(.jd_assign(pts, .jd_squares(bowtie_in = 9L)))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, 9L)
  msg <- conditionMessage(err)
  expect_match(msg, "refuses 1 of the 30 polygons of `polygons_sf` (row 9; \"Loop 0 is not valid: ",
               fixed = TRUE)
  expect_match(msg, "crosses edge", fixed = TRUE)
  expect_no_match(msg, "two corners", fixed = TRUE)
  expect_match(msg, "Transform `polygons_sf` to a projected CRS", fixed = TRUE)

  # Both kinds in one layer: the corners' cause for the one, s2's reason for
  # the other.
  err <- .jd_catch(.jd_assign(pts, .jd_squares(repeat_in = 7L, bowtie_in = 9L)))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, c(7L, 9L))
  msg <- conditionMessage(err)
  expect_match(msg, paste0("1 of them has two corners that s2 reads as one point, and s2 refuses ",
                           "the other 1 for another reason (row 9: \""), fixed = TRUE)
  expect_match(msg, "crosses edge", fixed = TRUE)
  expect_match(msg, "two cells meet there at two corners that differ only by rounding",
               fixed = TRUE)
})

test_that("polygons without a CRS that s2 refuses are named, with the CRS to set before the transform", {
  .jd_s2(TRUE)
  pts <- .jd_points()
  cells <- sf::st_set_crs(.jd_squares(repeat_in = c(7L, 20L)), NA)
  err <- .jd_catch(.jd_assign(pts, cells))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, c(7L, 20L))
  expect_match(conditionMessage(err), paste0(
    "`polygons_sf` has no CRS and was read as longitude and latitude: give it its CRS with ",
    "sf::st_set_crs(), transform it to a projected CRS, such as the one its cells were built ",
    "in, with sf::st_transform(), and assign again; the join then runs in that CRS."),
    fixed = TRUE)
  got <- .jd_assign(pts, sf::st_transform(sf::st_set_crs(cells, 4326), 32615))
  expect_identical(got$cell_id[match(seq_len(30L), got$row)], seq_len(30L))
})

test_that("an error of another kind from the join keeps its own message", {
  .jd_s2(TRUE)
  pts <- .jd_points()
  stopper <- function(x, y, ...) stop("the predicate stopped")
  for (cells in list(.jd_squares(), .jd_squares(repeat_in = c(7L, 20L)))) {
    err <- .jd_catch(.jd_assign(pts, cells, predicate = stopper))
    expect_s3_class(err, "simpleError")
    expect_false(inherits(err, "spatialkit_cells_degenerate"))
    expect_identical(conditionMessage(err), "the predicate stopped")
  }
  # The predicate's own condition comes out as it was raised: its class and
  # its call, not a copy.
  own <- errorCondition("the predicate stopped", class = "jd_predicate_error",
                        call = quote(stopper(x, y)))
  stopper <- function(x, y, ...) stop(own)
  for (cells in list(.jd_squares(), .jd_squares(repeat_in = c(7L, 20L)))) {
    err <- .jd_catch(.jd_assign(pts, cells, predicate = stopper))
    expect_identical(err, own)
  }
  # s2 refuses a feature and no polygon: the error does not blame the
  # polygons.
  err <- .jd_catch(.jd_assign(.jd_patches(repeated = TRUE), .jd_squares(), largest = FALSE))
  expect_s3_class(err, "error")
  expect_false(inherits(err, "spatialkit_cells_degenerate"))
  expect_no_match(conditionMessage(err), "polygons_sf", fixed = TRUE)
})

test_that("the largest-overlap join repairs the polygons s2 refuses, and names them if s2 still refuses one", {
  .jd_s2(TRUE)
  cells <- .jd_squares(repeat_in = c(7L, 20L))
  got <- NULL
  expect_warning(got <- .jd_assign(.jd_patches(), cells, largest = TRUE),
                 "0 of 30 feature(s) and 2 of 30 polygon(s) have invalid geometry", fixed = TRUE)
  expect_identical(got$cell_id[match(seq_len(30L), got$row)], seq_len(30L))
  # Without the repair the join meets s2's refusal: `largest = FALSE` would
  # meet it too, so the error names the cause and the transform.
  local_mocked_bindings(.safe_make_valid = function(g) g)
  err <- .jd_catch(.jd_assign(.jd_patches(), cells, largest = TRUE))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, c(7L, 20L))
  expect_no_match(conditionMessage(err), "Pass largest = FALSE", fixed = TRUE)
  expect_match(conditionMessage(err), "Transform `polygons_sf` to a projected CRS", fixed = TRUE)
  # Polygons without a CRS get the remedy that sets it first.
  err <- .jd_catch(.jd_assign(.jd_patches(), sf::st_set_crs(cells, NA), largest = TRUE))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, c(7L, 20L))
  expect_match(conditionMessage(err),
               "`polygons_sf` has no CRS and was read as longitude and latitude: give it its CRS",
               fixed = TRUE)
})

test_that("the other refusals of a polygon by s2 stop both joins with the error, with s2's reason and no corners' cause", {
  .jd_s2(TRUE)
  pts <- .jd_points()
  want <- seq_len(30L)
  check <- function(err, reason) {
    expect_s3_class(err, "spatialkit_cells_degenerate")
    expect_identical(err$rows, 9L)
    expect_identical(err$n, 30L)
    expect_match(err$reasons, reason, fixed = TRUE)
    msg <- conditionMessage(err)
    expect_match(msg, "refuses 1 of the 30 polygons of `polygons_sf` (row 9; \"", fixed = TRUE)
    expect_match(msg, paste0(reason, "[^\"]*\"\\), so the join cannot run in degrees\\."))
    expect_no_match(msg, "two corners", fixed = TRUE)
    expect_no_match(msg, "Pass largest = FALSE", fixed = TRUE)
    expect_match(msg, "Transform `polygons_sf` to a projected CRS", fixed = TRUE)
    expect_null(conditionCall(err))
    expect_s3_class(err$parent, "error")
    expect_match(conditionMessage(err$parent), reason, fixed = TRUE)
  }
  # A hole that crosses its outer ring: the plain join, and the largest-overlap
  # join, whose repair leaves two rings that s2 still refuses.
  cells <- .jd_holed(9L, "crossing")
  expect_identical(which(!sf::st_is_valid(cells)), 9L)
  check(.jd_catch(.jd_assign(pts, cells)), "Loop 0 edge 1 crosses loop 1 edge 0")
  err <- NULL
  expect_warning(err <- tryCatch(.jd_assign(.jd_patches(), cells, largest = TRUE),
                                 error = function(e) e),
                 "0 of 30 feature(s) and 1 of 30 polygon(s) have invalid geometry", fixed = TRUE)
  check(err, "Loop 0 edge 1 crosses loop 1 edge 0")
  # Beside a square with a repeated corner, which the repair makes one s2
  # accepts, the largest-overlap join names the hole's polygon alone.
  both <- cells
  geom <- sf::st_geometry(both)
  b7 <- sf::st_bbox(both[7L, ])
  geom[[7L]] <- .jd_square(b7[["xmin"]], b7[["ymin"]], repeated = TRUE)
  both <- sf::st_set_geometry(both, geom)
  expect_identical(which(!sf::st_is_valid(both)), c(7L, 9L))
  err <- .jd_catch(.jd_assign(.jd_patches(), both, largest = TRUE))
  check(err, "Loop 0 edge 1 crosses loop 1 edge 0")
  err <- .jd_catch(.jd_assign(pts, both))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, c(7L, 9L))
  # The transform the error names lets both joins run.
  got <- .jd_assign(pts, sf::st_transform(cells, 32615))
  expect_identical(got$cell_id[match(want, got$row)], want)
  got <- suppressWarnings(.jd_assign(.jd_patches(), sf::st_transform(cells, 32615), largest = TRUE))
  expect_identical(got$cell_id[match(want, got$row)], want)

  # Two identical holes: the plain join stops; the largest-overlap join
  # repairs the polygon into one s2 accepts, and runs.
  cells <- .jd_holed(9L, "twin")
  check(.jd_catch(.jd_assign(pts, cells)), "has duplicate near loop")
  got <- suppressWarnings(.jd_assign(.jd_patches(), cells, largest = TRUE))
  expect_identical(got$cell_id[match(want, got$row)], want)
  got <- .jd_assign(pts, sf::st_transform(cells, 32615))
  expect_identical(got$cell_id[match(want, got$row)], want)

  # Oriented polygons whose hole runs the way the outer ring runs.
  cells <- .jd_holed(9L, "same_way", oriented = TRUE)
  check(.jd_catch(.jd_assign(pts, cells)), "Inconsistent loop orientations detected")
})

# The 30 x 20 lattice at 30 m, its rows in EPSG:4326, and the k-means cells
# of the size the profile scored, returned in EPSG:4326 (built once).
.jd_kmeans_cache <- new.env(parent = emptyenv())
.jd_kmeans_cells <- function() {
  if (is.null(.jd_kmeans_cache$x)) {
    lat <- .jd_lattice(30L, 20L, 30)
    pts <- sf::st_transform(lat, 4326)
    prof <- suppressWarnings(resolution_profile(pts, "z", partition = "kmeans", nstart = 3,
                                                sac = "auto", seed = 123, estimate = "pooled",
                                                obs_per_cell = 4.5))
    sel <- select_resolution(prof, criterion = "size")
    seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
    tess <- suppressWarnings(build_tessellation(seeds, method = "voronoi", crs = 4326, quiet = TRUE))
    .jd_kmeans_cache$x <- list(lat = lat, pts = pts, cells = tess$cells)
  }
  .jd_kmeans_cache$x
}

test_that("the cells of a lattice returned in lon/lat without their record, or without their CRS, stop with the error; with it they join", {
  .jd_s2(TRUE)
  km <- .jd_kmeans_cells()
  lat <- km$lat
  pts <- km$pts
  cells <- km$cells
  rec <- attr(cells, "count_rule")
  expect_identical(rec$scored_cells, TRUE)
  bad <- which(!sf::st_is_valid(cells))
  skip_if(!length(bad), "s2 refuses none of these cells on this platform")

  # With the record the join runs in the record's CRS.
  scored <- .jd_assign(pts, cells)
  expect_identical(nrow(scored), 600L)
  expect_false(anyNA(scored$cell_id))

  # (1) Without the record.
  norec <- cells; attr(norec, "count_rule") <- NULL
  err <- .jd_catch(.jd_assign(pts, norec))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, bad)
  expect_identical(err$n, nrow(cells))
  expect_match(conditionMessage(err), "Transform `polygons_sf` to a projected CRS", fixed = TRUE)
  expect_match(conditionMessage(err), "Each has two corners that s2 reads as one point.",
               fixed = TRUE)
  got <- .jd_assign(pts, sf::st_transform(norec, rec$crs))
  expect_identical(got$cell_id[match(scored$row, got$row)], scored$cell_id)
  # Rows in the projected CRS are joined there.
  expect_identical(nrow(.jd_assign(lat, norec)), 600L)

  # (3) Without a CRS, with the record.
  bare <- sf::st_set_crs(cells, NA)
  expect_false(is.null(attr(bare, "count_rule")))
  err <- .jd_catch(.jd_assign(pts, bare))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, bad)
  expect_match(conditionMessage(err), "`polygons_sf` has no CRS and was read as longitude and latitude",
               fixed = TRUE)
})

test_that("cell_values_at() without seeds on lon/lat summary cells s2 refuses names `values` and its rows", {
  .jd_s2(TRUE)
  km <- .jd_kmeans_cells()
  pts <- km$pts
  cells <- km$cells
  bad <- which(!sf::st_is_valid(cells))
  skip_if(!length(bad), "s2 refuses none of these cells on this platform")
  asg <- .jd_assign(pts, cells)
  sm <- suppressWarnings(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = cells))
  expect_identical(which(!sf::st_is_valid(sm)), bad)
  rec <- attr(sm, "count_rule")
  nd <- pts[seq(1L, nrow(pts), by = 7L), ]
  at <- function(values) suppressWarnings(cell_values_at(nd, values, asg, "z", id_col = "cell_id"))
  # With the record the locations are placed in the record's CRS.
  ref <- at(sm)
  expect_identical(nrow(ref), nrow(nd))
  expect_false(anyNA(ref$cell_id))

  # Without the record they go to s2, and the error is this function's.
  norec <- sm; attr(norec, "count_rule") <- NULL
  err <- .jd_catch(at(norec))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, bad)
  expect_identical(err$n, nrow(sm))
  expect_length(err$reasons, length(bad))
  expect_s3_class(err$parent, "error")
  expect_null(conditionCall(err))
  msg <- conditionMessage(err)
  expect_match(msg, sprintf(paste0("cell_values_at(): s2, the engine sf uses for geometry in longitude ",
                                   "and latitude, refuses %d of the %d polygons of `values` (rows "),
                            length(bad), nrow(sm)), fixed = TRUE)
  expect_match(msg, "so the placement of the locations in its polygons cannot run in degrees.",
               fixed = TRUE)
  expect_match(msg, paste0("Transform `values` to a projected CRS, such as the one its cells were ",
                           "built in (sf::st_transform()) or the one ensure_projected() picks for ",
                           "it, and call cell_values_at() again; the locations are then placed in ",
                           "that CRS."), fixed = TRUE)
  expect_no_match(msg, "polygons_sf", fixed = TRUE)
  expect_no_match(msg, "assign_features_to_polygons", fixed = TRUE)
  expect_no_match(msg, "assign again", fixed = TRUE)
  # The transform it names places every location as the record did.
  got <- at(sf::st_transform(norec, rec$crs))
  expect_identical(got$cell_id, ref$cell_id)
  expect_identical(got$value, ref$value)

  # A row of `values` without an ID is no cell: the rows named are still
  # those of `values`.
  noid <- norec; noid$cell_id[2L] <- NA
  expect_false(2L %in% bad)
  err <- .jd_catch(at(noid))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, bad)
  expect_match(conditionMessage(err), sprintf(
    "refuses %d of the %d polygons of `values` that carry an ID (rows ", length(bad), nrow(sm) - 1L),
    fixed = TRUE)

  # Without a CRS: the CRS to set before the transform.
  err <- .jd_catch(at(sf::st_set_crs(sm, NA)))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, bad)
  expect_match(conditionMessage(err), paste0(
    "`values` has no CRS and was read as longitude and latitude: give it its CRS with ",
    "sf::st_set_crs(), transform it to a projected CRS, such as the one its cells were built ",
    "in, with sf::st_transform(), and call cell_values_at() again"), fixed = TRUE)
  got <- at(sf::st_transform(sf::st_set_crs(sf::st_set_crs(sm, NA), 4326), rec$crs))
  expect_identical(got$cell_id, ref$cell_id)
})

test_that("the cells of centres refitted to every point, returned in lon/lat, stop with the error where s2 refuses one", {
  .jd_s2(TRUE)
  pts <- sf::st_transform(.jd_lattice(50L, 40L, 250), 4326)
  prof <- suppressWarnings(resolution_profile(pts, "z", partition = "density", nstart = 3,
                                              seed = 123, estimate = NULL))
  sel <- suppressWarnings(select_resolution(prof))
  seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  tess <- suppressWarnings(build_tessellation(seeds, method = "voronoi", crs = 4326, quiet = TRUE))
  cells <- tess$cells
  rec <- attr(cells, "count_rule")
  expect_identical(rec$scored_cells, FALSE)
  bad <- which(!sf::st_is_valid(cells))
  skip_if(!length(bad), "s2 refuses none of these cells on this platform")
  # (2) The record does not say a profile scored them, so they are joined in
  # degrees.
  err <- .jd_catch(.jd_assign(pts, cells))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, bad)
  expect_match(conditionMessage(err), "Transform `polygons_sf` to a projected CRS", fixed = TRUE)
  got <- .jd_assign(pts, sf::st_transform(cells, rec$crs))
  expect_identical(nrow(got), nrow(pts))
  expect_false(anyNA(got$cell_id))
})

test_that("?assign_features_to_polygons, ?create_voronoi_polygons, ?cell_values_at and NEWS name the error", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "assign_features_to_polygons.Rd")
  expect_match(rd, paste0("In longitude and latitude, with s2 on (\\code{sf::sf_use_s2()}, sf's ",
                          "default), sf joins on the sphere with s2, which refuses a polygon that ",
                          "is not valid there, such as one with two corners it reads as one point ",
                          "or with edges that cross;"), fixed = TRUE)
  expect_match(rd, paste0("A join by \\code{predicate} then stops with an error of class ",
                          "\\code{\"spatialkit_cells_degenerate\"} that gives the rows of the ",
                          "polygons s2 refuses (the condition carries them as \\code{rows}), names ",
                          "the cause and says to transform \\code{polygons_sf} to a projected CRS, ",
                          "where the join runs. The largest-overlap join repairs such polygons ",
                          "first (see \\code{largest}), and stops with that error when s2 still ",
                          "refuses a polygon it repaired, such as one whose hole crosses its outer ",
                          "ring."), fixed = TRUE)
  expect_no_match(rd, "The largest-overlap join repairs such polygons instead", fixed = TRUE)
  rd <- read("man", "cell_values_at.Rd")
  expect_match(rd, paste0("Where that join runs on the sphere and s2 refuses a polygon of ",
                          "\\code{values} (see Details in ",
                          "\\code{\\link[=assign_features_to_polygons]{assign_features_to_polygons()}}), ",
                          "the call stops with an error of class ",
                          "\\code{\"spatialkit_cells_degenerate\"} that gives the rows of ",
                          "\\code{values} s2 refuses (\\code{rows} on the condition), names the ",
                          "cause and says to transform \\code{values} to a projected CRS."),
               fixed = TRUE)
  rd <- read("man", "create_voronoi_polygons.Rd")
  expect_match(rd, paste0("then stops with an error of class \\code{\"spatialkit_cells_degenerate\"} ",
                          "that names the cells s2 refuses, the cause and that transform; its ",
                          "largest-overlap join of polygon features repairs such cells instead, ",
                          "with a warning."), fixed = TRUE)
  news <- read("NEWS.md")
  expect_match(news, paste0("The join by `predicate` now stops with an error of class ",
                            "`spatialkit_cells_degenerate` that gives the number and rows of the ",
                            "polygons s2 refuses"), fixed = TRUE)
  expect_match(news, paste0("on a polygon s2 refuses even after the repair (a hole that crosses ",
                            "its outer ring) it stops with the same error, where it advised ",
                            "`largest = FALSE`, which meets the same refusal."), fixed = TRUE)
  expect_match(news, paste0("`cell_values_at()` without `seeds`, which places its locations with ",
                            "this join, gives the error in its own words: it names `values` and the ",
                            "rows of `values`."), fixed = TRUE)
})
