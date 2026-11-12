# tests/testthat/test-fix-join-features-degenerate.R
# ---------------------------------------------------------------------------
# assign_features_to_polygons() on layers in longitude and latitude.
#
# s2 refuses a polygon of `features_sf` as it refuses one of `polygons_sf`:
# the join by `predicate` stopped with s2's bare message, and the
# largest-overlap join, whose repair makes s2 accept a figure eight as two
# parts that touch at a point, stopped on a piece of an overlap and advised
# `largest = FALSE`, which stopped too.  Both joins now stop with an error
# that names the layer s2 refused, its rows and the transform:
# "spatialkit_features_degenerate" for the features, the cells' class for
# the polygons.  Whether s2 refuses that piece follows the last bits of the
# vertices it computes, which differ between platforms; where s2 accepts it
# the largest-overlap join runs (.jf_pieces_accepted()).
#
# The smallest-area tie-break read GEOS's area of a ring whose edges cross
# (its lobes cancel), and gave a point to a bowtie as the smaller polygon.
# With sf::sf_use_s2(FALSE) the tie-break and the largest-overlap join
# stopped with "package lwgeom required"; they now measure on s2 in degrees
# whatever the session says.
# ---------------------------------------------------------------------------

# s2 on (or off) for the length of the calling test.
.jf_s2 <- function(on = TRUE, env = parent.frame()) {
  old <- suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}

.jf_sq <- function(x0, y0, d) rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d),
                                    c(x0, y0 + d), c(x0, y0))
# Six 0.1 degree squares at 39.6N in EPSG:4326, numbered from the west.
.jf_squares <- function()
  sf::st_sf(cell_id = 1:6, geometry = sf::st_sfc(lapply(0:5, function(i)
    sf::st_polygon(list(.jf_sq(-94.3 + 0.1 * i, 39.6, 0.1)))), crs = 4326))
.jf_points <- function()
  sf::st_as_sf(data.frame(x = -94.3 + 0.1 * (0:5) + 0.04, y = 39.645, row = 1:6),
               coords = c("x", "y"), crs = 4326)
# A 0.02 degree patch around each point.
.jf_patches <- function() {
  xy <- sf::st_coordinates(.jf_points())
  sf::st_sf(row = 1:6, geometry = sf::st_sfc(lapply(1:6, function(i)
    sf::st_polygon(list(.jf_sq(xy[i, 1] - 0.01, xy[i, 2] - 0.01, 0.02)))), crs = 4326))
}
.jf_put <- function(layer, i, g) {
  gg <- sf::st_geometry(layer); gg[[i]] <- g; sf::st_set_geometry(layer, gg)
}
# Cell 3 a figure eight: one vertex listed twice, not one after the other.
.jf_cells8 <- function() {
  x0 <- -94.1; y0 <- 39.6
  .jf_put(.jf_squares(), 3L, sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + 0.05, y0 + 0.05), c(x0 + 0.1, y0), c(x0 + 0.1, y0 + 0.1),
    c(x0 + 0.05, y0 + 0.05), c(x0, y0 + 0.1), c(x0, y0)))))
}
# Patch 2 with a vertex listed twice, patch 4 a bowtie, patch 3 a figure
# eight inside cell 3, patch 3 a figure eight pinched on the edge of cells 2
# and 3, patch 5 with a hole that crosses its outer ring.
.jf_features <- function(kind = c("repeat", "bowtie", "fig8", "edge8", "hole")) {
  kind <- match.arg(kind)
  p <- .jf_patches()
  xy <- sf::st_coordinates(.jf_points())
  switch(kind,
    "repeat" = {
      m <- .jf_sq(xy[2, 1] - 0.01, xy[2, 2] - 0.01, 0.02)
      .jf_put(p, 2L, sf::st_polygon(list(rbind(m[1:2, ], m[2:5, ]))))
    },
    bowtie = .jf_put(p, 4L, sf::st_polygon(list(rbind(
      c(xy[4, 1] - 0.01, xy[4, 2] - 0.01), c(xy[4, 1] + 0.01, xy[4, 2] + 0.01),
      c(xy[4, 1] + 0.01, xy[4, 2] - 0.01), c(xy[4, 1] - 0.01, xy[4, 2] + 0.01),
      c(xy[4, 1] - 0.01, xy[4, 2] - 0.01))))),
    fig8 = .jf_put(p, 3L, sf::st_polygon(list(rbind(
      c(-94.07, 39.635), c(-94.06, 39.645), c(-94.05, 39.635), c(-94.05, 39.655),
      c(-94.06, 39.645), c(-94.07, 39.655), c(-94.07, 39.635))))),
    edge8 = .jf_put(p, 3L, sf::st_polygon(list(rbind(
      c(-94.11, 39.635), c(-94.10, 39.645), c(-94.09, 39.635), c(-94.09, 39.655),
      c(-94.10, 39.645), c(-94.11, 39.655), c(-94.11, 39.635))))),
    hole = {
      b <- sf::st_bbox(sf::st_geometry(p)[[5]])
      .jf_put(p, 5L, sf::st_polygon(list(.jf_sq(b[["xmin"]], b[["ymin"]], 0.02),
        rbind(c(b[["xmin"]] + 0.01, b[["ymin"]] + 0.004), c(b[["xmin"]] + 0.03, b[["ymin"]] + 0.004),
              c(b[["xmin"]] + 0.03, b[["ymin"]] + 0.016), c(b[["xmin"]] + 0.01, b[["ymin"]] + 0.016),
              c(b[["xmin"]] + 0.01, b[["ymin"]] + 0.004)))))
    })
}

.jf_assign <- function(x, cells, ...)
  assign_features_to_polygons(x, cells, polygon_id_col = "cell_id", quiet = TRUE, ...)
.jf_catch <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)
.jf_ids <- function(got) got$cell_id[order(got$row)]

# Whether s2 accepts every piece of the overlaps of two layers in degrees, as
# the largest-overlap join takes them: the rows s2 refuses repaired with
# sf::st_make_valid(), the features intersected with the polygons, and each
# piece's area taken on s2.  A figure eight is repaired into two parts that
# touch at its pinch, and where the pinch lies on the edge of a polygon of
# the other layer, s2 can compute that point twice in a row in the piece
# and then refuses it.  Whether it does follows the last bits of the
# vertices s2 computes: on macOS on arm64 (s2 1.1.11) the piece of `edge8`
# with cell 3 holds 39.644999999999996 twice and s2 refuses it, on Linux on
# x86-64 (s2 1.1.13) it holds 39.644999999999989 and then
# 39.644999999999996 and s2 accepts it, from the same repaired feature.  So
# the tests below read it on the platform they run on and pin each branch.
.jf_pieces_accepted <- function(feats, cells) spatialkit:::.with_s2({
  mend <- function(x) {
    g <- sf::st_geometry(x)
    bad <- !(sf::st_is_valid(g) %in% TRUE)
    if (any(bad)) g[bad] <- suppressWarnings(sf::st_make_valid(g[bad]))
    g
  }
  f <- mend(feats); p <- mend(cells)
  tryCatch({
    pieces <- suppressMessages(sf::st_intersection(f, p))
    all(is.finite(as.numeric(sf::st_area(pieces))))
  }, error = function(e) FALSE)
})

# Each feature's polygons of largest overlap, read without the package: the
# overlaps measured by GEOS in the equal-area EPSG:5070, on both layers
# repaired there.  The two lobes of `edge8` are mirror images across the
# edge of cells 2 and 3, and their areas there agree to 3e-14 relative, so
# areas within 1e-9 relative of the largest count as tied and every tied
# polygon is listed.
.jf_largest <- function(feats, cells) {
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)))
  f <- sf::st_make_valid(sf::st_transform(sf::st_geometry(feats), 5070))
  p <- sf::st_make_valid(sf::st_transform(sf::st_geometry(cells), 5070))
  lapply(seq_along(f), function(i) {
    a <- vapply(seq_along(p), function(j) {
      x <- suppressMessages(sf::st_intersection(f[i], p[j]))
      if (length(x)) sum(as.numeric(sf::st_area(x))) else 0
    }, numeric(1))
    cells$cell_id[a > 0 & a >= max(a) * (1 - 1e-9)]
  })
}

# Where s2 accepts the pieces, the largest-overlap join returns: each
# feature gets a polygon of its largest overlap, and for these layers the
# polygon the join by `predicate` gives (`ids`).
.jf_expect_joined <- function(got, feats, cells, ids) {
  expect_s3_class(got, "sf")
  expect_identical(nrow(got), nrow(feats))
  got_ids <- .jf_ids(got)
  expect_identical(got_ids, ids)
  lg <- .jf_largest(feats, cells)
  expect_true(all(vapply(seq_along(lg), function(i) got_ids[i] %in% lg[[i]], logical(1))))
}

test_that("a polygon feature s2 refuses stops the join by `predicate` with the features' error", {
  .jf_s2(TRUE)
  cells <- .jf_squares()
  cases <- list(list(kind = "repeat", row = 2L, fixable = TRUE),
                list(kind = "bowtie", row = 4L, fixable = FALSE),
                list(kind = "fig8",   row = 3L, fixable = TRUE),
                list(kind = "hole",   row = 5L, fixable = FALSE))
  for (cs in cases) {
    feats <- .jf_features(cs$kind)
    reason <- sf::st_is_valid(feats, reason = TRUE)[cs$row]
    expect_identical(which(!sf::st_is_valid(feats)), cs$row)
    err <- .jf_catch(.jf_assign(feats, cells, largest = FALSE))
    expect_s3_class(err, "spatialkit_features_degenerate")
    expect_false(inherits(err, "spatialkit_cells_degenerate"))
    expect_identical(err$rows, cs$row)
    expect_identical(err$n, 6L)
    expect_identical(err$reasons, unname(reason))
    expect_false(err$repaired)
    expect_s3_class(err$parent, "error")
    expect_null(conditionCall(err))
    msg <- conditionMessage(err)
    expect_match(msg, sprintf(paste0(
      "assign_features_to_polygons(): s2, the engine sf uses for geometry in longitude and ",
      "latitude, refuses 1 of the 6 features of `features_sf` (row %d; \"%s\"), so the join ",
      "cannot run in degrees."), cs$row, conditionMessage(err$parent)), fixed = TRUE)
    expect_match(msg, paste0("Transform `features_sf` to a projected CRS, such as the one ",
                             "ensure_projected() picks for it (sf::st_transform()), and assign ",
                             "again; the join then runs in that CRS."), fixed = TRUE)
    expect_no_match(msg, "polygons_sf", fixed = TRUE)
    expect_no_match(msg, "two corners", fixed = TRUE)
    expect_no_match(msg, "Pass largest = FALSE", fixed = TRUE)
    # The repair is offered where it makes s2 accept the row, and works.
    offer <- "Or repair it with sf::st_make_valid(), which s2 accepts, and assign again."
    if (cs$fixable) {
      expect_match(msg, offer, fixed = TRUE)
      expect_identical(.jf_ids(.jf_assign(sf::st_make_valid(feats), cells, largest = FALSE)), 1:6)
    } else {
      expect_no_match(msg, offer, fixed = TRUE)
    }
    # The transform it names lets the join run.
    expect_identical(.jf_ids(.jf_assign(sf::st_transform(feats, 32615), cells,
                                        largest = FALSE)), 1:6)
  }
})

test_that("features without a CRS that s2 refuses get the remedy that sets the CRS first", {
  .jf_s2(TRUE)
  feats <- sf::st_set_crs(.jf_features("repeat"), NA)
  err <- .jf_catch(.jf_assign(feats, .jf_squares(), largest = FALSE))
  expect_s3_class(err, "spatialkit_features_degenerate")
  expect_identical(err$rows, 2L)
  expect_match(conditionMessage(err), paste0(
    "`features_sf` has no CRS and was read as longitude and latitude: give it its CRS with ",
    "sf::st_set_crs(), transform it to a projected CRS with sf::st_transform(), and assign ",
    "again; the join then runs in that CRS."), fixed = TRUE)
  expect_no_match(conditionMessage(err), "Or repair", fixed = TRUE)
  got <- .jf_assign(sf::st_transform(sf::st_set_crs(feats, 4326), 32615), .jf_squares(),
                    largest = FALSE)
  expect_identical(.jf_ids(got), 1:6)
})

test_that("where s2 refuses polygons of both layers the cells' error quotes each row's own reason and names the features", {
  .jf_s2(TRUE)
  cells <- .jf_cells8()
  feats <- .jf_features("repeat")
  err <- .jf_catch(.jf_assign(feats, cells, largest = FALSE))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, 3L)
  expect_identical(err$reasons, unname(sf::st_is_valid(cells, reason = TRUE)[3L]))
  expect_identical(err$features_rows, 2L)
  expect_false(err$repaired)
  msg <- conditionMessage(err)
  expect_match(msg, sprintf("refuses 1 of the 6 polygons of `polygons_sf` (row 3: \"%s\")",
                            err$reasons), fixed = TRUE)
  # s2's message for the join is the feature's, not row 3's.
  expect_match(conditionMessage(err$parent), "Edge 1 is degenerate", fixed = TRUE)
  expect_no_match(msg, "(row 3; \"Loop 0 is not valid: Edge 1", fixed = TRUE)
  expect_match(msg, sprintf("s2's error reads \"%s\".", conditionMessage(err$parent)),
               fixed = TRUE)
  expect_match(msg, sprintf("s2 also refuses 1 of the 6 features of `features_sf` (row 2: \"%s\").",
                            sf::st_is_valid(feats, reason = TRUE)[2L]), fixed = TRUE)
  expect_match(msg, "Transform `polygons_sf` to a projected CRS", fixed = TRUE)
  # The cells in a projected CRS: the join runs there, the features moved.
  expect_identical(.jf_ids(.jf_assign(feats, sf::st_transform(cells, 32615), largest = FALSE)),
                   1:6)
})

test_that("a figure-eight cell the largest-overlap join repairs stops it with the cells' error, not the advice of `largest = FALSE`", {
  .jf_s2(TRUE)
  cells <- .jf_cells8()
  err <- NULL
  expect_warning(err <- tryCatch(.jf_assign(.jf_patches(), cells), error = function(e) e),
                 "0 of 6 feature(s) and 1 of 6 polygon(s) have invalid geometry", fixed = TRUE)
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, 3L)
  expect_true(err$repaired)
  expect_identical(err$features_rows, integer(0))
  msg <- conditionMessage(err)
  expect_match(msg, sprintf("refuses 1 of the 6 polygons of `polygons_sf` (row 3: \"%s\")",
                            err$reasons), fixed = TRUE)
  expect_match(msg, paste0("The largest-overlap join repaired it with sf::st_make_valid(), which ",
                           "s2 accepts, and s2 then refused a piece of an overlap of the two ",
                           "layers (\"Loop 0 is not valid: Edge 4 is degenerate (duplicate ",
                           "vertex)\"); the join by `predicate`, which repairs nothing, meets it ",
                           "as it arrived."), fixed = TRUE)
  expect_match(msg, "Transform `polygons_sf` to a projected CRS", fixed = TRUE)
  expect_no_match(msg, "Pass largest = FALSE", fixed = TRUE)
  expect_no_match(msg, "s2 also refuses", fixed = TRUE)
  # `largest = FALSE` meets the cell as it arrived.
  err <- .jf_catch(.jf_assign(.jf_patches(), cells, largest = FALSE))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, 3L)
  expect_false(err$repaired)
  # The cells in a projected CRS: both joins run.
  proj <- sf::st_transform(cells, 32615)
  expect_identical(.jf_ids(suppressWarnings(.jf_assign(.jf_patches(), proj))), 1:6)
  expect_identical(.jf_ids(.jf_assign(.jf_patches(), proj, largest = FALSE)), 1:6)
})

test_that("a figure-eight feature on a cell edge stops the largest-overlap join with the features' error where s2 refuses its piece, and is joined where s2 accepts it", {
  .jf_s2(TRUE)
  cells <- .jf_squares()
  feats <- .jf_features("edge8")
  err <- .jf_catch(.jf_assign(feats, cells))
  if (.jf_pieces_accepted(feats, cells)) {
    # s2 accepts the piece (Linux on x86-64): the join runs.
    .jf_expect_joined(err, feats, cells, c(1L, 2L, 2L, 4L, 5L, 6L))
  } else {
    # s2 refuses the piece (macOS on arm64): the features' error.
    expect_s3_class(err, "spatialkit_features_degenerate")
    expect_identical(err$rows, 3L)
    expect_true(err$repaired)
    msg <- conditionMessage(err)
    expect_match(msg, "refuses 1 of the 6 features of `features_sf` (row 3: \"", fixed = TRUE)
    expect_match(msg, "The largest-overlap join repaired it with sf::st_make_valid()", fixed = TRUE)
    expect_match(msg, "Transform `features_sf` to a projected CRS", fixed = TRUE)
    expect_no_match(msg, "Pass largest = FALSE", fixed = TRUE)
    expect_no_match(msg, "polygons_sf", fixed = TRUE)
  }
  # The join by `predicate` meets the feature as it arrived, which s2
  # refuses on every platform.
  err <- .jf_catch(.jf_assign(feats, cells, largest = FALSE))
  expect_s3_class(err, "spatialkit_features_degenerate")
  expect_identical(err$rows, 3L)
  expect_false(err$repaired)
  # The features in a projected CRS: both joins run.
  proj <- sf::st_transform(feats, 32615)
  expect_identical(nrow(suppressWarnings(.jf_assign(proj, cells))), 6L)
  expect_identical(.jf_ids(.jf_assign(proj, cells, largest = FALSE)), 1:6)
})

test_that("features the largest-overlap join's repair leaves refused stop it with the features' error", {
  .jf_s2(TRUE)
  cells <- .jf_squares()
  for (cs in list(list(kind = "bowtie", row = 4L), list(kind = "hole", row = 5L))) {
    feats <- .jf_features(cs$kind)
    err <- NULL
    expect_warning(err <- tryCatch(.jf_assign(feats, cells), error = function(e) e),
                   "1 of 6 feature(s) and 0 of 6 polygon(s) have invalid geometry", fixed = TRUE)
    expect_s3_class(err, "spatialkit_features_degenerate")
    expect_identical(err$rows, cs$row)
    expect_false(err$repaired)
    expect_no_match(conditionMessage(err), "Pass largest = FALSE", fixed = TRUE)
    expect_match(conditionMessage(err), "Transform `features_sf` to a projected CRS", fixed = TRUE)
    got <- suppressWarnings(.jf_assign(sf::st_transform(feats, 32615), cells))
    expect_identical(.jf_ids(got), 1:6)
  }
})

test_that("layers s2 accepts as they arrive keep the advice of `largest = FALSE`, which runs", {
  .jf_s2(TRUE)
  # The caller repaired the figure-eight cell, or the figure-eight feature:
  # s2 accepts both layers and the join by `predicate` runs.  Where s2
  # refuses a piece of an overlap (macOS on arm64) the largest-overlap join
  # stops with the advice; where it accepts the piece (Linux on x86-64) the
  # join runs.
  cells <- sf::st_make_valid(.jf_cells8())
  expect_true(all(sf::st_is_valid(cells)))
  err <- .jf_catch(.jf_assign(.jf_patches(), cells))
  if (.jf_pieces_accepted(.jf_patches(), cells)) {
    .jf_expect_joined(err, .jf_patches(), cells, 1:6)
  } else {
    expect_s3_class(err, "simpleError")
    expect_match(conditionMessage(err), paste0(
      "assign_features_to_polygons(): the largest-overlap join failed (Loop 0 is not valid: ",
      "Edge 4 is degenerate (duplicate vertex)). Pass largest = FALSE to assign every feature by ",
      "`predicate` and `tie_break` instead, or project both layers first (see ensure_projected())."),
      fixed = TRUE)
  }
  expect_identical(.jf_ids(.jf_assign(.jf_patches(), cells, largest = FALSE)), 1:6)
  feats <- sf::st_make_valid(.jf_features("edge8"))
  expect_true(all(sf::st_is_valid(feats)))
  err <- .jf_catch(.jf_assign(feats, .jf_squares()))
  if (.jf_pieces_accepted(feats, .jf_squares())) {
    .jf_expect_joined(err, feats, .jf_squares(), c(1L, 2L, 2L, 4L, 5L, 6L))
  } else {
    expect_s3_class(err, "simpleError")
    expect_match(conditionMessage(err), "Pass largest = FALSE", fixed = TRUE)
  }
  expect_identical(.jf_ids(.jf_assign(feats, .jf_squares(), largest = FALSE)), c(1L, 2L, 2L, 4L, 5L, 6L))
  # Layers s2 accepts join as before.
  expect_identical(.jf_ids(.jf_assign(.jf_points(), .jf_squares())), 1:6)
  expect_identical(.jf_ids(.jf_assign(.jf_patches(), .jf_squares())), 1:6)
  expect_identical(.jf_ids(.jf_assign(.jf_patches(), .jf_squares(), largest = FALSE)), 1:6)
})

# A 100 m square at the origin and a bowtie over it in EPSG:32615: GEOS
# gives the bowtie 0 m^2, its two lobes cancelling; repaired it covers
# 180,000 m^2.  The point at (60, 50) lies in both.
.jf_bowtie_cells <- function() {
  sq <- sf::st_polygon(list(.jf_sq(0, 0, 100)))
  bw <- sf::st_polygon(list(rbind(c(50, -200), c(650, 400), c(650, -200), c(50, 400), c(50, -200))))
  sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(sq, bw, crs = 32615))
}

test_that("the smallest-area tie-break measures a polygon whose ring crosses itself repaired, with a warning", {
  cells <- .jf_bowtie_cells()
  expect_identical(as.numeric(sf::st_area(cells))[2L], 0)
  pt <- sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_point(c(60, 50)), crs = 32615))
  got <- NULL
  expect_warning(got <- assign_features_to_polygons(pt, cells, quiet = TRUE), paste0(
    "assign_features_to_polygons(): 1 of the 2 polygon(s) whose areas the smallest-area ",
    "tie-break compares have invalid geometry (usually a self-intersecting ring), so their ",
    "areas are undefined; they were repaired with sf::st_make_valid() to measure them for ",
    "the tie-break."), fixed = TRUE)
  expect_identical(got$poly_id, 1L)
  expect_identical(attr(got, "ties")$which, 1L)
  # The geometry returned is the caller's.
  expect_identical(sf::st_geometry(got), sf::st_geometry(pt))
  # A polygon feature joined by `predicate` too.
  pg <- sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_polygon(list(.jf_sq(55, 45, 10))), crs = 32615))
  expect_warning(got <- assign_features_to_polygons(pg, cells, quiet = TRUE, largest = FALSE),
                 "whose areas the smallest-area tie-break compares", fixed = TRUE)
  expect_identical(got$poly_id, 1L)
  # The bowtie repaired by the caller, or a valid polygon in its place:
  # the same polygon, and no warning.
  fixed <- cells
  sf::st_geometry(fixed)[2] <- sf::st_geometry(sf::st_make_valid(cells[2, ]))
  expect_no_warning(got <- assign_features_to_polygons(pt, fixed, quiet = TRUE))
  expect_identical(got$poly_id, 1L)
  valid <- cells
  sf::st_geometry(valid)[2] <- sf::st_polygon(list(.jf_sq(50, -200, 600)))
  expect_no_warning(got <- assign_features_to_polygons(pt, valid, quiet = TRUE))
  expect_identical(got$poly_id, 1L)
  # With tie_break = "first" no area is taken, and no warning is raised.
  expect_no_warning(got <- assign_features_to_polygons(pt, cells, quiet = TRUE, tie_break = "first"))
  expect_identical(got$poly_id, 1L)
})

# Two overlapping 0.2 degree squares at 39N in EPSG:4326, three points (the
# second in both) and a 0.04 degree square across both.
.jf_overlap <- function() list(
  cells = sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(.jf_sq(-94, 39, 0.2))), sf::st_polygon(list(.jf_sq(-93.9, 39, 0.2))),
    crs = 4326)),
  pts = sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_point(c(-93.95, 39.1)), sf::st_point(c(-93.85, 39.1)), sf::st_point(c(-93.75, 39.1)),
    crs = 4326)),
  poly = sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_polygon(list(.jf_sq(-93.92, 39.08, 0.04))),
                                                  crs = 4326)))

test_that("with s2 off, lon/lat layers are assigned as with s2 on, and the session's setting stays", {
  ov <- .jf_overlap()
  run <- function() {
    a <- suppressMessages(assign_features_to_polygons(ov$pts, ov$cells, quiet = TRUE))
    b <- assign_features_to_polygons(ov$poly, ov$cells, quiet = TRUE)
    d <- suppressMessages(assign_features_to_polygons(ov$poly, ov$cells, quiet = TRUE,
                                                      largest = FALSE))
    list(pts = a$poly_id, pts_ties = attr(a, "ties")$which, poly = b$poly_id,
         poly_pred = d$poly_id)
  }
  .jf_s2(TRUE)
  on <- run()
  expect_identical(on$pts, c(1L, 1L, 2L))
  expect_identical(on$poly, 1L)
  suppressMessages(sf::sf_use_s2(FALSE))
  off <- run()
  expect_identical(off, on)
  expect_false(sf::sf_use_s2())
  # tie_break = "first" took no area and ran before; it still does.
  got <- suppressMessages(assign_features_to_polygons(ov$pts, ov$cells, quiet = TRUE,
                                                      tie_break = "first"))
  expect_identical(got$poly_id, c(1L, 1L, 2L))
  # Projected layers are measured by GEOS, as before.
  got <- assign_features_to_polygons(sf::st_transform(ov$pts, 32615),
                                     sf::st_transform(ov$cells, 32615), quiet = TRUE)
  expect_identical(got$poly_id, c(1L, 2L, 2L))
})

test_that("with s2 off, a tied polygon s2 refuses is measured repaired, or named when the repair leaves it refused", {
  ov <- .jf_overlap()
  .jf_s2(FALSE)
  # Cell 2 with a vertex listed twice, which GEOS accepts and s2 refuses:
  # measured on s2's repair, with the warning.
  m <- .jf_sq(-93.9, 39, 0.2)
  rep2 <- .jf_put(ov$cells, 2L, sf::st_polygon(list(rbind(m[1:2, ], m[2:5, ]))))
  got <- NULL
  expect_warning(got <- suppressMessages(assign_features_to_polygons(ov$pts, rep2, quiet = TRUE)),
                 "whose areas the smallest-area tie-break compares", fixed = TRUE)
  expect_identical(got$poly_id, c(1L, 1L, 2L))
  expect_false(sf::sf_use_s2())
  # Cell 2 a bowtie, which s2 refuses even repaired: the point in cell 1
  # and in the bowtie's western lobe has no area to compare.
  bow2 <- .jf_put(ov$cells, 2L, sf::st_polygon(list(rbind(
    c(-93.9, 39), c(-93.7, 39.2), c(-93.7, 39), c(-93.9, 39.2), c(-93.9, 39)))))
  err <- .jf_catch(suppressMessages(assign_features_to_polygons(ov$pts[2, ], bow2, quiet = TRUE)))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, 2L)
  expect_identical(err$n, 2L)
  expect_match(conditionMessage(err), paste0(
    "refuses 1 of the 2 polygons of `polygons_sf` (row 2; \"Loop 0 is not valid: Edge 0 ",
    "crosses edge 2\"), so the smallest-area tie-break cannot run in degrees."), fixed = TRUE)
  expect_match(conditionMessage(err), "Transform `polygons_sf` to a projected CRS", fixed = TRUE)
  # s2 is not the engine sf uses in this session: the message says so.
  expect_true(err$s2_forced)
  expect_match(conditionMessage(err), paste0(
    "assign_features_to_polygons(): s2, which the package uses here in longitude and latitude ",
    "whatever sf::sf_use_s2() says, refuses"), fixed = TRUE)
  expect_no_match(conditionMessage(err), "the engine sf uses", fixed = TRUE)
  expect_false(sf::sf_use_s2())
  # In a projected CRS GEOS repairs the bowtie into its two lobes, each a
  # quarter of cell 1, so the point goes to the bowtie as the smaller.
  got <- suppressWarnings(assign_features_to_polygons(ov$pts[2, ], sf::st_transform(bow2, 32615),
                                                      quiet = TRUE))
  expect_identical(got$poly_id, 2L)
})

# Cell 5 with a hole that crosses its outer ring, which s2 refuses even
# repaired.
.jf_hole_cells <- function() {
  cc <- .jf_squares(); b <- sf::st_bbox(sf::st_geometry(cc)[[5]])
  .jf_put(cc, 5L, sf::st_polygon(list(.jf_sq(b[["xmin"]], b[["ymin"]], 0.1),
    rbind(c(b[["xmin"]] + 0.05, b[["ymin"]] + 0.02), c(b[["xmin"]] + 0.15, b[["ymin"]] + 0.02),
          c(b[["xmin"]] + 0.15, b[["ymin"]] + 0.08), c(b[["xmin"]] + 0.05, b[["ymin"]] + 0.08),
          c(b[["xmin"]] + 0.05, b[["ymin"]] + 0.02)))))
}

test_that("with s2 off, the largest-overlap join's error for a polygon s2 refuses offers `largest = FALSE`, which runs with GEOS", {
  .jf_s2(FALSE)
  offer <- paste0("Or pass largest = FALSE to assign every feature by `predicate` and ",
                  "`tie_break` instead; with sf::sf_use_s2(FALSE), as in this session, that ",
                  "join runs with GEOS.")
  engine <- paste0("assign_features_to_polygons(): s2, which the package uses here in ",
                   "longitude and latitude whatever sf::sf_use_s2() says, refuses")
  cases <- list(
    list(feats = .jf_patches(), cells = .jf_cells8(), class = "spatialkit_cells_degenerate",
         rows = 3L, repaired = TRUE, ids = 1:6),
    list(feats = .jf_patches(), cells = .jf_hole_cells(), class = "spatialkit_cells_degenerate",
         rows = 5L, repaired = FALSE, ids = 1:6),
    list(feats = .jf_features("hole"), cells = .jf_squares(),
         class = "spatialkit_features_degenerate", rows = 5L, repaired = FALSE, ids = 1:6),
    list(feats = .jf_features("edge8"), cells = .jf_squares(),
         class = "spatialkit_features_degenerate", rows = 3L, repaired = TRUE,
         ids = c(1L, 2L, 2L, 4L, 5L, 6L)))
  # Whether s2 refuses the piece of the repaired figure eight `edge8`
  # follows the platform (.jf_pieces_accepted()); where it accepts it, the
  # largest-overlap join runs, on s2 in degrees whatever the session says.
  # The other three cases stop on every platform.
  for (cs in cases) {
    err <- .jf_catch(.jf_assign(cs$feats, cs$cells))
    if (cs$repaired && .jf_pieces_accepted(cs$feats, cs$cells)) {
      .jf_expect_joined(err, cs$feats, cs$cells, cs$ids)
    } else {
      expect_s3_class(err, cs$class)
      expect_identical(err$rows, cs$rows)
      expect_identical(err$repaired, cs$repaired)
      expect_true(err$s2_forced)
      msg <- conditionMessage(err)
      expect_match(msg, engine, fixed = TRUE)
      expect_match(msg, "so the largest-overlap join cannot run in degrees.", fixed = TRUE)
      expect_match(msg, offer, fixed = TRUE)
      expect_no_match(msg, "meets it as it arrived", fixed = TRUE)
      expect_no_match(msg, "the engine sf uses", fixed = TRUE)
    }
    expect_false(sf::sf_use_s2())
    # The join it offers runs, with GEOS.
    got <- suppressMessages(suppressWarnings(.jf_assign(cs$feats, cs$cells, largest = FALSE)))
    expect_identical(.jf_ids(got), cs$ids)
    expect_false(sf::sf_use_s2())
  }
  # With s2 on the join by `predicate` meets the same polygons, and the
  # offer is not made.
  .jf_s2(TRUE)
  for (cs in cases) {
    err <- .jf_catch(.jf_assign(cs$feats, cs$cells))
    if (cs$repaired && .jf_pieces_accepted(cs$feats, cs$cells)) {
      .jf_expect_joined(err, cs$feats, cs$cells, cs$ids)
      next
    }
    expect_s3_class(err, cs$class)
    expect_false(err$s2_forced)
    expect_no_match(conditionMessage(err), "largest = FALSE", fixed = TRUE)
    expect_match(conditionMessage(err), "s2, the engine sf uses for geometry in longitude",
                 fixed = TRUE)
  }
})

test_that("with s2 off, the largest-overlap join checks and repairs on s2 a feature only s2 refuses", {
  .jf_s2(FALSE)
  # A vertex listed twice: GEOS accepts it, s2 refuses it and accepts its
  # repair.
  feats <- .jf_features("repeat")
  expect_true(all(sf::st_is_valid(feats)))
  got <- NULL
  expect_warning(got <- .jf_assign(feats, .jf_squares()),
                 "1 of 6 feature(s) and 0 of 6 polygon(s) have invalid geometry", fixed = TRUE)
  expect_identical(.jf_ids(got), 1:6)
  expect_identical(sf::st_geometry(got)[order(got$row)], sf::st_geometry(feats))
  expect_false(sf::sf_use_s2())
})

test_that("a feature s2 refuses even repaired, beside a figure-eight cell, stops the largest-overlap join with the features' error", {
  .jf_s2(TRUE)
  cells <- .jf_cells8()
  for (cs in list(list(kind = "hole", row = 5L), list(kind = "bowtie", row = 4L))) {
    feats <- .jf_features(cs$kind)
    err <- NULL
    expect_warning(err <- tryCatch(.jf_assign(feats, cells), error = function(e) e),
                   "1 of 6 feature(s) and 1 of 6 polygon(s) have invalid geometry", fixed = TRUE)
    expect_s3_class(err, "spatialkit_features_degenerate")
    expect_identical(err$rows, cs$row)
    expect_false(err$repaired)
    expect_identical(err$reasons, unname(sf::st_is_valid(feats, reason = TRUE)[cs$row]))
    msg <- conditionMessage(err)
    expect_match(msg, sprintf("refuses 1 of the 6 features of `features_sf` (row %d; \"%s\")",
                              cs$row, conditionMessage(err$parent)), fixed = TRUE)
    expect_no_match(msg, "piece of an overlap", fixed = TRUE)
    expect_no_match(msg, "polygons_sf", fixed = TRUE)
    expect_match(msg, "Transform `features_sf` to a projected CRS", fixed = TRUE)
    # The transform it names lets the join run, the cell repaired as before.
    got <- suppressWarnings(.jf_assign(sf::st_transform(feats, 32615), cells))
    expect_identical(.jf_ids(got), 1:6)
  }
})

test_that("seven refused rows are quoted six at a time with the count of the rest", {
  .jf_s2(TRUE)
  xy <- sf::st_coordinates(.jf_points())
  g <- lapply(c(1:6, 1L), function(i) {
    m <- .jf_sq(xy[i, 1] - 0.01, xy[i, 2] - 0.01, 0.02)
    sf::st_polygon(list(rbind(m[1:2, ], m[2:5, ])))
  })
  feats <- sf::st_sf(row = 1:7, geometry = sf::st_sfc(g, crs = 4326))
  err <- .jf_catch(.jf_assign(feats, .jf_cells8(), largest = FALSE))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$features_rows, 1:7)
  why <- sf::st_is_valid(feats, reason = TRUE)
  expect_match(conditionMessage(err), sprintf(
    "s2 also refuses 7 of the 7 features of `features_sf` (%s; and 1 more).",
    paste(sprintf("row %d: \"%s\"", 1:6, why[1:6]), collapse = "; ")), fixed = TRUE)
})

test_that("cell_values_at() with s2 off names s2 as the engine the package uses there", {
  ov <- .jf_overlap()
  .jf_s2(FALSE)
  # Cell 2 a bowtie, which s2 refuses even repaired; the location lies in
  # cell 1 and in the bowtie's western lobe.
  bow2 <- .jf_put(ov$cells, 2L, sf::st_polygon(list(rbind(
    c(-93.9, 39), c(-93.7, 39.2), c(-93.7, 39), c(-93.9, 39.2), c(-93.9, 39)))))
  bow2$resp_mean_z <- c(1.5, 3)
  asg <- sf::st_sf(poly_id = c(1L, 1L, 2L), z = c(1, 2, 3), geometry = sf::st_geometry(ov$pts))
  err <- .jf_catch(suppressMessages(cell_values_at(ov$pts[2, ], bow2, asg, "z",
                                                   id_col = "poly_id")))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, 2L)
  expect_true(err$s2_forced)
  expect_match(conditionMessage(err), paste0(
    "cell_values_at(): s2, which the package uses here in longitude and latitude whatever ",
    "sf::sf_use_s2() says, refuses 1 of the 2 polygons of `values` (row 2;"), fixed = TRUE)
  expect_false(sf::sf_use_s2())
})

test_that("the largest-overlap join does not advise `largest = FALSE` for an error that join meets too", {
  real_st_join <- sf::st_join
  local_mocked_bindings(
    st_join = function(x, y, ..., largest = FALSE) {
      if (isTRUE(largest)) stop("package lwgeom required, please install it first")
      real_st_join(x, y, ..., largest = largest)
    },
    .package = "sf"
  )
  .jf_s2(TRUE)
  ov <- .jf_overlap()
  err <- .jf_catch(assign_features_to_polygons(ov$poly, ov$cells, quiet = TRUE))
  expect_identical(conditionMessage(err), paste0(
    "assign_features_to_polygons(): the largest-overlap join failed (package lwgeom required, ",
    "please install it first). Project both layers first (see ensure_projected())."))
  err <- .jf_catch(assign_features_to_polygons(sf::st_transform(ov$poly, 32615),
                                               sf::st_transform(ov$cells, 32615), quiet = TRUE))
  expect_identical(conditionMessage(err), paste0(
    "assign_features_to_polygons(): the largest-overlap join failed (package lwgeom required, ",
    "please install it first)."))
})

test_that("?assign_features_to_polygons says what the joins do with these layers", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "assign_features_to_polygons.Rd")
  expect_match(rd, paste0("A polygon of \\code{features_sf} that s2 refuses stops the join by ",
                          "\\code{predicate} with an error of class ",
                          "\\code{\"spatialkit_features_degenerate\"}"), fixed = TRUE)
  expect_match(rd, paste0("the join stops with the error for the layer whose polygons s2 refused ",
                          "as they arrived; with s2 on (sf's default), the join by ",
                          "\\code{predicate} meets those polygons too."),
               fixed = TRUE)
  expect_match(rd, paste0("In longitude and latitude with s2 on (sf's default), where s2 refuses ",
                          "a polygon of either layer as it arrived, the join by \\code{predicate} ",
                          "meets that refusal too, and the error names the polygons and the ",
                          "transform instead"),
               fixed = TRUE)
  expect_match(rd, paste0("with s2 off the join by \\code{predicate} runs with GEOS, and the ",
                          "error that names the polygons also offers \\code{largest = FALSE}."),
               fixed = TRUE)
  expect_no_match(rd, "which the join by \\code{predicate} meets too.", fixed = TRUE)
  expect_match(rd, "Areas in longitude and latitude are measured on the sphere with s2, whatever ",
               fixed = TRUE)
  expect_match(rd, "is measured on a copy repaired with \\code{sf::st_make_valid()}, with a warning",
               fixed = TRUE)
})

test_that("the help says whether s2 refuses a piece of a repaired figure eight can follow the platform", {
  root <- testthat::test_path("..", "..")
  # The page's source where it is here, else the installed page as Rd text.
  rd <- if (file.exists(file.path(root, "man", "assign_features_to_polygons.Rd"))) {
    paste(readLines(file.path(root, "man", "assign_features_to_polygons.Rd"), warn = FALSE),
          collapse = "\n")
  } else {
    db <- tools::Rd_db("spatialkit")
    if (is.null(db[["assign_features_to_polygons.Rd"]])) skip("the help is not here to read")
    paste(as.character(db[["assign_features_to_polygons.Rd"]]), collapse = "")
  }
  rd <- gsub("[[:space:]]+", " ", rd)
  expect_match(rd, paste0(
    "Whether s2 refuses such a piece can differ between platforms and s2 versions. ",
    "Where s2 accepts the piece, the join runs and gives each feature the polygon of its ",
    "largest overlap."), fixed = TRUE)
  # Most figure eights pinched on that edge get the same outcome on both
  # platforms (some refused, some accepted on each); only some differ, so
  # the help does not say that every one does.
  expect_no_match(rd, "gives a piece that s2 refuses", fixed = TRUE)
})
