# ===========================================================================
# ensure_stable_poly_id() on a feature that is valid in its own CRS and that
# s2 refuses in lon/lat, because one of its rings crosses itself there.
#
# The sort key is measured on a lon/lat copy with s2.  st_transform() moves
# vertices only, and s2 joins them with great-circle arcs, so a long straight
# edge that passes about a metre from another vertex of the same ring can end
# up on the other side of that vertex.  The repair of the lon/lat copy
# (st_make_valid(), which is s2's rebuild there) does not split crossing
# edges, and st_centroid() and st_area() then stopped the function with
# "Loop 1 is not valid: Edge 36 crosses edge 52": one of 291 Voronoi cells of
# Texas clipped to the state outline, so build_tessellation(method =
# "voronoi") returned nothing for the whole state.
#
# The fixture holds that cell and the six cells touching it, with the
# coordinates create_voronoi_polygons() produced.  Whether s2 refuses the cell
# depends on the s2 and PROJ in use, so each test that needs the refusal
# checks for it first and is skipped, saying so, when s2 takes the cell.
# ===========================================================================

.sce_bad_cell <- 136L

.sce_cells <- function() {
  d <- utils::read.delim(
    testthat::test_path("fixtures", "texas-voronoi-crossing-edge.tsv"),
    comment.char = "#", stringsAsFactors = FALSE)
  sf::st_sf(src = d$cell, geometry = sf::st_as_sfc(d$wkt, crs = 5070))
}

# Evaluate `code` with sf_use_s2() set to `on`, restoring the session's value.
.sce_with_s2 <- function(on, code) {
  was <- suppressMessages(sf::sf_use_s2(on))
  on.exit(suppressMessages(sf::sf_use_s2(was)), add = TRUE)
  force(code)
}

# The lon/lat copy the key is measured on, as the function has always made
# it: repaired in the layer's own CRS, transformed, repaired again.
.sce_sort_copy <- function(lyr) {
  .sce_with_s2(TRUE, {
    own <- suppressWarnings(sf::st_make_valid(sf::st_geometry(lyr)))
    suppressWarnings(sf::st_make_valid(sf::st_transform(own, 4326)))
  })
}

# Positions s2 refuses, asked one feature at a time: deliberately not the
# package's own (bisecting) search, so that the two can be compared.
.sce_refused <- function(g) {
  which(vapply(seq_along(g), function(i)
    inherits(tryCatch(sf::st_as_s2(g[i]), error = function(e) e), "error"),
    logical(1)))
}

# The rank the unguarded key gives each feature of a lon/lat copy s2 accepts
# throughout: centroid x, centroid y and area, rounded as the function rounds
# them, then position.
.sce_old_rank <- function(ll) {
  .sce_with_s2(TRUE, {
    xy <- sf::st_coordinates(sf::st_centroid(ll))
    a  <- as.numeric(sf::st_area(ll))
    r  <- integer(length(ll))
    r[order(round(xy[, 1], 7L), round(xy[, 2], 7L), signif(a, 9L),
            seq_along(ll))] <- seq_along(ll)
    r
  })
}

# Run the function, keeping its result, every warning it raises and every
# line it logs.
.sce_run <- function(lyr, ...) {
  raised <- character(0)
  logged <- capture_spatialkit_log(
    out <- withCallingHandlers(
      ensure_stable_poly_id(lyr, ...),
      warning = function(w) {
        raised <<- c(raised, conditionMessage(w))
        invokeRestart("muffleWarning")
      }))
  list(out = out, warnings = raised, log = logged)
}

# The result alone, with the warning and its log line kept off the console.
.sce_ids <- function(lyr, ...) .sce_run(lyr, ...)$out

# The precondition of every test below.  The layer is valid where it was
# built; after the transform and the repair the function has always made, s2
# refuses the feature at `bad` and no other.
.sce_expect_refused <- function(lyr, bad) {
  expect_true(all(sf::st_is_valid(lyr)))
  refused <- .sce_refused(.sce_sort_copy(lyr))
  if (!(bad %in% refused))
    skip(sprintf(paste0(
      "the s2 in use (s2 %s, sf %s, PROJ %s) accepts this feature after the ",
      "transform to lon/lat and the first repair, so the guard has nothing ",
      "to do here"),
      tryCatch(as.character(utils::packageVersion("s2")),
               error = function(e) "unknown"),
      as.character(utils::packageVersion("sf")),
      sf::sf_extSoftVersion()[["PROJ"]]))
  expect_identical(refused, bad)
}


# ---------------------------------------------------------------------------
# The Texas cell
# ---------------------------------------------------------------------------

test_that("a cell whose ring crosses itself in lon/lat gets its ID, and its neighbours keep theirs", {
  cells <- .sce_cells()
  expect_equal(nrow(cells), 7L)
  bad <- which(cells$src == .sce_bad_cell)
  .sce_expect_refused(cells, bad)

  # The failure as it was: both measurements stop on the repaired copy.
  ll <- .sce_sort_copy(cells)
  .sce_with_s2(TRUE, {
    expect_error(sf::st_centroid(ll), "crosses edge")
    expect_error(sf::st_area(ll), "crosses edge")
  })

  res <- .sce_run(cells)
  out <- res$out
  expect_s3_class(out, "sf")
  expect_equal(nrow(out), nrow(cells))
  expect_equal(out$poly_id, seq_len(nrow(cells)))     # complete, unique, in order
  expect_setequal(out$src, cells$src)

  # It says what it did, once: one R warning and one line in the log.
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "s2 refused 1 of 7 feature\\(s\\)")
  expect_match(res$warnings, "crosses edge")
  expect_match(res$warnings, "second repair of that copy with crossing edges split: 1\\.")
  expect_match(res$warnings, "plane geometry in the layer's own CRS instead: 0\\.")
  expect_match(res$warnings, "close to the spherical ones but not exactly equal")
  expect_equal(sum(grepl("s2 refused", res$log)), 1L)

  # The six neighbours: the unguarded key, computed here with sf alone, ranks
  # them 1 to 6.  With the seventh cell in the layer each keeps that rank,
  # moved up by one where it sorts after the newcomer.
  old <- .sce_old_rank(ll[-bad])
  id  <- out$poly_id[match(cells$src, out$src)]
  expect_equal(id[-bad], old + as.integer(old >= id[bad]))
  expect_equal(out$src, c(113L, 118L, 136L, 129L, 146L, 161L, 153L))

  # What comes back is the caller's geometry (repaired in its own CRS, as on
  # any call), never the split copy: the cell still has its two parts, where
  # the split copy has four.
  expect_equal(sf::st_crs(out), sf::st_crs(cells))
  own <- suppressWarnings(sf::st_make_valid(sf::st_geometry(cells)))
  expect_equal(sf::st_geometry(out), own[match(out$src, cells$src)])
  expect_equal(length(sf::st_cast(sf::st_geometry(out)[out$src == .sce_bad_cell],
                                  "POLYGON")), 2L)
})


test_that("the guarded IDs do not depend on row order, on the call or on sf_use_s2()", {
  cells <- .sce_cells()
  .sce_expect_refused(cells, which(cells$src == .sce_bad_cell))

  ref <- .sce_ids(cells)
  expect_identical(.sce_ids(cells), ref)

  set.seed(4)
  for (i in 1:3) {
    shuffled <- cells[sample(nrow(cells)), , drop = FALSE]
    expect_equal(.sce_ids(shuffled)$src, ref$src)
  }

  # s2 is switched on for the key and handed back as it was found, on the
  # guarded route as on any other.
  off <- .sce_with_s2(FALSE, {
    o <- .sce_ids(cells)
    expect_false(sf::sf_use_s2())
    o
  })
  expect_equal(off$src, ref$src)
  .sce_with_s2(TRUE, {
    .sce_ids(cells)
    expect_true(sf::sf_use_s2())
  })
})


test_that("the same cells get the same IDs whichever projection they arrive in", {
  cells <- .sce_cells()
  .sce_expect_refused(cells, which(cells$src == .sce_bad_cell))
  ref <- .sce_ids(cells)

  # Three different routes to the same IDs.  In Texas Centric Albers
  # (EPSG:3083) the cell is valid and s2 refuses its lon/lat copy, as in
  # EPSG:5070.  In UTM zone 14N the straight edge runs on the other side of
  # the coastline vertex, so the cell is invalid on arrival, the repair in
  # its own CRS mends it, and s2 has nothing to refuse.  In lon/lat it
  # arrives already crossing itself.
  for (crs in c(3083, 32614, 4326)) {
    moved <- sf::st_transform(cells, crs)
    got   <- .sce_ids(moved)
    expect_equal(got$poly_id, seq_len(nrow(cells)), info = paste("EPSG", crs))
    expect_equal(got$src, ref$src, info = paste("EPSG", crs))
  }
})


test_that("a feature s2 refuses after the second repair is measured in the plane", {
  cells <- .sce_cells()
  bad   <- which(cells$src == .sce_bad_cell)
  .sce_expect_refused(cells, bad)
  ref <- .sce_ids(cells)

  # make_valid = FALSE asks for no repair, so the second one is not tried
  # either and the refused feature goes straight to the plane.  Nothing is
  # repaired on the way out: the geometry is the input's, to the last digit.
  raw <- .sce_run(cells, make_valid = FALSE)
  expect_length(raw$warnings, 1L)
  expect_match(raw$warnings, "crossing edges split: 0 \\(none is tried with make_valid = FALSE\\)\\.")
  n_refused <- as.integer(sub(".*s2 refused ([0-9]+) of 7 .*", "\\1", raw$warnings))
  expect_gte(n_refused, 1L)
  expect_match(raw$warnings, sprintf("own CRS instead: %d\\.", n_refused))
  expect_equal(raw$out$poly_id, seq_len(nrow(cells)))
  expect_equal(raw$out$src, ref$src)
  expect_identical(sf::st_as_text(sf::st_geometry(raw$out), digits = 17),
                   sf::st_as_text(sf::st_geometry(cells)[match(raw$out$src, cells$src)],
                                  digits = 17))

  # With repairs on, the plane is where a feature goes when the second
  # repair is of no use.  Simulated here by a second repair that gives up.
  local_mocked_bindings(.split_crossing_edges = function(g1) NULL,
                        .package = "spatialkit")
  flat <- .sce_run(cells)
  expect_length(flat$warnings, 1L)
  expect_match(flat$warnings, "s2 refused 1 of 7 feature\\(s\\)")
  expect_match(flat$warnings, "second repair of that copy with crossing edges split: 0\\.")
  expect_match(flat$warnings, "plane geometry in the layer's own CRS instead: 1\\.")
  expect_equal(flat$out$poly_id, seq_len(nrow(cells)))
  expect_equal(flat$out$src, ref$src)

  # A layer that arrives in lon/lat has no plane of its own to fall back to
  # but its degrees; the plane route must not send it back to s2 (the call
  # that stopped) or to lwgeom (not a dependency) to measure it.
  ll_in <- .sce_run(sf::st_transform(cells, 4326))
  expect_length(ll_in$warnings, 1L)
  expect_match(ll_in$warnings, "plane geometry in the layer's own CRS instead: 1\\.")
  expect_equal(ll_in$out$src, ref$src)
})


test_that("features s2 accepts keep exactly the key they had, and a refused one gets a key close to it", {
  cells <- .sce_cells()
  bad   <- which(cells$src == .sce_bad_cell)
  .sce_expect_refused(cells, bad)
  keyed <- spatialkit:::.sort_key_with_refusals

  .sce_with_s2(TRUE, {
    own     <- sf::st_sf(geometry = suppressWarnings(sf::st_make_valid(sf::st_geometry(cells))))
    ll      <- .sce_sort_copy(cells)
    sort_sf <- sf::st_sf(geometry = ll)
    # What ensure_stable_poly_id() holds when it calls for help: two errors.
    pts_err  <- tryCatch(sf::st_centroid(ll), error = function(e) e)
    area_err <- tryCatch(sf::st_area(ll), error = function(e) e)
    expect_s3_class(pts_err, "error")
    expect_s3_class(area_err, "error")

    logged <- capture_spatialkit_log(
      key <- suppressWarnings(keyed(sort_sf, own, "centroid", pts_err, area_err)))
    expect_equal(sum(grepl("s2 refused 1 of 7", logged)), 1L)
    expect_s3_class(key$points, "sfc_POINT")
    expect_equal(sf::st_crs(key$points), sf::st_crs(4326))
    expect_length(key$points, nrow(cells))
    expect_length(key$area, nrow(cells))

    # The six: bit for bit what the unguarded calls give them.
    expect_identical(unname(sf::st_coordinates(key$points[-bad])),
                     unname(sf::st_coordinates(sf::st_centroid(ll[-bad]))))
    expect_identical(key$area[-bad], as.numeric(sf::st_area(ll[-bad])))

    # The seventh, by the second repair: within metres of the centroid taken
    # in the plane (20 m here; the two differ by 3 to 23 m for the six
    # neighbours), and the cell's area, not a lobe of it and not the rest of
    # the sphere.  EPSG:5070 is an equal-area projection, so the planar area
    # is the ellipsoid's, and the sphere's is within a fifth of a percent.
    planar_pt   <- sf::st_transform(sf::st_centroid(sf::st_geometry(own)[bad]), 4326)
    planar_area <- as.numeric(sf::st_area(sf::st_geometry(own)[bad]))
    expect_lt(as.numeric(sf::st_distance(key$points[bad], planar_pt)), 50)
    expect_equal(key$area[bad] / planar_area, 1, tolerance = 0.003)

    # The plane route gives exactly those planar numbers, and still leaves
    # the six alone.
    local_mocked_bindings(.split_crossing_edges = function(g1) NULL,
                          .package = "spatialkit")
    capture_spatialkit_log(
      flat <- suppressWarnings(keyed(sort_sf, own, "centroid", pts_err, area_err)))
    expect_equal(unname(sf::st_coordinates(flat$points[bad])),
                 unname(sf::st_coordinates(planar_pt)))
    expect_equal(flat$area[bad], planar_area)
    expect_identical(unname(sf::st_coordinates(flat$points[-bad])),
                     unname(sf::st_coordinates(key$points[-bad])))
    expect_identical(flat$area[-bad], key$area[-bad])
  })
})


test_that("the surface-point and box-centre keys are guarded too", {
  cells <- .sce_cells()
  bad   <- which(cells$src == .sce_bad_cell)
  .sce_expect_refused(cells, bad)

  # Neither point is taken with s2, but the area, the third key of every
  # method, is: both methods stopped in st_area() on this layer.
  for (m in c("surface_point", "bbox_center")) {
    res <- .sce_run(cells, method = m)
    expect_equal(sum(grepl("s2 refused 1 of 7", res$warnings)), 1L, info = m)
    expect_equal(sum(grepl("s2 refused", res$log)), 1L, info = m)
    expect_equal(res$out$poly_id, seq_len(nrow(cells)), info = m)
    expect_setequal(res$out$src, cells$src)

    # The six neighbours on their own need no guard, and within the seven
    # they keep the order they have on their own.
    alone <- .sce_run(cells[-bad, , drop = FALSE], method = m)
    expect_false(any(grepl("s2 refused", alone$warnings)), info = m)
    expect_equal(res$out$src[res$out$src != .sce_bad_cell], alone$out$src,
                 info = m)
  }
})


# ---------------------------------------------------------------------------
# The helpers
# ---------------------------------------------------------------------------

test_that(".s2_refused() finds every refused feature, wherever it sits", {
  cells <- .sce_cells()
  bad   <- which(cells$src == .sce_bad_cell)
  .sce_expect_refused(cells, bad)
  refused <- spatialkit:::.s2_refused
  ll <- .sce_sort_copy(cells)

  expect_identical(refused(ll), bad)
  expect_identical(refused(ll[-bad]), integer(0))
  expect_identical(refused(ll[integer(0)]), integer(0))
  # First, last, adjacent and alone, against the one-at-a-time search.
  for (pick in list(c(bad, 1L, 2L, 3L, 5L, 6L, 7L),
                    c(1L, 2L, 3L, 5L, 6L, 7L, bad),
                    c(bad, bad, 1L, 2L, bad, 3L, 5L, 6L, bad, bad),
                    bad, c(bad, bad))) {
    expect_identical(refused(ll[pick]), .sce_refused(ll[pick]))
    expect_identical(refused(ll[pick]), which(pick == bad))
  }
})


test_that("an error that is not s2's is raised again as it was", {
  cells <- .sce_cells()
  bad   <- which(cells$src == .sce_bad_cell)
  keyed <- spatialkit:::.sort_key_with_refusals
  boom  <- simpleError("not from s2")

  # A sort copy that is not lon/lat never goes to s2.
  expect_error(keyed(cells, cells, "centroid", boom, boom), "not from s2")
  # Nor does a lon/lat one in which s2 refuses nothing have anything to guard.
  ok <- sf::st_sf(geometry = .sce_sort_copy(cells[-bad, , drop = FALSE]))
  expect_error(keyed(ok, cells[-bad, , drop = FALSE], "centroid", boom, boom),
               "not from s2")
  # The surface point is GEOS's, so its error is not one to work around
  # even when s2 does refuse a feature.
  ll <- sf::st_sf(geometry = .sce_sort_copy(cells))
  expect_error(keyed(ll, cells, "surface_point", boom, boom), "not from s2")
})


# ---------------------------------------------------------------------------
# The same thing built from scratch
# ---------------------------------------------------------------------------

test_that("a long edge a metre from a vertex of its own ring, built from scratch", {
  # Two triangles joined by a neck one metre wide under a 40 km east-west
  # edge, in EPSG:5070 near 27N.  There the great circle through the edge's
  # ends runs south of the straight edge, by between 4 and 6 m at
  # mid-length (a neck of 4 m is refused, one of 6 m is not), so in lon/lat
  # the vertex at the neck pokes through the edge above it.
  x0 <- -170000; y0 <- 420000; len <- 40000
  bow <- sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + len, y0), c(x0 + len, y0 - 10000),
    c(x0 + len / 2, y0 - 1), c(x0, y0 - 10000), c(x0, y0))))
  sq <- function(x) sf::st_polygon(list(rbind(
    c(x, y0), c(x + 10000, y0), c(x + 10000, y0 - 10000), c(x, y0 - 10000),
    c(x, y0))))
  lyr <- sf::st_sf(src = c("east", "bow", "west"),
                   geometry = sf::st_sfc(sq(x0 + len + 20000), bow,
                                         sq(x0 - 30000), crs = 5070))
  .sce_expect_refused(lyr, 2L)

  res <- .sce_run(lyr)
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "s2 refused 1 of 3 feature\\(s\\)")
  expect_match(res$warnings, "crossing edges split: 1\\.")
  expect_equal(res$out$poly_id, 1:3)
  expect_equal(res$out$src, c("west", "bow", "east"))
  # The two squares in the order the unguarded key gives them.
  expect_equal(.sce_old_rank(.sce_sort_copy(lyr)[-2L]), c(2L, 1L))
  # The bow comes back as it went in: one ring of six points.
  expect_equal(nrow(sf::st_coordinates(res$out[res$out$src == "bow", ])), 6L)

  # The plane route agrees.
  flat <- .sce_run(lyr, make_valid = FALSE)
  expect_match(flat$warnings, "own CRS instead: 1\\.")
  expect_equal(flat$out$src, res$out$src)
})
