# tests/testthat/test-fix-stable-ids-unplaced.R
# ---------------------------------------------------------------------------
# ensure_stable_poly_id() on rows its sort copy cannot place or measure, and
# on rows its repair leaves with no polygon.
#
# 1. PROJ returns a row of a projected layer it cannot move to the sort CRS
#    (EPSG:4326 by default) as an empty row, or a ring without the vertices
#    it cannot place, instead of stopping.  A UTM layer shifted by 1e8, 1e12
#    or 1e15 had a NaN key for every row and was numbered by the row order.
#    Such a move now takes the documented fallback: the whole layer is
#    keyed in its own CRS (.rows_unplaced()), on keys rounded to a step that
#    grows with the magnitude of the coordinates.
# 2. A feature whose point or area the sort copy cannot give for a reason
#    other than s2's refusal (GEOS stops on a ring with fewer than three
#    distinct vertices) is keyed by the centre of its box with area 0,
#    where GEOS's error ended the call.
# 3. With make_valid = TRUE an empty polygon came back as LINESTRING EMPTY,
#    and a ring with no area as a GEOMETRYCOLLECTION, without a word.  An
#    empty row is not repaired now; a row the repair leaves with no polygon
#    comes back as given, with a warning of class
#    "spatialkit_cells_degenerate" that names its rows.
# ---------------------------------------------------------------------------

.su_q <- function(x) suppressWarnings(suppressMessages(x))

# Twelve squares in UTM 17N, 150 m to 400 m on a side, each in a column of
# its own 350 m from the next, so the order by longitude at the origin and
# by easting at an offset is the same.  Coordinates are multiples of 0.5 m,
# so a shift by 1e15 is exact.  The rows are not in the sort order.
.su_layer <- function(off = 0) {
  x0 <- 6e5 + 350 * (0:11)
  y0 <- 3.99e6 + c(0, 1200, 400, 2000, 800, 0, 1600, 2400, 300, 1000, 2200, 600)
  s  <- c(150, 300, 200, 400, 250, 180, 320, 220, 260, 340, 160, 280)
  geoms <- lapply(seq_along(x0), function(i) {
    x <- x0[i] + off; y <- y0[i] + off
    sf::st_polygon(list(rbind(c(x, y), c(x + s[i], y), c(x + s[i], y + s[i]),
                              c(x, y + s[i]), c(x, y))))
  })
  ord <- c(7, 2, 11, 4, 9, 1, 12, 5, 3, 10, 6, 8)
  sf::st_sf(name = sprintf("p%02d", ord), geometry = sf::st_sfc(geoms[ord], crs = 32617))
}

# Each square's ID, by name.
.su_ids <- function(out) stats::setNames(out$poly_id, out$name)[sort(out$name)]

test_that("a layer far from its CRS's origin keeps the IDs it has at the origin", {
  for (m in c("centroid", "surface_point", "bbox_center")) {
    for (mv in c(TRUE, FALSE)) {
      ref <- .su_ids(ensure_stable_poly_id(.su_layer(), method = m, make_valid = mv))
      # The order at the origin is the columns west to east.
      expect_identical(names(ref)[order(ref)], sprintf("p%02d", 1:12))
      for (off in c(1e8, 1e12, 1e15)) {
        lay <- .su_layer(off)
        expect_no_warning(out <- ensure_stable_poly_id(lay, method = m, make_valid = mv))
        expect_identical(.su_ids(out), ref,
                         info = sprintf("%s, make_valid = %s, offset %g", m, mv, off))
        # The same square gets the same ID whatever order the rows come in.
        rev <- ensure_stable_poly_id(lay[nrow(lay):1, ], method = m, make_valid = mv)
        expect_identical(.su_ids(rev), ref)
      }
    }
  }
})

test_that("the fallback is logged with the reason, and a placed layer takes the sort CRS", {
  far <- capture_spatialkit_log(ensure_stable_poly_id(.su_layer(1e8)))
  expect_true(log_has(far, "the sort CRS cannot place every row"))
  near <- capture_spatialkit_log(ensure_stable_poly_id(.su_layer()))
  expect_false(log_has(near, "could not transform to the sort CRS"))
})

test_that("the rows of a valid layer are the repair of every row, in ID order", {
  lay <- .su_layer()
  for (m in c("centroid", "surface_point", "bbox_center")) {
    out <- ensure_stable_poly_id(lay, method = m)
    expect_identical(sf::st_geometry(out)[order(out$name)],
                     .make_valid_rows_all(sf::st_geometry(lay))[order(lay$name)])
  }
})

test_that("create_grid_polygons_cached() numbers the cells far out as at the origin", {
  sq <- function(off) {
    p <- sf::st_polygon(list(rbind(c(599750, 3989500), c(604250, 3989500),
                                   c(604250, 3993500), c(599750, 3993500),
                                   c(599750, 3989500))))
    sf::st_sf(geometry = sf::st_sfc(p + c(off, off), crs = 32617))
  }
  centres <- function(g, off) {
    xy <- sf::st_coordinates(.su_q(sf::st_centroid(sf::st_geometry(g))))
    cbind(xy[, 1] - off, xy[, 2] - off)[order(g$poly_id), , drop = FALSE]
  }
  for (type in c("square", "hex")) {
    ref <- centres(create_grid_polygons_cached(sq(0), target_cells = 12, type = type,
                                               quiet = TRUE, cache_env = new.env()), 0)
    for (off in c(1e8, 1e12)) {
      expect_no_warning(g <- create_grid_polygons_cached(
        sq(off), target_cells = 12, type = type, quiet = TRUE, cache_env = new.env()))
      # The cells are laid at 1e12 to about 1e-4 m (the spacing of doubles
      # there), so their centres are compared to 1 cm, not to the bit.
      got <- centres(g, off)
      expect_identical(dim(got), dim(ref))
      expect_lt(max(abs(got - ref)), 0.01)
    }
  }
})

test_that("a ring with fewer than three distinct vertices is keyed by its box", {
  r <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
  g <- sf::st_sfc(sf::st_polygon(list(r(5, 5, 5, 5))),
                  sf::st_polygon(list(r(0, 0, 1, 0, 1, 1, 0, 1, 0, 0))),
                  sf::st_polygon(list(r(9, 0, 10, 0, 10, 1, 9, 1, 9, 0))),
                  crs = 32617)
  for (m in c("centroid", "surface_point")) {
    out <- ensure_stable_poly_id(g, method = m, make_valid = FALSE,
                                 transform_for_sort = NULL)
    # Keyed at (5, 5), between the squares at x = 0.5 and x = 9.5.
    first <- vapply(sf::st_geometry(out), function(s) s[[1]][1, 1], numeric(1))
    expect_identical(first, c(0, 5, 9), info = m)
  }
})

test_that("a ring PROJ leaves with two vertices no longer stops the call", {
  r <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
  g <- sf::st_sfc(sf::st_polygon(list(r(0, 0, 1e30, 0, 1e30, 1e30, 0, 0))),
                  sf::st_polygon(list(r(6e5, 4e6, 6e5 + 10, 4e6, 6e5 + 10, 4e6 + 10,
                                        6e5, 4e6 + 10, 6e5, 4e6))),
                  crs = 32617)
  for (mv in c(TRUE, FALSE)) {
    out <- .su_q(ensure_stable_poly_id(g, method = "surface_point", make_valid = mv))
    expect_identical(out$poly_id, 1:2)
    rev <- .su_q(ensure_stable_poly_id(g[2:1], method = "surface_point", make_valid = mv))
    expect_identical(sf::st_as_text(sf::st_geometry(rev)),
                     sf::st_as_text(sf::st_geometry(out)))
  }
})

test_that("an empty polygon comes back as an empty polygon of its type, with its ID", {
  lay <- .su_layer()[1:6, ]
  g <- sf::st_geometry(lay)
  g[2] <- sf::st_sfc(sf::st_polygon(), crs = 32617)
  sf::st_geometry(lay) <- g
  for (m in c("centroid", "bbox_center")) {
    expect_no_warning(out <- ensure_stable_poly_id(lay, method = m))
    types <- as.character(sf::st_geometry_type(out, by_geometry = TRUE))
    expect_identical(types, rep("POLYGON", 6L))
    # Keyed after every polygon with a location.
    expect_identical(out$name[6], lay$name[2])
    expect_true(sf::st_is_empty(sf::st_geometry(out)[[6]]))
  }
  mp <- sf::st_sf(name = c("a", "b"), geometry = sf::st_sfc(
    sf::st_multipolygon(),
    sf::st_multipolygon(list(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))),
    crs = 32617))
  out <- ensure_stable_poly_id(mp)
  expect_identical(out$name, c("b", "a"))
  expect_identical(sf::st_geometry(out)[[2]], sf::st_multipolygon())
})

# Two lon/lat rings with a vertex every 0.25 degrees, every vertex at
# latitude 87: rings with no height.
.su_flat <- function() {
  ring <- function(x0) {
    b <- sf::st_polygon(list(rbind(c(x0, 86), c(x0 + 1, 86), c(x0 + 1, 88),
                                   c(x0, 88), c(x0, 86))))
    m <- sf::st_coordinates(sf::st_segmentize(sf::st_sfc(b), 0.25))[, 1:2]
    m[, 2] <- 87
    sf::st_polygon(list(unname(m)))
  }
  sf::st_sf(name = c("q", "p"), geometry = sf::st_sfc(ring(2), ring(0), crs = 4326))
}

test_that("a ring with no area comes back as given, with a classed warning naming it", {
  flat <- .su_flat()
  for (m in c("centroid", "surface_point", "bbox_center")) {
    w <- NULL
    out <- withCallingHandlers(
      ensure_stable_poly_id(flat, method = m),
      spatialkit_cells_degenerate = function(cnd) {
        w <<- cnd
        invokeRestart("muffleWarning")
      })
    expect_s3_class(w, "spatialkit_cells_degenerate")
    expect_identical(w$rows, 1:2)
    expect_identical(w$reasons, c("has no area", "has no area"))
    expect_match(conditionMessage(w), "rows 1 and 2 of the 2 of `polygons_sf`", fixed = TRUE)
    expect_identical(as.character(sf::st_geometry_type(out, by_geometry = TRUE)),
                     c("POLYGON", "POLYGON"))
    expect_identical(sf::st_geometry(out)[order(out$name)],
                     sf::st_geometry(flat)[order(flat$name)])
    expect_identical(out$name, c("p", "q"), info = m)
  }
  # Without the repair the rows came back as given, and still do.
  out <- .su_q(ensure_stable_poly_id(flat, make_valid = FALSE))
  expect_identical(out$name, c("p", "q"))
  expect_identical(sf::st_geometry(out)[order(out$name)],
                   sf::st_geometry(flat)[order(flat$name)])
})

test_that("the warning names the rows of polygons_sf as the caller gave them", {
  flat <- .su_flat()
  pt <- sf::st_sf(name = "pt", geometry = sf::st_sfc(sf::st_point(c(5, 40)), crs = 4326))
  sq <- sf::st_sf(name = "s", geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(5, 40), c(6, 40), c(6, 41), c(5, 41), c(5, 40)))), crs = 4326))
  lay <- rbind(pt, flat[1, ], sq, flat[2, ])
  w <- NULL
  out <- withCallingHandlers(
    ensure_stable_poly_id(lay),
    spatialkit_cells_degenerate = function(cnd) { w <<- cnd; invokeRestart("muffleWarning") },
    warning = function(cnd) invokeRestart("muffleWarning"))
  expect_identical(w$rows, c(2L, 4L))
  expect_identical(nrow(out), 3L)
})

test_that("a row the repair turns into a collection keeps its polygonal part", {
  r <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
  mp <- sf::st_multipolygon(list(list(r(0, 0, 10, 0, 10, 10, 0, 10, 0, 0)),
                                 list(r(20, 0, 30, 0, 20, 0, 20, 0))))
  g <- sf::st_sfc(mp, crs = 32617)
  expect_identical(as.character(sf::st_geometry_type(.safe_make_valid(g))),
                   "GEOMETRYCOLLECTION")
  expect_no_warning(out <- ensure_stable_poly_id(g))
  expect_identical(as.character(sf::st_geometry_type(out)), "POLYGON")
  expect_equal(as.numeric(sf::st_area(out)), 100)
})

test_that("valid layers come back as the repair of every row gives them", {
  lay <- .su_layer()
  ll  <- sf::st_transform(lay, 4326)
  for (s2 in c(TRUE, FALSE)) {
    old <- sf::sf_use_s2()
    suppressMessages(sf::sf_use_s2(s2))
    for (x in list(lay, ll, sf::st_geometry(lay), sf::st_set_crs(lay, NA))) {
      expect_identical(.make_valid_rows(x), .make_valid_rows_all(x))
      rep <- .make_valid_rows(x, report = TRUE)
      expect_identical(rep$no_area, integer(0))
    }
    suppressMessages(sf::sf_use_s2(old))
  }
})

# ---------------------------------------------------------------------------
# Review of batch grid-bug-217.
#
# R1. The key of a layer keyed in its own CRS was rounded at each feature's
#     own power-of-two step, so the cells of one column whose northings
#     straddle 2^40 took two steps and were ordered by the rounding, not by
#     y.  The first key now takes one step for the whole layer.
# R2. A ring of two points stopped the default repair with GEOS's "Invalid
#     number of points in LinearRing found 2", and in degrees sf read it as
#     the whole globe and every call stopped with "To set the crs to NA,
#     first remove the full polygons".
# R3. A layer without a CRS, or with transform_for_sort = NULL, kept the
#     7-decimal key at any magnitude, so a hex grid 1e12 m out took its
#     order from the last bits of its centres.
# ---------------------------------------------------------------------------

# A 4 x 4 grid of squares of side w laid from corner (o, o); the cells of one
# column share their x exactly.
.su_grid <- function(o, w = 100.11, crs = 32617) {
  cells <- list()
  lab <- character(0)
  for (i in 0:3) for (j in 0:3) {
    x <- o + i * w
    y <- o + j * w
    cells[[length(cells) + 1L]] <- sf::st_polygon(list(rbind(
      c(x, y), c(x + w, y), c(x + w, y + w), c(x, y + w), c(x, y))))
    lab <- c(lab, sprintf("c%dr%d", i, j))
  }
  sf::st_sf(name = lab, geometry = sf::st_sfc(cells, crs = crs))
}

test_that("a column that straddles a power of two is ordered by y, as at the origin", {
  # At o = 2^40 - 250 the northings of a column run from below 2^40 to above
  # it.  Rounded at each cell's own step (1/16 below, 1/8 above), the same x
  # gave two keys and the column came out r2 r3 r0 r1.
  ref <- sprintf("c%dr%d", rep(0:3, each = 4), rep(0:3, 4))
  expect_identical(ensure_stable_poly_id(.su_grid(0), transform_for_sort = NULL)$name, ref)
  for (o in c(2^40 - 250, 2^40 - 1e6, 2^40 + 1e6)) {
    for (m in c("centroid", "bbox_center")) {
      lay <- .su_grid(o)
      expect_identical(ensure_stable_poly_id(lay, method = m)$name, ref,
                       info = sprintf("%s, o = %.0f", m, o))
      expect_identical(ensure_stable_poly_id(lay[16:1, ], method = m)$name, ref)
      expect_identical(ensure_stable_poly_id(sf::st_set_crs(lay, NA), method = m)$name, ref)
    }
  }
})

test_that("a hex grid 1e12 out without a CRS, or sorted in its own CRS, keeps its IDs", {
  p <- sf::st_polygon(list(rbind(c(599750, 3989500), c(604250, 3989500),
                                 c(604250, 3993500), c(599750, 3993500),
                                 c(599750, 3989500))))
  centres <- function(off, crs) {
    b <- sf::st_sf(geometry = sf::st_sfc(p + c(off, off), crs = crs))
    g <- .su_q(create_grid_polygons(b, target_cells = 12, type = "hex", quiet = TRUE))
    g <- ensure_stable_poly_id(g, transform_for_sort = NULL)
    xy <- sf::st_coordinates(.su_q(sf::st_centroid(sf::st_geometry(g))))
    cbind(xy[, 1] - off, xy[, 2] - off)[order(g$poly_id), , drop = FALSE]
  }
  for (crs in list(sf::NA_crs_, sf::st_crs(32617))) {
    ref <- centres(0, crs)
    got <- centres(1e12, crs)
    expect_identical(dim(got), dim(ref))
    # The cells are laid at 1e12 to about 1e-4 m, so the centres are
    # compared to 1 cm, not to the bit.
    expect_lt(max(abs(got - ref)), 0.01)
  }
})

test_that("polygons that tie at the layer's step are ordered by their own keys", {
  # One polygon 1e12 m out sets the layer's step to 1/16 m.  Two squares near
  # the origin 1 cm apart in x tie there, and are told apart by keys rounded
  # at their own step, not by the row order.
  sq <- function(x, y, s) sf::st_polygon(list(rbind(c(x, y), c(x + s, y),
                                                    c(x + s, y + s), c(x, y + s), c(x, y))))
  for (crs in list(sf::st_crs(32617), sf::NA_crs_)) {
    lay <- sf::st_sf(name = c("far", "b", "a"), geometry = sf::st_sfc(
      sq(-1e12, 4e6, 1000), sq(600000.01, 4e6, 1), sq(600000, 4e6, 1), crs = crs))
    for (rows in list(1:3, 3:1, c(2, 1, 3))) {
      out <- ensure_stable_poly_id(lay[rows, ])
      expect_identical(out$name, c("far", "a", "b"))
    }
  }
})

test_that("a key in the layer's unit is rounded to 7 decimals near the origin", {
  # Centroids 4e-8 apart in x round to one key, so y orders them.
  sq <- function(x, y) sf::st_polygon(list(rbind(c(x, y), c(x + 1, y),
                                                 c(x + 1, y + 1), c(x, y + 1), c(x, y))))
  lay <- sf::st_sf(name = c("west_high", "east_low"), geometry = sf::st_sfc(
    sq(-4e-8, 0), sq(0, -5)))
  expect_identical(ensure_stable_poly_id(lay)$name, c("east_low", "west_high"))
  expect_identical(ensure_stable_poly_id(lay[2:1, ])$name, c("east_low", "west_high"))
})

test_that("a transform_for_sort PROJ cannot read gives the IDs of NULL, without GDAL's warning", {
  lay <- .su_layer()
  ref <- ensure_stable_poly_id(lay, transform_for_sort = NULL)
  expect_no_warning(log <- capture_spatialkit_log(
    out <- ensure_stable_poly_id(lay, transform_for_sort = "+proj=bogus")))
  expect_true(log_has(log, "could not transform to the sort CRS"))
  expect_identical(out$name, ref$name)
})

test_that("a ring of two points is keyed by its box, with or without the repair, in any CRS", {
  r <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
  for (crs in c(32617, 4326)) {
    g <- sf::st_sf(name = c("ring2", "west", "east"), geometry = sf::st_sfc(
      sf::st_polygon(list(r(5, 5, 5, 5))),
      sf::st_polygon(list(r(0, 0, 1, 0, 1, 1, 0, 1, 0, 0))),
      sf::st_polygon(list(r(9, 0, 10, 0, 10, 1, 9, 1, 9, 0))),
      crs = crs))
    for (tf in list(4326, NULL)) for (m in c("centroid", "surface_point", "bbox_center")) {
      info <- sprintf("crs %d, %s, transform_for_sort %s", crs, m, format(tf))
      # With the repair: the ring is returned as given, named by the warning.
      w <- NULL
      out <- withCallingHandlers(
        ensure_stable_poly_id(g, method = m, transform_for_sort = tf),
        spatialkit_cells_degenerate = function(cnd) {
          w <<- cnd
          invokeRestart("muffleWarning")
        })
      expect_identical(w$rows, 1L, info = info)
      # Keyed at (5, 5), between the squares at x = 0.5 and x = 9.5.
      expect_identical(out$name, c("west", "ring2", "east"), info = info)
      expect_identical(sf::st_geometry(out)[[2]], sf::st_geometry(g)[[1]])
      # Without the repair: the same order, no warning.
      expect_no_warning(out <- ensure_stable_poly_id(g, method = m, make_valid = FALSE,
                                                     transform_for_sort = tf))
      expect_identical(out$name, c("west", "ring2", "east"), info = info)
      expect_identical(ensure_stable_poly_id(g[3:1, ], method = m, make_valid = FALSE,
                                             transform_for_sort = tf)$name,
                       c("west", "ring2", "east"), info = info)
    }
  }
})

test_that("a short ring in a polygon with area is dropped by the repair, not the polygon", {
  r <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)
  shell <- r(0, 0, 10, 0, 10, 10, 0, 10, 0, 0)
  g <- sf::st_sfc(
    sf::st_polygon(list(shell, r(5, 5, 5, 5))),
    sf::st_multipolygon(list(list(shell + 20), list(r(40, 0, 40, 0)))),
    crs = 32617)
  expect_no_warning(out <- ensure_stable_poly_id(g))
  expect_true(all(sf::st_geometry_type(out, by_geometry = TRUE) %in%
                    c("POLYGON", "MULTIPOLYGON")))
  expect_equal(as.numeric(sf::st_area(out)), c(100, 100))
  expect_true(all(sf::st_is_valid(out)))
  # The polygon keeps its shell and loses the ring of two points.
  expect_identical(length(sf::st_geometry(out)[[1]]), 1L)
})
