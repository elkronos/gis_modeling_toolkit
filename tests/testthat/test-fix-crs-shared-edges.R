# tests/testthat/test-fix-crs-shared-edges.R
# ---------------------------------------------------------------------------
# A lon/lat layer of fewer than 40 features is scored on its outlines'
# vertices too, and when they are few, on points added along each edge.  An
# edge two features share is listed once by each, from opposite ends, so its
# added points were computed twice, once from each end, and the two copies
# of a point could differ in their last bits.  unique() kept both, and the
# pair, about 1e-9 m apart, was measured on rounding noise: every candidate
# CRS scored on it, the CRS was chosen on it, and the warning quoted it.  The
# same happens where a point added along one feature's edge falls on a
# vertex of another (a T-junction, or a corner touching the edge), or on a
# point added along an edge that overlaps it in part.  Each edge is now
# filled from the same end whichever feature lists it, and an added point
# within about 1e-9 degrees of a point listed before it is left out.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
box <- function(x0, y0, w, h)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
boxes <- function(nx, ny, x0, y0, w, h)
  sf::st_sf(id = seq_len(nx * ny), geometry = sf::st_sfc(unlist(lapply(0:(ny - 1), function(j)
    lapply(0:(nx - 1), function(i) box(x0 + i * w, y0 + j * h, w, h))), recursive = FALSE),
    crs = 4326))
layer <- function(...)
  sf::st_sf(id = seq_along(list(...)), geometry = sf::st_sfc(list(...), crs = 4326))
# Boxes from x0 to x1, one above the other, split at the latitudes ys.
stack <- function(x0, x1, ys) do.call(layer, lapply(seq_len(length(ys) - 1L), function(j)
  sf::st_polygon(list(rbind(c(x0, ys[j]), c(x1, ys[j]), c(x1, ys[j + 1]), c(x0, ys[j + 1]),
                            c(x0, ys[j]))))))
# A box of w by 2h beside one of w by h: the corner of the second lies on
# the edge of the first, half way along it.
tee <- function(x0, y0, w, h) layer(box(x0, y0, w, 2 * h), box(x0 + w, y0, w, h))
# A box and a triangle whose apex lies on the middle of the box's top edge.
# The two share no edge.
apex <- function(x0, y0, w, h) {
  ap <- c(x0 + w / 2, y0 + h)
  layer(box(x0, y0, w, h),
        sf::st_polygon(list(rbind(ap, c(x0 + w, y0 + h + 4), c(x0, y0 + h + 4), ap))))
}
# The layer with every vertex rounded to 9 decimal places, where a corner on
# another feature's edge and the point added there are the same double.
snapped <- function(lay) {
  g <- lapply(sf::st_geometry(lay), function(s) {
    s[] <- rapply(unclass(s), function(m) round(m, 9), how = "replace"); s })
  sf::st_geometry(lay) <- sf::st_sfc(g, crs = sf::st_crs(lay))
  lay
}
# The layer's outline as one box, with no edge inside it.
outline_of <- function(lay) sf::st_sf(id = 1, geometry = sf::st_as_sfc(sf::st_bbox(lay)))
shared <- list(
  "two boxes of 12 by 15 degrees, one above the other" = boxes(1, 2, -110, 30, 12, 15),
  "a 2 by 2 grid of boxes of 3.7 by 2.9 degrees"       = boxes(2, 2, -104.3, 33.1, 3.7, 2.9),
  "two boxes of 3.4 by 9.1 degrees, side by side"      = boxes(2, 1, -110.6, 37.9, 3.4, 9.1),
  # Corners given to 8 decimals: with three points added per edge, a point
  # added along a shared edge lies half way between two steps of the 9th
  # decimal, and its two copies were rounded to either side.
  "four boxes one above the other, corners to 8 decimals" = stack(
    -104.02908244, -101.05384711,
    c(33.23366968, 33.89516779, 34.64330048, 35.27466651, 35.87610188)))
# The two figures, in percent, of the line .pick_local_projected_crs() logs,
# or nothing when it logs no line.
quoted <- function(lines) {
  line <- grep(".pick_local_projected_crs()", lines, fixed = TRUE, value = TRUE)
  stopifnot(length(line) <= 1L)
  if (!length(line)) return(numeric(0))
  as.numeric(sub("%", "", regmatches(line, gregexpr("[0-9]+\\.[0-9]+%", line))[[1L]], fixed = TRUE))
}

test_that("features that share edges score as their outline as one box does", {
  # The two boxes scored 0.769 in UTM zone 13 and 0.868 in EPSG:5070, where
  # their outline scores 0.005 and 0.109; the grid 0.272 and 0.096, against
  # 0.007 and 0.010; the boxes side by side 0.318 and 1, against 0.003 and
  # 0.010.  The four boxes scored 0.112 in both; rounded to 9 decimals they
  # still scored 0.0029 in UTM zone 13, against 0.0013 for their outline.
  for (nm in names(shared)) {
    lay <- shared[[nm]]
    expect_lt(nrow(unique(sf::st_coordinates(lay)[, 1:2])), 40)   # densified
    for (code in c(32613, 5070)) {
      e  <- .crs_distance_error(lay, sf::st_crs(code))
      e1 <- .crs_distance_error(outline_of(lay), sf::st_crs(code))
      expect_true(is.finite(e) && is.finite(e1))
      expect_lt(abs(e - e1), 1e-3, label = paste(nm, "in EPSG", code))
    }
  }
})

test_that("an edge shared across the antimeridian scores as the outline does", {
  # Two boxes from 175.3 to -175.3, one above the other, scored 0.792 in a
  # Lambert azimuthal CRS on 180, where their outline scores 0.0025.
  am <- function(y0, h) sf::st_polygon(list(rbind(c(175.3, y0), c(-175.3, y0), c(-175.3, y0 + h),
                                                  c(175.3, y0 + h), c(175.3, y0))))
  laea <- sf::st_crs("+proj=laea +lat_0=-12 +lon_0=180 +datum=WGS84 +units=m +no_defs")
  e  <- .crs_distance_error(layer(am(-20, 8), am(-12, 8)), laea)
  e1 <- .crs_distance_error(layer(am(-20, 16)), laea)
  expect_true(is.finite(e) && is.finite(e1))
  expect_lt(abs(e - e1), 1e-3)
})

test_that("a ring listed the other way round scores the same", {
  # Turned clockwise, the lower box lists the shared edge in the direction
  # the upper box does, the points added along it are computed from the
  # same end, and there was no noise to measure: the layer scored 0.005 in
  # UTM zone 13 with one ring turned and 0.769 with neither.
  lay <- shared[[1L]]
  rev_ring <- lay
  sf::st_geometry(rev_ring)[1] <- sf::st_sfc(sf::st_polygon(list(
    sf::st_geometry(lay)[[1]][[1]][5:1, ])), crs = 4326)
  for (code in c(32613, 5070, 3857))
    expect_equal(.crs_distance_error(rev_ring, sf::st_crs(code)),
                 .crs_distance_error(lay, sf::st_crs(code)), tolerance = 1e-12,
                 label = paste("EPSG", code))
})

test_that("the score does not depend on the order of the features or where a ring starts", {
  # Each edge is filled from the same end, and of the points that fall on
  # one another the first in coordinate order is kept, so the points
  # measured are the same doubles whatever the order of the listing.
  turned <- function(lay) {
    g <- lapply(sf::st_geometry(lay), function(s) sf::st_polygon(list(s[[1]][nrow(s[[1]]):1, ])))
    sf::st_geometry(lay) <- sf::st_sfc(g, crs = 4326)
    lay
  }
  started <- function(lay) {
    g <- lapply(sf::st_geometry(lay), function(s) {
      m <- s[[1]][-nrow(s[[1]]), , drop = FALSE]
      m <- rbind(m[-1L, , drop = FALSE], m[1L, ])
      sf::st_polygon(list(rbind(m, m[1L, ])))
    })
    sf::st_geometry(lay) <- sf::st_sfc(g, crs = 4326)
    lay
  }
  # Without the sort the two overlapping boxes, listed the other way round,
  # scored 0.0097751547341586864 against 0.0097751547341604628 in EPSG:5070.
  lays <- c(shared, list(tee = tee(-109.7, 31.3, 6.1, 7.3), apex = apex(-105.7, 31.4, 6, 6.8),
                         overlap = layer(box(-100.9, 35.1, 3.4, 4.6),
                                         box(-100.9, 35.1, 2 * 3.4, 4.6 / 2))))
  for (nm in names(lays)) {
    lay <- lays[[nm]]
    for (code in c(32613, 5070)) {
      crs <- sf::st_crs(code)
      e <- .crs_distance_error(lay, crs)
      expect_identical(.crs_distance_error(lay[rev(seq_len(nrow(lay))), ], crs), e,
                       label = paste(nm, "with its rows reversed, EPSG", code))
      expect_identical(.crs_distance_error(turned(lay), crs), e,
                       label = paste(nm, "with its rings turned, EPSG", code))
      expect_identical(.crs_distance_error(started(lay), crs), e,
                       label = paste(nm, "with its rings started at the next vertex, EPSG", code))
    }
  }
})

test_that("ensure_projected() quotes the outline's figures for features that share edges", {
  # The two boxes kept UTM zone 13 with "76.92% vs 76.92% worst-case over
  # sampled pairs", where their outline gives 1.04% and 0.51%; the grid moved
  # to a Lambert azimuthal CRS with "11.22% against the zone's 13.50%", where
  # its outline gives 0.06% and 0.26%.  The four boxes keep UTM zone 13 and
  # log no line; with their copies rounded to 9 decimals, attr(, "crs_choice")
  # gave the zone 0.29% where the outline gives 0.13%.
  for (nm in names(shared)) {
    lay <- shared[[nm]]
    one <- outline_of(lay)
    x  <- NULL; x1 <- NULL
    l  <- capture_spatialkit_log(x  <- quiet(ensure_projected(lay)))
    l1 <- capture_spatialkit_log(x1 <- quiet(ensure_projected(one)))
    q <- quoted(l); q1 <- quoted(l1)
    expect_length(q, length(q1))
    if (length(q1)) {
      expect_length(q, 2L)
      expect_lt(max(abs(q - q1)), 0.1, label = nm)
    }
    ch <- attr(x, "crs_choice"); ch1 <- attr(x1, "crs_choice")
    expect_identical(ch$name[ch$chosen], ch1$name[ch1$chosen], label = nm)
    expect_lt(max(abs(ch$distance_error - ch1$distance_error)), 1e-3, label = nm)
  }
})

test_that("a corner on another feature's edge scores as the features' union and the snapped layer", {
  # The point added half way along the long box's edge of a T falls on the
  # corner of the short box, the point added on the middle of the box's top
  # edge on the triangle's apex, and points added along two edges that
  # overlap in part fall on one another.  The corner is a computed double
  # (-109.7 + 6.1 is -103.60000000000001).  Rounding the added point alone
  # kept the two 1e-14 degrees apart: the first T scored 0.058 in UTM zone 13
  # and 0.152 in EPSG:5070, the apex 0.318 and 1, where the snapped layers
  # score 0.0058, 0.0098, 0.0026 and 0.0098.  The second and third T scored
  # 0.561 and 1, and 0.318 and 1, both before the change and with the added
  # points rounded, where their unions score 0.0049 and 0.0069, and 0.0012
  # and 0.0098.  The overlapping boxes scored 0.547 and 0.538, where their
  # union scores 0.012 and 0.010.
  union_of <- function(lay) sf::st_sf(id = 1, geometry = sf::st_set_crs(
    suppressMessages(sf::st_union(sf::st_set_crs(sf::st_geometry(lay), NA))), 4326))
  lays <- list(tee = tee(-109.7, 31.3, 6.1, 7.3), "second tee" = tee(-111.7, 28.3, 2.7, 2.6),
               "third tee" = tee(-106.7, 33.9, 2.8, 4.3),
               overlap = layer(box(-100.9, 35.1, 3.4, 4.6), box(-100.9, 35.1, 2 * 3.4, 4.6 / 2)),
               apex = apex(-105.7, 31.4, 6, 6.8))
  for (nm in names(lays)) {
    lay <- lays[[nm]]
    for (code in c(32613, 5070)) {
      e  <- .crs_distance_error(lay, sf::st_crs(code))
      e1 <- .crs_distance_error(snapped(lay), sf::st_crs(code))
      expect_true(is.finite(e) && is.finite(e1))
      expect_lt(abs(e - e1), 1e-6, label = paste(nm, "against the snapped layer in EPSG", code))
      # The apex touches the box at one point, so its union is no one polygon.
      if (nm != "apex")
        expect_lt(abs(e - .crs_distance_error(union_of(lay), sf::st_crs(code))), 1e-3,
                  label = paste(nm, "against the union in EPSG", code))
    }
  }
  # ensure_projected() chose UTM zone 13 for the first T with "9.30% vs
  # 5.77%", and a Lambert azimuthal CRS for the apex with "27.68% against the
  # zone's 31.79%"; the snapped layers get a Lambert azimuthal CRS with
  # "0.30% against the zone's 0.58%" and "0.15% against the zone's 0.26%".
  for (nm in c("tee", "apex")) {
    lay <- lays[[nm]]
    x <- NULL; x1 <- NULL
    l  <- capture_spatialkit_log(x  <- quiet(ensure_projected(lay)))
    l1 <- capture_spatialkit_log(x1 <- quiet(ensure_projected(snapped(lay))))
    ch <- attr(x, "crs_choice"); ch1 <- attr(x1, "crs_choice")
    expect_match(ch$name[ch$chosen], "^Lambert azimuthal", label = nm)
    expect_identical(ch$name[ch$chosen], ch1$name[ch1$chosen], label = nm)
    q <- quoted(l); q1 <- quoted(l1)
    expect_length(q, 2L)
    expect_equal(q, q1, label = nm)
  }
})

test_that("two points a few bits apart either side of a half step of the 9th decimal are one", {
  # Two T layers with every corner 5e-10 degrees off a 9-decimal value: the
  # corner of the short box lies half way between two steps of the 9th
  # decimal, and the point added on it a few bits to one side.  The first
  # scored 0.557 in UTM zone 13 and 1 in EPSG:5070, and as much when the
  # added points are compared with the vertices on one key of 9 decimals;
  # the second 0.122 and 0.322, and as much on one grid of 1e-9 degrees.
  # Moved back by 5e-10 degrees, they score 0.0058 and 0.0098, and 0.0090
  # and 0.0191.
  pairs <- list(
    list(tee(-109.7000000005, 31.3000000005, 6.1, 7.3), tee(-109.7, 31.3, 6.1, 7.3)),
    list(tee(-108.863847813 + 5e-10, 35.738197212 + 5e-10, 6.8, 7.3),
         tee(-108.863847813, 35.738197212, 6.8, 7.3)))
  for (k in seq_along(pairs))
    for (code in c(32613, 5070))
      expect_lt(abs(.crs_distance_error(pairs[[k]][[1L]], sf::st_crs(code)) -
                    .crs_distance_error(pairs[[k]][[2L]], sf::st_crs(code))), 1e-6,
                label = paste("layer", k, "in EPSG", code))
})

test_that("layers with no added point on another feature's point score as before", {
  # Two boxes apart and two that touch at a corner, where no point added
  # along an edge can fall on another point, and the T and the apex, where
  # the added point and the corner it falls on were the same double before
  # the change.  The figures are those recorded before the change, with the
  # geodesic's stopping test measured against the pair's longitude
  # difference (test-fix-geod-short-pairs.R), which moved the Albers
  # figures of `apart` and `corner` by 2.6e-9.
  aea <- sf::st_crs(paste("+proj=aea +lat_1=35 +lat_2=55 +lat_0=40 +lon_0=-104",
                          "+datum=WGS84 +units=m +no_defs"))
  apart  <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(box(-110, 30, 5, 5), box(-100, 40, 4, 3),
                                                      crs = 4326))
  corner <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(box(-110, 30, 5, 5), box(-105, 35, 5, 5),
                                                      crs = 4326))
  expect_equal(.crs_distance_error(apart, sf::st_crs(32613)), 0.0068189865459915566,
               tolerance = 1e-9)
  expect_equal(.crs_distance_error(apart, aea), 0.015469862162994286, tolerance = 1e-9)
  expect_equal(.crs_distance_error(corner, sf::st_crs(32613)), 0.0024495482006103853,
               tolerance = 1e-9)
  expect_equal(.crs_distance_error(corner, aea), 0.015470787824463095, tolerance = 1e-9)
  t <- tee(-109.7, 31.3, 6.1, 7.3)
  a <- apex(-105.7, 31.4, 6, 6.8)
  expect_equal(.crs_distance_error(t, sf::st_crs(32613)), 0.0058229730051715389, tolerance = 1e-9)
  expect_equal(.crs_distance_error(t, aea), 0.015141462752655821, tolerance = 1e-9)
  expect_equal(.crs_distance_error(a, sf::st_crs(32613)), 0.0026202396204950684, tolerance = 1e-9)
  expect_equal(.crs_distance_error(a, aea), 0.013235802505612737, tolerance = 1e-9)
})
