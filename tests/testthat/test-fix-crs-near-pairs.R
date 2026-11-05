# tests/testthat/test-fix-crs-near-pairs.R
# ---------------------------------------------------------------------------
# .crs_distance_error() takes the worst |d_planar / d_geodesic - 1| over the
# sampled pairs.  Two listings of one place that are not the same double (a
# corner two features share, computed two ways; a point listed twice a few
# bits apart; a point at longitude 180 and its copy at -180) made a pair
# about 1e-9 m long, whose ratio is rounding noise or, across the
# antimeridian, the width of the world, and that pair set the score of every
# candidate.  Pairs under 1 mm are now left out.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
box <- function(x0, y0, w, h)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
layer <- function(...)
  sf::st_sf(id = seq_along(list(...)), geometry = sf::st_sfc(list(...), crs = 4326))
points_of <- function(m)
  sf::st_sf(id = seq_len(nrow(m)), geometry = sf::st_sfc(lapply(seq_len(nrow(m)), function(i)
    sf::st_point(m[i, ])), crs = 4326))
# The two figures, in percent, of the line .pick_local_projected_crs() logs.
quoted <- function(lines) {
  line <- grep(".pick_local_projected_crs()", lines, fixed = TRUE, value = TRUE)
  stopifnot(length(line) == 1L)
  as.numeric(sub("%", "", regmatches(line, gregexpr("[0-9]+\\.[0-9]+%", line))[[1L]], fixed = TRUE))
}
# -110.6 + 3.4 is 1.4e-14 off -107.2.
side_lit <- layer(box(-110.6, 37.9, 3.4, 9.1), box(-107.2, 37.9, 3.4, 9.1))
side_one <- layer(box(-110.6, 37.9, 6.8, 9.1))

test_that("corners typed two ways score as the outline and as the corner computed", {
  # The two boxes scored 0.137 in UTM zone 13 and 0.245 in EPSG:5070, where
  # their outline scores 0.0026 and 0.0097; the T with its corner typed
  # -103.6 scored 0.137 in both, where the corner computed gives 0.0058 and
  # 0.0098.
  tee <- function(corner) layer(box(-109.7, 31.3, 6.1, 14.6), box(corner, 31.3, 6.1, 7.3))
  for (code in c(32613, 5070)) {
    crs <- sf::st_crs(code)
    expect_lt(abs(.crs_distance_error(side_lit, crs) - .crs_distance_error(side_one, crs)), 1e-3,
              label = paste("two boxes in EPSG", code))
    expect_lt(abs(.crs_distance_error(tee(-103.6), crs) - .crs_distance_error(tee(-109.7 + 6.1), crs)),
              1e-3, label = paste("T in EPSG", code))
  }
})

test_that("ensure_projected() quotes the outline's figures for corners typed two ways", {
  # The two boxes got an Albers CRS with "10.17% against the zone's 13.70%";
  # their outline gets a Lambert azimuthal CRS with "0.09% against the zone's
  # 0.26%".
  x <- NULL; x1 <- NULL
  l  <- capture_spatialkit_log(x  <- quiet(ensure_projected(side_lit)))
  l1 <- capture_spatialkit_log(x1 <- quiet(ensure_projected(side_one)))
  ch <- attr(x, "crs_choice"); ch1 <- attr(x1, "crs_choice")
  expect_identical(ch$name[ch$chosen], ch1$name[ch1$chosen])
  expect_match(ch$name[ch$chosen], "^Lambert azimuthal")
  q <- quoted(l); q1 <- quoted(l1)
  expect_length(q, 2L)
  expect_lt(max(abs(q - q1)), 0.1 + 1e-9)
})

test_that("a thinned border written two ways scores as the border written one way", {
  # 206 distinct vertices, so the outline is thinned.  0.125 and 0.189,
  # against 0.0011 and 0.0098.
  yy <- seq(37.9, 47, length.out = 101)
  west <- rbind(c(-110.6, 37.9), cbind(-107.2, yy), c(-110.6, 47), c(-110.6, 37.9))
  east <- function(x) rbind(c(x, 37.9), c(-103.8, 37.9), c(-103.8, 47), cbind(x, rev(yy)))
  two <- layer(sf::st_polygon(list(west)), sf::st_polygon(list(east(-110.6 + 3.4))))
  one <- layer(sf::st_polygon(list(west)), sf::st_polygon(list(east(-107.2))))
  v <- sf::st_coordinates(two)[, 1:2]
  expect_gt(sum(!duplicated(v)), 160)
  for (code in c(32613, 5070))
    expect_lt(abs(.crs_distance_error(two, sf::st_crs(code)) - .crs_distance_error(one, sf::st_crs(code))),
              1e-4, label = paste("EPSG", code))
})

test_that("a point listed twice a few bits apart scores as the point listed once", {
  # 0.089 in UTM zone 13 and 0.180 in EPSG:5070, where the 30 points score
  # 0.0004 and 0.0097.
  set.seed(7)
  pts <- cbind(stats::runif(30, -107.5, -102.5), stats::runif(30, 33, 41))
  for (code in c(32613, 5070)) {
    crs <- sf::st_crs(code)
    e0 <- .crs_distance_error(points_of(pts), crs)
    expect_identical(.crs_distance_error(points_of(rbind(pts, pts[5, ] + c(1e-14, 0))), crs), e0,
                     label = paste("EPSG", code))
    expect_identical(.crs_distance_error(points_of(rbind(pts, pts[5, ])), crs), e0)
  }
})

test_that("a layer split at the antimeridian scores as the layer listed one way", {
  # The split boxes reported 1 for the one candidate; the ring round the
  # South Pole quoted 1.5e19 percent for Web Mercator.
  split2 <- layer(box(170, -20, 10, 8), box(-180, -20, 10, 8))
  split1 <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_multipolygon(list(
    box(170, -20, 10, 8), box(-180, -20, 10, 8))), crs = 4326))
  wrote <- layer(sf::st_polygon(list(rbind(c(170, -20), c(-180, -20), c(-180, -12), c(170, -12),
                                           c(170, -20)))), box(-180, -20, 10, 8))
  e_w <- attr(quiet(ensure_projected(wrote)), "crs_choice")$distance_error
  expect_length(e_w, 1L)
  for (lay in list(split2, split1)) {
    ch <- attr(quiet(ensure_projected(lay)), "crs_choice")
    expect_equal(nrow(ch), 1L)
    expect_lt(abs(ch$distance_error - e_w), 1e-4)
  }
  ring <- layer(box(-180, -80, 120, 15), box(-60, -80, 120, 15), box(60, -80, 120, 15))
  ch <- attr(quiet(ensure_projected(ring)), "crs_choice")
  expect_match(ch$name[ch$chosen], "South Pole")
  expect_lt(ch$distance_error[ch$name == "Web Mercator (EPSG:3857)"], 1e3)
})

test_that("pairs under 1 mm are left out and pairs from 1 mm are measured", {
  # North-south, where .geod_distance() is exact at these lengths.
  near <- points_of(rbind(c(-104, 40), c(-104, 40 + 4e-9)))   # about 0.44 mm
  far  <- points_of(rbind(c(-104, 40), c(-104, 40 + 2e-8)))   # about 2.2 mm
  just <- points_of(rbind(c(-104, 40), c(-104, 40 + 1e-8)))   # about 1.1 mm
  expect_true(is.na(.crs_distance_error(near, sf::st_crs(32613))))
  expect_true(is.finite(.crs_distance_error(far, sf::st_crs(32613))))
  expect_gt(.geod_distance(-104, 40, -104, 40 + 1e-8), 1e-3)
  expect_lt(.geod_distance(-104, 40, -104, 40 + 1e-8), 1.15e-3)
  expect_true(is.finite(.crs_distance_error(just, sf::st_crs(32613))))
})

test_that("layers with no pair under 1 mm score as before", {
  # The figures recorded before the change, to the last bits the change of
  # the geodesic's stopping test moves (test-fix-geod-short-pairs.R).
  # BUG-23's pins in test-fix-crs-shared-edges.R cover the polygon layers.
  set.seed(7)
  pts <- points_of(cbind(stats::runif(30, -107.5, -102.5), stats::runif(30, 33, 41)))
  aea <- sf::st_crs(paste("+proj=aea +lat_1=35 +lat_2=55 +lat_0=40 +lon_0=-104",
                          "+datum=WGS84 +units=m +no_defs"))
  expect_equal(.crs_distance_error(pts, sf::st_crs(32613)), 0.0003982212900482951, tolerance = 1e-8)
  expect_equal(.crs_distance_error(pts, aea), 0.010731689399051114, tolerance = 1e-8)
})
