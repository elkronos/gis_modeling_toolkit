# tests/testthat/test-fix-crs-area-copies.R
# ---------------------------------------------------------------------------
# The area check (.crs_area_error()) takes the spread of planar-to-geodesic
# area ratios over a layer's polygons.  It counted rows: copies of a polygon
# weighted the median by their number, took places in the sample of 200
# polygons, and two copies of one polygon passed as two probes.  So
# ensure_projected(purpose = "area") warned, and summarize_by_cell(area =
# TRUE) refused, by the number of copies.  The probes are now the distinct
# polygons, the first copy of each kept; a layer whose polygons are all
# copies of one is measured as one polygon is, on the grid over its box.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
square <- function(x0, y0, w, h = w)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
# Three squares of 50 km in UTM zone 14, far apart east to west.
utm3 <- list(square(3e5, 4.4e6, 5e4), square(7.5e5, 4.4e6, 5e4), square(1.15e6, 4.4e6, 5e4))
lay3 <- function(idx) sf::st_sf(id = seq_along(idx), geometry = sf::st_sfc(utm3[idx], crs = 32614))
# One box from 100W to 90W and 30N to 60N in Web Mercator.
merc_box <- function() sf::st_transform(sf::st_sfc(square(-100, 30, 10, 30), crs = 4326), 3857)

test_that("copies of a polygon do not move the figure, the warning or the refusal", {
  once <- .crs_area_error(lay3(1:3))
  expect_true(is.finite(once))
  expect_lt(once, .area_error_tol)
  expect_equal(once, 0.009331, tolerance = 1e-3)
  for (idx in list(c(1, 1, 1, 2, 3), c(1, 2, 3, 3, 3), c(1, 2, 3, 1, 2, 3), c(3, 3, 2, 1, 1),
                   c(2, 1, 3, 2, 2, 2, 2)))
    expect_identical(.crs_area_error(lay3(idx)), once, info = paste(idx, collapse = ","))
  # ensure_projected(purpose = "area"): no distortion warning, once or with
  # the left square listed three times (it warned with 1.0% then).
  for (idx in list(1:3, c(1, 1, 1, 2, 3))) {
    lines <- capture_spatialkit_log(out <- ensure_projected(lay3(idx), purpose = "area"))
    expect_false(log_has(lines, "distorts areas"), info = paste(idx, collapse = ","))
    expect_identical(out, lay3(idx))
  }
  # summarize_by_cell(area = TRUE): the left square listed three times as
  # three cells runs, with the figure of the three squares (it was refused).
  pts <- sf::st_sf(v = 1:3, poly_id = 1:3,
                   geometry = sf::st_centroid(sf::st_sfc(utm3, crs = 32614)))
  for (idx in list(1:3, c(1, 1, 1, 2, 3))) {
    cs <- sf::st_sf(poly_id = seq_along(idx), geometry = sf::st_sfc(utm3[idx], crs = 32614))
    s <- quiet(summarize_by_cell(pts, cells_sf = cs, area = TRUE, response_var = "v",
                                estimate = "mean"))
    expect_identical(attr(s, "area_error"), once, info = paste(idx, collapse = ","))
    expect_equal(s$cell_area, rep(2.5e9, length(idx)), tolerance = 1e-9)
  }
})

test_that("a polygon listed several times is measured as the polygon alone", {
  b <- merc_box()
  g <- .crs_area_error(sf::st_sf(a = 1, geometry = b), grid = TRUE)
  expect_gt(g, 0.6)
  for (k in c(1, 2, 5)) {
    x <- sf::st_sf(a = seq_len(k), geometry = rep(b, k))
    # Unmeasured on its polygons, as one polygon is (two copies gave 0).
    expect_true(is.na(.crs_area_error(x)), info = k)
    expect_identical(.crs_area_error(x, one_by_grid = TRUE), g, info = k)
    lines <- capture_spatialkit_log(ensure_projected(x, purpose = "area"))
    expect_true(log_has(lines, sprintf("distorts areas across this extent by up to %.1f%%",
                                       100 * g)), info = k)
    cs <- sf::st_sf(cell_id = seq_len(k), geometry = rep(b, k))
    pts <- sf::st_sf(cell_id = seq_len(k), z = seq_len(k),
                     geometry = rep(sf::st_centroid(b), k))
    expect_error(quiet(summarize_by_cell(pts, "z", id_col = "cell_id", cells_sf = cs,
                                         estimate = "mean", area = TRUE)),
                 sprintf("distorts areas across the cells by up to %.1f%%", 100 * g),
                 fixed = TRUE, info = k)
  }
})

test_that("above 200 polygons the sample is taken over the distinct polygons", {
  set.seed(4)
  sqs <- lapply(1:150, function(i) square(runif(1, -130, -71), runif(1, 20, 49), 1))
  g150 <- sf::st_transform(sf::st_sfc(sqs, crs = 4326), 32615)
  f0 <- .crs_area_error(sf::st_sf(i = 1:150, geometry = g150))
  expect_equal(f0, 0.2498, tolerance = 1e-3)
  # Square 1 listed 300 more times (450 rows, so the sample of 200 ran on
  # rows) measured 0.3079; every square three times measured 0.2498.
  expect_identical(.crs_area_error(sf::st_sf(i = 1:450, geometry = c(g150, rep(g150[1], 300)))), f0)
  expect_identical(.crs_area_error(sf::st_sf(i = 1:450, geometry = rep(g150, each = 3))), f0)
  expect_identical(.crs_area_error(sf::st_sf(i = 1:450, geometry = rep(g150, 3))), f0)
  # Copies interleaved in reverse order.
  expect_identical(.crs_area_error(sf::st_sf(i = 1:300, geometry = c(rev(g150), g150))), f0)
})

test_that("polygons that tie on the sort keys without being copies are both measured", {
  # A polygon whose vertices average to the centre of its box, and the same
  # polygon turned half a circle about that centre, agree on all nine sort
  # keys.  They are not copies, so both are probes: dropped as a copy, the
  # turned one would leave one probe and the figure NA.
  p <- sf::st_polygon(list(rbind(c(0, 0), c(4, 0), c(4, 1), c(1, 1), c(1, 3), c(3, 3),
                                 c(3, 4), c(0, 4), c(0, 0))))
  m <- unclass(p)[[1]]
  c0 <- c(2, 2)
  q <- sf::st_polygon(list(cbind(2 * c0[1] - m[, 1], 2 * c0[2] - m[, 2])))
  at <- function(poly, s = 1e5) sf::st_polygon(list(unclass(poly)[[1]] * s + cbind(rep(3e5, nrow(m)), rep(4.4e6, nrow(m)))))
  P <- at(p); Q <- at(q)
  g <- sf::st_sfc(list(P, Q), crs = 32614)
  k <- .polygon_sort_keys(g)
  expect_identical(k[1, ], k[2, ])
  expect_false(identical(P, Q))
  two <- .crs_area_error(sf::st_sf(i = 1:2, geometry = g))
  expect_true(is.finite(two))
  # Copies of P beside Q, in any order, measure as P and Q once.
  for (ord in list(c(1, 2, 1), c(2, 1, 1), c(1, 1, 2, 2, 1)))
    expect_identical(.crs_area_error(sf::st_sf(i = seq_along(ord), geometry = g[ord])), two,
                     info = paste(ord, collapse = ","))
})

test_that("layers without a polygon listed twice get the figure they got", {
  # The figure on distinct polygons is the figure over every row: the same
  # three squares measured one at a time, by hand, as the check does.
  x <- lay3(1:3)
  by_hand <- function(x) {
    g <- sf::st_geometry(x)
    pb <- sf::st_bbox(g)
    ext <- max(pb["xmax"] - pb["xmin"], pb["ymax"] - pb["ymin"])
    gd <- sf::st_segmentize(g, dfMaxLength = as.numeric(ext) / 100)
    r <- as.numeric(sf::st_area(gd)) /
      spatialkit:::.with_s2(as.numeric(sf::st_area(sf::st_transform(gd, 4326))))
    max(abs(r / stats::median(r) - 1))
  }
  expect_identical(.crs_area_error(x), by_hand(x))
  four <- sf::st_sf(id = 1:4, geometry = sf::st_sfc(lapply(0:3, function(k)
    square(4e5 + 1e4 * k, 4.4e6, 1e4)), crs = 32615))
  expect_identical(.crs_area_error(four), by_hand(four))
})
