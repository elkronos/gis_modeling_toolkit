# tests/testthat/test-fix-lonlat-surface-warning.R
# ---------------------------------------------------------------------------
# coerce_to_points() takes the point on the surface of a polygon with GEOS,
# on the plane in degrees for a layer in longitude and latitude, and passed
# on sf's warning "st_point_on_surface may not give correct results for
# longitude/latitude data" in sf's words, so every function that reduces
# such a layer raised it, while summarize_by_cell() and cell_values_at()
# muffled it (BUG-222).  On a polygon written across longitude 180 the
# warning was true: GEOS read it the long way round the globe and its point
# lay about 12,000 km from it, and the summaries read it there with nothing
# said (BUG-255).  The warning is now muffled in coerce_to_points() alone,
# and such a polygon gets a point on it.
# ---------------------------------------------------------------------------

# Thirty squares of 1 km in EPSG:32617, moved to EPSG:4326.
.fls_layer <- function() {
  set.seed(21, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  i <- 0:29
  g <- lapply(i, function(k) {
    x0 <- 500000 + (k %% 6) * 1000; y0 <- 4e6 + (k %/% 6) * 1000
    sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 1000, y0), c(x0 + 1000, y0 + 1000),
                              c(x0, y0 + 1000), c(x0, y0))))
  })
  sf::st_transform(sf::st_sf(id = i + 1L, y = sin(i / 3) + stats::rnorm(30, 0, 0.2),
                             x1 = stats::rnorm(30), geometry = sf::st_sfc(g, crs = 32617)),
                   4326)
}

.fls_sf_text <- "st_point_on_surface may not give correct results"

# The value of `expr`; fails when sf's warning is among its warnings.  The
# package's own warnings are let pass without a word.
.fls_no_sf_warning <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  expect_false(any(grepl(.fls_sf_text, w, fixed = TRUE)))
  invisible(v)
}

# Squares of `side` degrees in a row of `n` from longitude `x0`, latitude
# `y0`, each vertex written in -180 to 180: those across 180 are split.
.fls_across <- function(n = 36L, x0 = 179.43, side = 0.08) {
  set.seed(5, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  i <- seq_len(n) - 1L
  g <- lapply(i, function(k) {
    a <- x0 + (k %% 12) * 0.1; b <- 10 + (k %/% 12) * 0.1
    m <- rbind(c(a, b), c(a + side, b), c(a + side, b + side), c(a, b + side), c(a, b))
    m[, 1] <- ifelse(m[, 1] > 180, m[, 1] - 360, m[, 1])
    sf::st_polygon(list(m))
  })
  cell <- (i %% 12) %/% 3 + 4 * (i %/% 12) + 1
  sf::st_sf(poly_id = cell, y = sin(i) + (i %% 12) / 4 + stats::rnorm(n, 0, 0.3),
            geometry = sf::st_sfc(g, crs = 4326))
}

# Whether each point lies on its own polygon, read on the sphere.
.fls_on_own <- function(p, x) {
  gp <- sf::st_geometry(p); gx <- sf::st_geometry(x)
  vapply(seq_along(gx), function(i)
    .with_s2(lengths(sf::st_intersects(gp[i], gx[i])) > 0L), logical(1))
}

# The points of `x`'s squares taken with their longitudes read from 0 to 360,
# brought back into -180 to 180: points on them.
.fls_points_on <- function(x) {
  g <- sf::st_geometry(x)
  un <- sf::st_sfc(lapply(unclass(g), function(s) {
    m <- unclass(s)[[1L]]; m[, 1] <- ifelse(m[, 1] < 0, m[, 1] + 360, m[, 1])
    sf::st_polygon(list(m))
  }))
  p <- sf::st_point_on_surface(un)
  sf::st_set_geometry(x, sf::st_sfc(lapply(p, function(q) {
    if (q[1L] > 180) q[1L] <- q[1L] - 360
    q
  }), crs = 4326))
}

test_that("coerce_to_points() raises no sf warning on lon/lat polygons, and its points do not move", {
  x <- .fls_layer()
  for (s2 in c(TRUE, FALSE)) {
    old <- suppressMessages(sf::sf_use_s2(s2))
    on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
    want <- sf::st_set_geometry(x, suppressWarnings(sf::st_point_on_surface(sf::st_geometry(x))))
    for (m in c("auto", "point_on_surface", "surface")) {
      expect_no_warning(got <- coerce_to_points(x, m))
      expect_identical(got, want, label = paste(s2, m))
    }
    suppressMessages(sf::sf_use_s2(old))
  }
})

test_that("the functions that reduce lon/lat polygons to points raise no sf warning", {
  skip_if_not_installed("ranger")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  x <- .fls_layer()
  for (s2 in c(TRUE, FALSE)) {
    old <- suppressMessages(sf::sf_use_s2(s2))
    on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
    .fls_no_sf_warning(prep_model_data(x, "y", "x1"))
    .fls_no_sf_warning(make_folds(x, k = 3, method = "block_kfold", seed = 1))
    .fls_no_sf_warning(try(estimate_sac_range(x, "y"), silent = TRUE))
    .fls_no_sf_warning(fit_rf_model(x, "y", "x1", num_trees = 20L))
    suppressMessages(sf::sf_use_s2(old))
  }
})

test_that("the package's own warnings on lon/lat polygons are still raised", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  x <- .fls_layer()[1:4, ]
  g <- sf::st_geometry(x)
  g3 <- lapply(unclass(g), function(s)
    structure(lapply(unclass(s), function(m) cbind(m, 1)), class = c("XYZ", "POLYGON", "sfg")))
  mixed <- sf::st_set_geometry(x, sf::st_sfc(c(g3[1], unclass(g)[-1]), crs = 4326))
  w <- list()
  withCallingHandlers(coerce_to_points(mixed, "auto"), warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  cls <- unlist(lapply(w, function(c) class(c)[1L]))
  expect_true("spatialkit_zm_dropped" %in% cls)
  expect_false(any(grepl(.fls_sf_text, vapply(w, conditionMessage, ""), fixed = TRUE)))
})

test_that("a polygon written across longitude 180 gets a point on it", {
  am <- .fls_across()
  span <- vapply(sf::st_geometry(am), function(s) diff(range(unclass(s)[[1L]][, 1L])), 0)
  expect_true(any(span > 180))       # some squares are split
  gos <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(am)))
  expect_false(all(.fls_on_own(gos, am)))  # GEOS's points are not all on them
  for (m in c("auto", "point_on_surface")) {
    expect_no_warning(p <- coerce_to_points(am, m))
    expect_true(all(.fls_on_own(p, am)), label = m)
    # The rows that do not cross keep GEOS's point, to the bit.
    expect_identical(sf::st_geometry(p)[span <= 180], gos[span <= 180], label = m)
    expect_identical(sf::st_geometry(p), sf::st_geometry(.fls_points_on(am)), label = m)
  }
})

test_that("a polygon across 180 whose point cannot be shown on it is said, by row", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  # A bowtie across 180, which s2 cannot read, beside a square.
  x <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(10, 0), c(11, 0), c(11, 1), c(10, 1), c(10, 0)))),
    sf::st_polygon(list(rbind(c(179, 0), c(-179, 1), c(-179, 0), c(179, 1), c(179, 0)))),
    crs = 4326))
  w <- NULL
  p <- withCallingHandlers(coerce_to_points(x, "auto"),
                           spatialkit_surface_point_off = function(c) {
                             w <<- c; invokeRestart("muffleWarning")
                           })
  expect_s3_class(w, "spatialkit_surface_point_off")
  expect_identical(w$rows, 2L)
  expect_match(conditionMessage(w), "^coerce_to_points\\(\\): row 2 of `x` spans more than 180")
  expect_identical(sf::st_geometry(p),
                   suppressWarnings(sf::st_point_on_surface(sf::st_geometry(x))))
})

# --- BUG-255: the summaries' path ---------------------------------------------

test_that("summarize_by_cell() and cell_values_at() read polygons across 180 at points on them", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  am  <- .fls_across()
  pts <- .fls_points_on(am)
  # .pooled_xy(), the path both read: the points, moved back to degrees, lie
  # on their polygons.
  xy <- .pooled_xy(am, NULL, "summarize_by_cell", "`assigned_points_sf`")
  back <- sf::st_transform(sf::st_as_sf(as.data.frame(xy$xy), coords = 1:2, crs = xy$crs), 4326)
  expect_true(all(.fls_on_own(back, am)))
  quiet <- function(expr) suppressWarnings(expr)
  .fls_no_sf_warning(summarize_by_cell(am, "y", estimate = "pooled"))
  expect_identical(quiet(summarize_by_cell(am, "y", estimate = "pooled")),
                   quiet(summarize_by_cell(pts, "y", estimate = "pooled")))
  cells <- sf::st_sf(poly_id = 1:12, geometry = sf::st_sfc(lapply(1:12, function(k) {
    a <- 179.42 + ((k - 1) %% 4) * 0.3; b <- 9.99 + ((k - 1) %/% 4) * 0.1
    m <- rbind(c(a, b), c(a + 0.3, b), c(a + 0.3, b + 0.1), c(a, b + 0.1), c(a, b))
    m[, 1] <- ifelse(m[, 1] > 180, m[, 1] - 360, m[, 1])
    sf::st_polygon(list(m))
  }), crs = 4326))
  vals <- quiet(summarize_by_cell(pts, "y", estimate = "pooled", cells_sf = cells))
  nd <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(179.95, 10.1)),
                                        sf::st_point(c(-179.6, 10.15)), crs = 4326))
  .fls_no_sf_warning(cell_values_at(nd, vals, am, "y"))
  expect_identical(quiet(cell_values_at(nd, vals, am, "y")),
                   quiet(cell_values_at(nd, vals, pts, "y")))
})

test_that("plot_tessellation_map() labels cells across 180 on them", {
  skip_if_not_installed("ggplot2")
  am <- .fls_across(n = 12L, side = 0.08)
  am$cell_id <- seq_len(nrow(am))
  p <- plot_tessellation_map(am, labels = TRUE, label_col = "cell_id")
  lab <- Filter(function(d) "label" %in% names(d), ggplot2::ggplot_build(p)$data)[[1L]]
  at <- sf::st_as_sf(as.data.frame(lab[, c("x", "y")]), coords = c("x", "y"), crs = 4326)
  expect_true(all(.fls_on_own(at, am[match(lab$label, am$cell_id), ])))
})

test_that("the rows said are rows of `x` in a layer that mixes lines and polygons", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  # A line first, so the polygons are rows 2 and 3 of `x` and rows 1 and 2
  # of the polygons alone.
  x <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(10, 5), c(11, 6))),
    sf::st_polygon(list(rbind(c(10, 0), c(11, 0), c(11, 1), c(10, 1), c(10, 0)))),
    sf::st_polygon(list(rbind(c(179, 0), c(-179, 1), c(-179, 0), c(179, 1), c(179, 0)))),
    crs = 4326))
  w <- NULL
  withCallingHandlers(coerce_to_points(x, "auto"),
                      spatialkit_surface_point_off = function(c) {
                        w <<- c; invokeRestart("muffleWarning")
                      })
  expect_identical(w$rows, 3L)
  expect_match(conditionMessage(w), "row 3 of `x`", fixed = TRUE)
})

# --- A layer without a CRS read as EPSG:4326 ----------------------------------

test_that("polygons across 180 in a layer without a CRS read as degrees get points on them", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  am  <- .fls_across(n = 12L)
  pts <- .fls_points_on(am)
  span <- vapply(sf::st_geometry(am), function(s) diff(range(unclass(s)[[1L]][, 1L])), 0)
  moved <- which(span > 180 & !.fls_on_own(
    suppressWarnings(sf::st_point_on_surface(sf::st_geometry(am))), am))
  expect_true(length(moved) > 0L)
  na <- sf::st_set_crs(am, sf::NA_crs_)
  # coerce_to_points(): the points of the layer with EPSG:4326, without its
  # CRS, and a warning that names the rows the reading moves.
  w <- character(0)
  p <- withCallingHandlers(coerce_to_points(na, "auto"), warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  expect_identical(sf::st_geometry(p), sf::st_set_crs(sf::st_geometry(pts), sf::NA_crs_))
  expect_length(w, 1L)
  expect_match(w, "has no CRS and its coordinates look like lon/lat", fixed = TRUE)
  expect_match(w, paste0("row ", moved, " of `x`"), fixed = TRUE)
  # .pooled_xy() (summarize_by_cell(), cell_values_at()) and
  # prep_model_data(), which read such a layer as EPSG:4326: its points,
  # moved back to degrees, lie on their polygons.
  xy <- suppressWarnings(.pooled_xy(na, NULL, "summarize_by_cell", "`assigned_points_sf`"))
  back <- sf::st_transform(sf::st_as_sf(as.data.frame(xy$xy), coords = 1:2, crs = xy$crs), 4326)
  expect_true(all(.fls_on_own(back, am)))
  pm <- suppressWarnings(prep_model_data(cbind(na, x1 = seq_len(nrow(na))), "y", "x1"))
  expect_true(all(.fls_on_own(sf::st_transform(pm, 4326), am)))
})

test_that("a layer without a CRS that is not read as degrees keeps GEOS's points", {
  # A planar layer with a row 400 units wide, whose box leaves the lon/lat
  # envelope: not taken for degrees, so no row is read on the sphere, and
  # nothing is said.
  sq <- function(x0, y0, w, h) sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
  planar <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(
    sq(-100, 10, 400, 10), sq(500, 10, 10, 10)))
  expect_no_warning(p <- coerce_to_points(planar, "auto"))
  expect_identical(sf::st_geometry(p), sf::st_point_on_surface(sf::st_geometry(planar)))
})

test_that("plot_tessellation_map() passes on no sf warning about its label points", {
  skip_if_not_installed("ggplot2")
  # The label points were taken under suppressWarnings(); taken with
  # .surface_points(), sf's "bounding box has potentially an invalid value
  # range for longlat data" for a latitude beyond 90 reached the caller.
  # Built whole, so the box sf stores on the column reads latitude 95.
  sq <- function(x0, y0) sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + 1, y0), c(x0 + 1, y0 + 1), c(x0, y0 + 1), c(x0, y0))))
  g <- lapply(0:5, function(k) sq(-80 + k, 40))
  m <- unclass(g[[2L]])[[1L]]; m[3L, 2L] <- 95
  g[[2L]] <- sf::st_polygon(list(m))
  x <- sf::st_sf(cell_id = 1:6, geometry = sf::st_sfc(g, crs = 4326))
  w <- character(0)
  withCallingHandlers(
    p <- plot_tessellation_map(x, labels = TRUE, label_col = "cell_id"),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("invalid value range for longlat data", w, fixed = TRUE)))
})
