# tests/testthat/test-fix-lonlat-beyond-90.R
# ---------------------------------------------------------------------------
# BUG-223.  A row of a lon/lat data layer with a vertex at a latitude beyond
# 90 is no place on the Earth.  The readers of such a layer read it as
# given: sf's bare "bounding box has potentially an invalid value range for
# longlat data" reached the caller from coerce_to_points(), prep_model_data()
# and the fits on polygons, the k-means seeders and the join's features, and
# a polygon reduced to its box centre (65 N for a square at 35 N) or keyed
# for a stable ID kept the vertex without a word.  These readers now read
# such a row as a row without a location, before anything else reads it, in
# the words they use for one.
# ---------------------------------------------------------------------------

.lb_sq <- function(x0, y0, h = 0.1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + h, y0), c(x0 + h, y0 + h),
                            c(x0, y0 + h), c(x0, y0))))

# 30 points near (-80, 35); `lat95` moves row 1 to latitude 95, `empty`
# gives row 1 an empty point.
.lb_points <- function(row1 = c("given", "lat95", "empty")) {
  row1 <- match.arg(row1)
  set.seed(7)
  n <- 30L
  x <- sf::st_sf(key = seq_len(n), y = stats::rnorm(n), x1 = stats::rnorm(n),
                 x2 = stats::rnorm(n),
                 geometry = sf::st_sfc(lapply(seq_len(n), function(i)
                   sf::st_point(c(-80 + stats::runif(1, 0, 0.5), 35 + stats::runif(1, 0, 0.5)))),
                   crs = 4326))
  g <- sf::st_geometry(x)
  if (row1 == "lat95") g[1] <- sf::st_sfc(sf::st_point(c(g[[1]][1], 95)), crs = 4326)
  if (row1 == "empty") g[1] <- sf::st_sfc(sf::st_point(), crs = 4326)
  sf::st_geometry(x) <- g
  x
}

# 25 squares of 0.1 degree from (-80, 35); `lat95` moves the third vertex
# of row 1 to latitude 95, `empty` gives row 1 an empty polygon.
.lb_polygons <- function(row1 = c("given", "lat95", "empty")) {
  row1 <- match.arg(row1)
  xy <- expand.grid(i = 0:4, j = 0:4)
  set.seed(8)
  x <- sf::st_sf(key = seq_len(nrow(xy)), y = stats::rnorm(nrow(xy)),
                 x1 = stats::rnorm(nrow(xy)), cell_id = rep(1:5, each = 5L),
                 geometry = sf::st_sfc(lapply(seq_len(nrow(xy)), function(k)
                   .lb_sq(-80 + 0.1 * xy$i[k], 35 + 0.1 * xy$j[k])), crs = 4326))
  g <- sf::st_geometry(x)
  if (row1 == "lat95") {
    m <- g[[1]][[1]]; m[3, 2] <- 95
    g[1] <- sf::st_sfc(sf::st_polygon(list(m)), crs = 4326)
  }
  if (row1 == "empty") g[1] <- sf::st_sfc(sf::st_polygon(), crs = 4326)
  sf::st_geometry(x) <- g
  x
}

.LB_SF_TEXT <- "bounding box has potentially an invalid value range for longlat data"

# Every warning an expression raises, and its value; the warnings are
# muffled.
.lb_catch <- function(expr) {
  ws <- list()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws,
       texts = vapply(ws, conditionMessage, character(1)),
       classes = vapply(ws, function(w) class(w)[1L], character(1)))
}
.lb_no_sf_text <- function(r, info = NULL)
  expect_false(any(grepl(.LB_SF_TEXT, r$texts, fixed = TRUE)), info = info)
.lb_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
# sf's caveat that st_point_on_surface() may not give correct results for
# lon/lat data, which ensure_stable_poly_id(method = "surface_point") passes
# on for every lon/lat layer, the latitude aside (BUG-217 owns it), taken
# out of a .lb_catch() result.
.lb_drop_sf_surface <- function(r) {
  keep <- !grepl("st_point_on_surface may not give correct results for longitude/latitude data",
                 r$texts, fixed = TRUE)
  r$warnings <- r$warnings[keep]
  r$texts <- r$texts[keep]
  r$classes <- r$classes[keep]
  r
}

test_that("coerce_to_points() gives a row with a latitude beyond 90 an EMPTY POINT in every mode, and says so", {
  x  <- .lb_polygons("lat95")
  x0 <- .lb_polygons("empty")
  for (s2 in c(TRUE, FALSE)) {
    old <- suppressMessages(sf::sf_use_s2(s2))
    on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
    for (mode in c("auto", "centroid", "point_on_surface", "line_midpoint", "bbox_center")) {
      info <- sprintf("%s, s2 %s", mode, s2)
      r <- .lb_catch(coerce_to_points(x, mode))
      .lb_no_sf_text(r, info)
      expect_identical(r$classes, "spatialkit_lat_beyond_90", info = info)
      w <- r$warnings[[1L]]
      expect_identical(w$rows, 1L, info = info)
      expect_identical(w$n, 25L, info = info)
      expect_match(conditionMessage(w), paste0(
        "^coerce_to_points\\(\\): row 1 of the 25 of `x` has a vertex at a latitude beyond ",
        "90 degrees, which is no place on the Earth, so it has no location and gets an EMPTY POINT"),
        info = info)
      expect_true(sf::st_is_empty(r$value[1, ]), info = info)
      # The row is read as the same row empty: the result is the one the
      # layer with that row empty gives, but for the warning.
      r0 <- .lb_catch(coerce_to_points(x0, mode))
      expect_length(r0$warnings, 0L)
      expect_identical(sf::st_geometry(r$value), sf::st_geometry(r0$value), info = info)
    }
  }
  # Layers with every latitude within 90 raise nothing.
  for (mode in c("auto", "centroid", "point_on_surface", "bbox_center"))
    expect_length(.lb_catch(coerce_to_points(.lb_polygons(), mode))$warnings, 0L)
})

test_that("prep_model_data() drops such a row in its own words, points and polygons, every pointize", {
  for (what in c("points", "polygons")) {
    mk <- if (what == "points") .lb_points else .lb_polygons
    preds <- if (what == "points") c("x1", "x2") else "x1"
    modes <- if (what == "points") "auto" else c("auto", "centroid", "point_on_surface", "bbox_center")
    for (pz in modes) {
      info <- paste(what, pz)
      r  <- .lb_catch(prep_model_data(mk("lat95"), "y", preds, pointize = pz))
      r0 <- .lb_catch(prep_model_data(mk("empty"), "y", preds, pointize = pz))
      .lb_no_sf_text(r, info)
      # The drop line only: no projection warning that it cannot place the
      # row, no spatialkit_lat_beyond_90 from the coerce.
      expect_identical(r$classes, "spatialkit_rows_dropped", info = info)
      expect_match(r$texts, paste0(
        "1 without a location (an empty point, a coordinate that is not finite or lies ",
        "beyond 1e150, or a latitude beyond 90 degrees)"), fixed = TRUE, info = info)
      d <- attr(r$value, "dropped")
      expect_identical(d$which, 1L, info = info)
      expect_identical(d$n_geometry, 1L, info = info)
      expect_identical(d$reason, "geometry", info = info)
      # Read as the layer with that row empty.
      expect_identical(r$value, r0$value, info = info)
      expect_identical(r$texts, r0$texts, info = info)
      # Prediction mode logs the drop and raises nothing.
      rr <- .lb_catch(prep_model_data(mk("lat95"), "y", preds, pointize = pz,
                                      require_response = FALSE))
      expect_length(rr$warnings, 0L)
      expect_identical(attr(rr$value, "dropped")$which, 1L, info = info)
    }
  }
  # A projected layer keeps its words: no latitude is named.
  p <- sf::st_transform(.lb_points("empty"), 32617)
  r <- .lb_catch(prep_model_data(p, "y", c("x1", "x2")))
  expect_match(r$texts, "1 without a location (an empty point, or a coordinate that is not finite or lies beyond 1e150)",
               fixed = TRUE)
  # Within 90, nothing is said.
  expect_length(.lb_catch(prep_model_data(.lb_points(), "y", c("x1", "x2")))$warnings, 0L)
  expect_length(.lb_catch(prep_model_data(.lb_polygons(), "y", "x1"))$warnings, 0L)
})

test_that("the fits on such polygons drop the row with prep_model_data()'s words", {
  skip_if_not_installed("ranger")
  r  <- .lb_catch(fit_rf_model(.lb_polygons("lat95"), "y", "x1", num_trees = 50))
  r0 <- .lb_catch(fit_rf_model(.lb_polygons("empty"), "y", "x1", num_trees = 50))
  .lb_no_sf_text(r)
  expect_true("spatialkit_rows_dropped" %in% r$classes)
  expect_identical(r$texts, r0$texts)
  expect_identical(stats::fitted(r$value), stats::fitted(r0$value))
})

test_that("the k-means seeders drop such a point from the cloud with their own line", {
  pts  <- .lb_points("lat95")
  pts0 <- .lb_points("empty")
  r  <- .lb_catch(voronoi_seeds_kmeans(pts, k = 4, set_seed = 1))
  r0 <- .lb_catch(voronoi_seeds_kmeans(pts0, k = 4, set_seed = 1))
  .lb_no_sf_text(r)
  expect_length(r$warnings, 1L)
  expect_identical(r$texts, paste0(
    "voronoi_seeds_kmeans(): dropping 1 point(s) without a location (an empty point, a ",
    "coordinate that is not finite or lies beyond 1e150, or a latitude beyond 90 degrees)."))
  expect_identical(r$texts, r0$texts)
  expect_identical(r$value, r0$value)

  r  <- .lb_catch(get_voronoi_seeds(method = "kmeans", n = 4, sample_points = pts, set_seed = 1))
  r0 <- .lb_catch(get_voronoi_seeds(method = "kmeans", n = 4, sample_points = pts0, set_seed = 1))
  .lb_no_sf_text(r)
  expect_length(r$warnings, 1L)
  expect_identical(r$texts, paste0(
    "get_voronoi_seeds(kmeans): dropping 1 point(s) without a location (an empty point, a ",
    "coordinate that is not finite or lies beyond 1e150, or a latitude beyond 90 degrees)."))
  expect_identical(r$value, r0$value)
  # A projected cloud keeps its words.
  r <- .lb_catch(voronoi_seeds_kmeans(sf::st_transform(pts0, 32617), k = 4, set_seed = 1))
  expect_identical(r$texts, paste0(
    "voronoi_seeds_kmeans(): dropping 1 point(s) without a location (an empty point, or a ",
    "coordinate that is not finite or lies beyond 1e150)."))
})

test_that("assign_features_to_polygons() leaves such a feature unassigned for want of a location", {
  pts  <- .lb_points("lat95")
  pts0 <- .lb_points("empty")
  # Two cells that cover every point but row 1.
  rect <- function(x0, x1) sf::st_polygon(list(rbind(c(x0, 34.9), c(x1, 34.9), c(x1, 35.6),
                                                     c(x0, 35.6), c(x0, 34.9))))
  cells <- sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(rect(-80.1, -79.75), rect(-79.75, -79.4),
                                                          crs = 4326))
  for (keep in c(FALSE, TRUE)) for (largest in c(TRUE, FALSE)) {
    info <- sprintf("keep_unassigned %s, largest %s", keep, largest)
    r  <- .lb_catch(assign_features_to_polygons(pts, cells, "poly_id",
                                                keep_unassigned = keep, largest = largest))
    r0 <- .lb_catch(assign_features_to_polygons(pts0, cells, "poly_id",
                                                keep_unassigned = keep, largest = largest))
    .lb_no_sf_text(r, info)
    expect_identical(r$texts, r0$texts, info = info)
    if (!keep) {
      expect_identical(r$classes, "spatialkit_unassigned", info = info)
      expect_match(r$texts, paste0(
        "They have no location (an empty point, a coordinate that is not finite or lies ",
        "beyond 1e150, or a latitude beyond 90 degrees)"), fixed = TRUE, info = info)
      expect_identical(r$value, r0$value, info = info)
    } else {
      expect_length(r$warnings, 0L)
      # The row keeps the geometry it arrived with; the rest is as with the
      # row empty.
      expect_identical(sf::st_drop_geometry(r$value), sf::st_drop_geometry(r0$value), info = info)
      expect_identical(attributes(r$value)[c("unassigned", "ties")],
                       attributes(r0$value)[c("unassigned", "ties")], info = info)
    }
  }
  # Cells in a projected CRS: the features are moved to it, and the move of
  # the row with the latitude beyond 90 is not reported as well, beside the
  # unassigned text (the join's copy of the row is emptied before the move).
  for (crs in c(32617, 3857)) for (keep in c(FALSE, TRUE)) {
    info <- sprintf("cells in EPSG:%s, keep_unassigned %s", crs, keep)
    cp <- sf::st_transform(cells, crs)
    r  <- .lb_catch(assign_features_to_polygons(pts, cp, "poly_id", keep_unassigned = keep))
    r0 <- .lb_catch(assign_features_to_polygons(pts0, cp, "poly_id", keep_unassigned = keep))
    .lb_no_sf_text(r, info)
    expect_false("spatialkit_rows_unplaced" %in% r$classes, info = info)
    expect_identical(r$texts, r0$texts, info = info)
    expect_identical(sf::st_drop_geometry(r$value), sf::st_drop_geometry(r0$value), info = info)
    if (!keep) expect_identical(r$classes, "spatialkit_unassigned", info = info)
    else expect_length(r$warnings, 0L)
  }
})

test_that("the pooled readers leave such a row out of the pooled columns in their own words", {
  x  <- .lb_polygons("lat95")
  x0 <- .lb_polygons("empty")
  r  <- .lb_catch(summarize_by_cell(x, "y", id_col = "cell_id", estimate = "pooled"))
  r0 <- .lb_catch(summarize_by_cell(x0, "y", id_col = "cell_id", estimate = "pooled"))
  .lb_no_sf_text(r)
  w <- .lb_of(r, "spatialkit_pool_rows_skipped")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$n, 1L)
  expect_match(conditionMessage(w[[1L]]), paste0(
    "have no location (an empty point, a coordinate that is not finite or lies beyond ",
    "1e150, or a latitude beyond 90 degrees)"), fixed = TRUE)
  expect_identical(r$texts, r0$texts)
  expect_identical(sf::st_drop_geometry(r$value), sf::st_drop_geometry(r0$value))
  # The coordinates the pooled columns are read from: NA for row 1, the
  # rest as for the layer with that row empty.
  p  <- .lb_catch(.pooled_xy(x))
  .lb_no_sf_text(p)
  expect_identical(p$value, .pooled_xy(x0))
  expect_true(all(is.na(p$value$xy[1, ])))
})

test_that("ensure_stable_poly_id() keys such a row as an empty one, returns it as it came, and says so once", {
  x  <- .lb_polygons("lat95")
  x0 <- .lb_polygons("empty")
  for (args in list(list(), list(method = "bbox_center"), list(method = "surface_point"),
                    list(make_valid = FALSE),
                    list(method = "bbox_center", make_valid = FALSE),
                    list(method = "surface_point", make_valid = FALSE))) {
    info <- paste(names(args), unlist(args), collapse = ", ")
    r  <- .lb_drop_sf_surface(.lb_catch(do.call(ensure_stable_poly_id, c(list(x), args))))
    r0 <- .lb_drop_sf_surface(.lb_catch(do.call(ensure_stable_poly_id, c(list(x0), args))))
    .lb_no_sf_text(r, info)
    expect_identical(r$classes, "spatialkit_lat_beyond_90", info = info)
    expect_length(r0$warnings, 0L)
    w <- r$warnings[[1L]]
    expect_identical(w$rows, 1L, info = info)
    expect_identical(w$n, 25L, info = info)
    expect_match(conditionMessage(w), paste0(
      "^ensure_stable_poly_id\\(\\): row 1 of the 25 of `polygons_sf` has a vertex at a ",
      "latitude beyond 90 degrees"), info = info)
    # Keyed as the empty row is: last, and every ID as for the empty row.
    expect_identical(r$value$key, r0$value$key, info = info)
    expect_identical(r$value$poly_id, r0$value$poly_id, info = info)
    expect_identical(r$value$key[nrow(r$value)], 1L, info = info)
    # Returned as it came.
    expect_identical(sf::st_geometry(r$value)[[nrow(r$value)]], sf::st_geometry(x)[[1L]],
                     info = info)
    expect_identical(sf::st_geometry(r$value)[-nrow(r$value)],
                     sf::st_geometry(r0$value)[-nrow(r0$value)], info = info)
  }
  # One polygon, the one with the vertex beyond 90: one ID, with the warning,
  # by every method.  sf::st_point_on_surface() of the sf layer read
  # st_dimension(), NA for the emptied row, and "surface_point" stopped with
  # base R's "missing value where TRUE/FALSE needed".
  for (m in c("centroid", "surface_point", "bbox_center")) for (mv in c(TRUE, FALSE)) {
    info <- sprintf("one row, %s, make_valid %s", m, mv)
    r <- .lb_drop_sf_surface(.lb_catch(ensure_stable_poly_id(x[1, ], method = m, make_valid = mv)))
    expect_identical(r$classes, "spatialkit_lat_beyond_90", info = info)
    expect_identical(r$value$poly_id, 1L, info = info)
    expect_identical(sf::st_geometry(r$value)[[1L]], sf::st_geometry(x)[[1L]], info = info)
  }
  # Within 90, nothing is said.
  expect_length(.lb_catch(ensure_stable_poly_id(.lb_polygons()))$warnings, 0L)
})

test_that("ensure_stable_poly_id() gives IDs to a layer whose every polygon is empty", {
  # sf::st_centroid() of the sf layer read st_dimension(), NA for an empty
  # polygon, and stopped with base R's "missing value where TRUE/FALSE
  # needed".
  # sf::st_point_on_surface() of the sf layer did the same.
  for (crs in list(4326, 32617)) {
    e <- sf::st_sf(key = 1:2, geometry = sf::st_sfc(sf::st_polygon(), sf::st_polygon(), crs = crs))
    for (m in c("centroid", "surface_point", "bbox_center")) for (mv in c(TRUE, FALSE)) {
      info <- sprintf("EPSG:%s, %s, make_valid %s", crs, m, mv)
      out <- suppressWarnings(ensure_stable_poly_id(e, method = m, make_valid = mv))
      expect_identical(out$poly_id, 1:2, info = info)
      expect_identical(out$key, 1:2, info = info)
    }
  }
})

test_that("the latitude is read on the coordinates of a geographic layer only", {
  g <- sf::st_geometry(.lb_polygons("lat95"))
  expect_identical(.unlocated_features(g), rep(FALSE, 25L))
  expect_identical(.unlocated_features(g, lat = TRUE), c(TRUE, rep(FALSE, 24L)))
  # A projected layer, or one without a CRS, with a coordinate of 95 is
  # read as numbers.
  expect_identical(.unlocated_features(suppressWarnings(sf::st_set_crs(g, 32617)), lat = TRUE),
                   rep(FALSE, 25L))
  expect_identical(.unlocated_features(suppressWarnings(sf::st_set_crs(g, NA)), lat = TRUE),
                   rep(FALSE, 25L))
  # A point edited in place: sf keeps the box it stored on the column, which
  # the test does not read.
  p <- sf::st_geometry(.lb_points())
  p[[3]] <- sf::st_point(c(-80, -95))
  expect_true(sf::st_bbox(p)[["ymin"]] > 30)
  expect_identical(which(.unlocated_features(p, lat = TRUE)), 3L)
  expect_identical(which(.beyond_pole_features(p)), 3L)
  # A row already without a location (a longitude beyond 1e150) is not
  # named among the rows a latitude beyond 90 leaves without one, though it
  # has such a latitude too; it still gets an EMPTY POINT.
  q <- .lb_points()[1:4, ]
  gq <- sf::st_geometry(q)
  gq[2] <- sf::st_sfc(sf::st_point(c(-80, 95)), crs = 4326)
  gq[3] <- sf::st_sfc(sf::st_point(c(1e200, 95)), crs = 4326)
  sf::st_geometry(q) <- gq
  expect_identical(.beyond_pole_rows(q), 2L)
  r <- .lb_catch(coerce_to_points(q))
  expect_identical(r$classes, "spatialkit_lat_beyond_90")
  expect_identical(r$warnings[[1L]]$rows, 2L)
  expect_identical(sf::st_is_empty(r$value), c(FALSE, TRUE, TRUE, FALSE))
  # A layer with every latitude within 90 is returned as it is.
  ok <- .lb_polygons()
  expect_identical(.empty_unlocated(ok, lat = TRUE), ok)
  # The movers still compare such a row before and after a move, and say
  # the move cannot place it.
  expect_warning(ensure_projected(.lb_points("lat95")), class = "spatialkit_rows_unplaced")
})
