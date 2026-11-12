# tests/testthat/test-fix-unlocated-na.R
# ---------------------------------------------------------------------------
# A row with a coordinate that is not finite or lies beyond 1e150 has no
# location.  The test of a row (.unlocated_features()) read the layer's
# bounding box first and read no vertex when it was finite.  sf::st_bbox()
# skips NA, so a ring or a line with a vertex at (NA, y) passed as located
# (BUG-220); and sf keeps the box it stored on a geometry column when a row
# is replaced with `[[<-`, so a row edited to Inf or 1e200 passed as well
# (BUG-253).  coerce_to_points() gave a polygon with such a vertex a point,
# where a line got an EMPTY POINT, and assign_features_to_polygons() joined
# such a polygon of `polygons_sf`, or stopped GEOS on it (BUG-221).
# kriging_adequacy() repaired such a cell of `cells_sf` into the rest of its
# ring and kriged it.  Each such row is now read as the same row empty.
# ---------------------------------------------------------------------------

.fun_sq <- function(x0, y0, s = 1000)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                            c(x0, y0 + s), c(x0, y0))))

# Twelve squares of 1 km in EPSG:32617.  `bad`: the x of vertex 2 of the ring
# of row `row`; "empty": that row empty.  sf::st_polygon() refuses NA, so the
# vertex is written into the ring's matrix, as a file read with `nan` gives.
.fun_squares <- function(bad = NULL, row = 2L, multi = FALSE) {
  g <- lapply(0:11, function(i) .fun_sq(500000 + (i %% 4) * 1000, 4e6 + (i %/% 4) * 1000))
  if (identical(bad, "empty")) {
    g[[row]] <- sf::st_polygon()
  } else if (!is.null(bad)) {
    r <- unclass(g[[row]]); r[[1L]][2L, 1L] <- bad
    g[[row]] <- structure(r, class = class(g[[row]]))
  }
  if (multi)
    g <- lapply(g, function(s)
      structure(if (length(s)) list(unclass(s)) else list(),
                class = c("XY", "MULTIPOLYGON", "sfg")))
  set.seed(11, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  sf::st_sf(poly_id = 1:12, y = stats::rnorm(12), x1 = stats::rnorm(12),
            geometry = sf::st_sfc(g, crs = 32617))
}

.fun_bad <- list(na = NA_real_, nan = NaN, inf = Inf, far = 1e200)

# Value, warning messages and classes, without letting them through.
.fun_run <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  list(value = v, msgs = vapply(w, conditionMessage, ""),
       classes = unlist(lapply(w, function(c) class(c)[1L])))
}

# --- BUG-220: .unlocated_features() -----------------------------------------

test_that(".unlocated_features() marks a ring or a line with a vertex at NA", {
  sq <- .fun_squares(NA_real_)
  expect_true(all(is.finite(sf::st_bbox(sq))))   # the box skips the NA
  expect_identical(which(.unlocated_features(sf::st_geometry(sq))), 2L)
  expect_identical(which(.unlocated_features(sf::st_geometry(.fun_squares(NA_real_, multi = TRUE)))), 2L)
  # NA in y, and NaN, read from a file as `nan`.
  r <- unclass(.fun_sq(0, 0)); r[[1L]][3L, 2L] <- NA_real_
  ny <- sf::st_sfc(.fun_sq(0, 0), structure(r, class = c("XY", "POLYGON", "sfg")))
  expect_identical(.unlocated_features(ny), c(FALSE, TRUE))
  wkt <- sf::st_as_sfc(c("POLYGON ((0 0, 1 0, 1 1, 0 0))", "POLYGON ((0 0, nan 0, 1 1, 0 0))"))
  expect_identical(.unlocated_features(wkt), c(FALSE, TRUE))
  ln  <- structure(rbind(c(0, 0), c(NA, 1), c(2, 2)), class = c("XY", "LINESTRING", "sfg"))
  mln <- structure(list(rbind(c(0, 0), c(1, 1)), rbind(c(5, 5), c(6, NA))),
                   class = c("XY", "MULTILINESTRING", "sfg"))
  g <- sf::st_sfc(ln, sf::st_linestring(rbind(c(0, 0), c(1, 1))), mln, crs = 32617)
  expect_identical(.unlocated_features(g), c(TRUE, FALSE, TRUE))
  # An empty point holds NA coordinates and is not counted, beside a located
  # point or alone; a layer with no such coordinate has every row located.
  expect_identical(.unlocated_features(sf::st_sfc(sf::st_point(c(1, 2)), sf::st_point())),
                   c(FALSE, FALSE))
  expect_identical(.unlocated_features(sf::st_sfc(sf::st_point())), FALSE)
  # A layer of empty polygons and lines holds no coordinate at all.
  expect_identical(.unlocated_features(sf::st_sfc(sf::st_polygon(), sf::st_linestring())),
                   c(FALSE, FALSE))
  expect_identical(.unlocated_features(sf::st_geometry(.fun_squares())), logical(12))
  # A layer of XY and XYZ points is read point by point (st_coordinates()
  # stops on it).
  mixed <- sf::st_sfc(sf::st_point(c(1, 2)), sf::st_point(c(3, 4, 5)), sf::st_point(),
                      sf::st_point(c(1e200, 1)))
  expect_identical(.unlocated_features(mixed), c(FALSE, FALSE, FALSE, TRUE))
})

test_that("prep_model_data() drops a polygon with a vertex at NA as the same row empty", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  ctl <- .fun_run(prep_model_data(.fun_squares("empty"), "y", "x1"))
  expect_identical(nrow(ctl$value), 11L)
  for (nm in names(.fun_bad)) {
    r <- .fun_run(prep_model_data(.fun_squares(.fun_bad[[nm]]), "y", "x1"))
    expect_identical(r$value, ctl$value, label = nm)
    expect_identical(r$msgs, ctl$msgs, label = nm)
  }
  expect_true("spatialkit_rows_dropped" %in% ctl$classes)
})

test_that("fit_rf_model() fits the rows with a location only", {
  skip_if_not_installed("ranger")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  ctl <- suppressWarnings(fit_rf_model(.fun_squares("empty"), "y", "x1", num_trees = 20L))
  r   <- suppressWarnings(fit_rf_model(.fun_squares(NA_real_), "y", "x1", num_trees = 20L))
  expect_identical(nrow(r$data_sf), 11L)
  expect_identical(stats::fitted(r), stats::fitted(ctl))
})

test_that("make_folds() and fold_separation() read a vertex at NA as an empty row", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  ctl_b <- suppressWarnings(make_folds(.fun_squares("empty"), k = 3, method = "block_kfold", seed = 1))
  ctl_r <- suppressWarnings(make_folds(.fun_squares("empty"), k = 3, method = "random_kfold", seed = 1))
  expect_false(2L %in% ctl_b$assignment$row_id)
  for (nm in names(.fun_bad)) {
    b <- suppressWarnings(make_folds(.fun_squares(.fun_bad[[nm]]), k = 3,
                                     method = "block_kfold", seed = 1))
    r <- suppressWarnings(make_folds(.fun_squares(.fun_bad[[nm]]), k = 3,
                                     method = "random_kfold", seed = 1))
    expect_identical(b$assignment, ctl_b$assignment, label = nm)
    expect_identical(b$folds, ctl_b$folds, label = nm)
    expect_identical(r$assignment, ctl_r$assignment, label = nm)
  }
  # Prediction polygons: the row with the NA vertex takes no part.
  set.seed(4)
  pts <- sf::st_as_sf(data.frame(x = 500000 + stats::runif(40, 0, 4000),
                                 y = 4e6 + stats::runif(40, 0, 3000)),
                      coords = c("x", "y"), crs = 32617)
  f <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
  ctl <- fold_separation(f, pts, sac = 800, prediction_points = .fun_squares("empty"))
  got <- fold_separation(f, pts, sac = 800, prediction_points = .fun_squares(NA_real_))
  expect_identical(got, ctl)
})

test_that("harmonize_crs() returns a polygon with a vertex at NA empty when it moves it", {
  ll <- sf::st_transform(.fun_squares(), 4326)
  h  <- suppressWarnings(harmonize_crs(.fun_squares(NA_real_), ll, prefer = "b"))
  h0 <- suppressWarnings(harmonize_crs(.fun_squares("empty"), ll, prefer = "b"))
  expect_true(sf::st_is_empty(h$a)[2L])
  expect_identical(h$a, h0$a)
})

# --- kriging_adequacy(): a cell without a location ----------------------------

.fun_ka_points <- function() {
  set.seed(3, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  i <- rep(0:11, each = 5)
  p <- sf::st_as_sf(data.frame(x = 500500 + (i %% 4) * 1000 + stats::runif(60, -300, 300),
                               y0 = 4000500 + (i %/% 4) * 1000 + stats::runif(60, -300, 300),
                               poly_id = i + 1L),
                    coords = c("x", "y0"), crs = 32617)
  xy <- sf::st_coordinates(p)
  p$y <- sin(xy[, 1] / 700) + cos(xy[, 2] / 900) + stats::rnorm(60, 0, 0.1)
  p
}

test_that("kriging_adequacy() leaves out a cell with a vertex at NA, Inf or 1e200, as an empty cell", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  pts <- .fun_ka_points()
  cells <- function(bad) .fun_squares(bad)[, "poly_id"]
  ctl <- .fun_run(kriging_adequacy(pts, "y", cells("empty"), k = 3L, nmax = 20L))
  expect_true(is.na(ctl$value$kr_pred[2L]))
  expect_true(any(grepl("1 cell\\(s\\) of `cells_sf` left out .*no location", ctl$msgs)))
  # The empty cell is not described as thin for its bounding box.
  expect_false(any(grepl("max_box_ratio", ctl$msgs)))
  # The report wraps at getOption("width"); print it on lines that do not wrap.
  expect_output(withr::with_options(list(width = 10000), print(ctl$value)),
                "1 left out as empty or without a location")
  for (nm in names(.fun_bad)) {
    r <- .fun_run(kriging_adequacy(pts, "y", cells(.fun_bad[[nm]]), k = 3L, nmax = 20L))
    expect_identical(r$value, ctl$value, label = nm)
    expect_identical(r$msgs, ctl$msgs, label = nm)
  }
})

# --- BUG-253: a box sf keeps after `[[<-` ------------------------------------

.fun_points <- function() {
  set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  sf::st_as_sf(data.frame(x = 5e5 + stats::runif(40, 0, 2000),
                          y = 4e6 + stats::runif(40, 0, 2000),
                          v = stats::rnorm(40), grp = rep(1:8, 5)),
               coords = c("x", "y"), crs = 32617)
}

# The layer with row i edited in place: its stored box is the old one.
.fun_edit <- function(layer, i, s) {
  g <- sf::st_geometry(layer)
  g[[i]] <- s
  sf::st_geometry(layer) <- g
  layer
}

# The same layer built whole: sf computes its box from the geometries.
.fun_whole <- function(layer) {
  g <- sf::st_geometry(layer)
  sf::st_set_geometry(layer, sf::st_sfc(unclass(g), crs = sf::st_crs(g)))
}

test_that("a row edited in place to Inf or 1e200 is read as without a location", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  pts <- .fun_points()
  for (at in list(c(1e200, 4e6), c(Inf, 4e6), c(5e5, -1e200))) {
    lab <- paste(at, collapse = ", ")
    ed <- .fun_edit(pts, 7L, sf::st_point(at))
    expect_identical(sf::st_bbox(ed), sf::st_bbox(pts), label = lab)   # stale
    expect_identical(which(.unlocated_features(sf::st_geometry(ed))), 7L, label = lab)
    wh <- .fun_whole(ed)
    for (m in c("random_kfold", "block_kfold", "leave_location_out")) {
      a <- suppressWarnings(make_folds(ed, k = 4, method = m, seed = 1,
                                       group_var = if (m == "leave_location_out") "grp"))
      b <- suppressWarnings(make_folds(wh, k = 4, method = m, seed = 1,
                                       group_var = if (m == "leave_location_out") "grp"))
      expect_false(7L %in% a$assignment$row_id, label = paste(lab, m))
      expect_identical(a$assignment, b$assignment, label = paste(lab, m))
      expect_identical(a$folds, b$folds, label = paste(lab, m))
    }
  }
  # A polygon edited in place the same way.
  sq <- .fun_squares()
  r  <- unclass(sf::st_geometry(sq)[[5L]]); r[[1L]][2L, 1L] <- 1e200
  ed <- .fun_edit(sq, 5L, structure(r, class = c("XY", "POLYGON", "sfg")))
  expect_true(all(is.finite(sf::st_bbox(ed))) && sf::st_bbox(ed)[["xmax"]] < 1e7)
  expect_identical(which(.unlocated_features(sf::st_geometry(ed))), 5L)
  expect_identical(sf::st_is_empty(coerce_to_points(ed, "auto"))[5L], TRUE)
})

# --- BUG-221: coerce_to_points() and assign_features_to_polygons() ------------

test_that("coerce_to_points() gives a polygon without a location an EMPTY POINT in every mode", {
  modes <- c("auto", "centroid", "point_on_surface", "line_midpoint", "bbox_center")
  for (multi in c(FALSE, TRUE)) for (m in modes) {
    ctl <- coerce_to_points(.fun_squares("empty", multi = multi), m)
    expect_true(sf::st_is_empty(ctl)[2L])
    for (nm in names(.fun_bad)) {
      got <- coerce_to_points(.fun_squares(.fun_bad[[nm]], multi = multi), m)
      expect_identical(got, ctl, label = paste(multi, m, nm))
    }
  }
  # A point at Inf or beyond 1e150 is emptied too, under "auto" as well.
  p <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(sf::st_point(c(1, 2)), sf::st_point(c(Inf, 2)),
                                                 sf::st_point(c(3, 1e200)), crs = 32617))
  for (m in c("auto", "bbox_center"))
    expect_identical(sf::st_is_empty(coerce_to_points(p, m)), c(FALSE, TRUE, TRUE), label = m)
  # The other rows keep their points to the bit.
  expect_identical(sf::st_geometry(coerce_to_points(.fun_squares(), "auto"))[-2L],
                   sf::st_geometry(coerce_to_points(.fun_squares(Inf), "auto"))[-2L])
})

.fun_assign_points <- function() {
  set.seed(3, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  i <- rep(0:11, each = 5)
  sf::st_as_sf(data.frame(x = 500500 + (i %% 4) * 1000 + stats::runif(60, -200, 200),
                          y = 4000500 + (i %/% 4) * 1000 + stats::runif(60, -200, 200)),
               coords = c("x", "y"), crs = 32617)
}

test_that("assign_features_to_polygons() reads a polygon without a location as an empty one", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  pts <- .fun_assign_points()
  for (keep in c(FALSE, TRUE)) for (tb in c("smallest_area", "first")) {
    run <- function(bad) .fun_run(assign_features_to_polygons(
      pts, .fun_squares(bad), keep_unassigned = keep, tie_break = tb, quiet = TRUE))
    ctl <- run("empty")
    expect_identical(sum(ctl$value$poly_id == 2L, na.rm = TRUE), 0L)
    expect_identical(attr(ctl$value, "unassigned")$n, 5L)
    for (nm in names(.fun_bad)) {
      r <- run(.fun_bad[[nm]])
      lab <- paste(keep, tb, nm)
      expect_identical(r$value, ctl$value, label = lab)
      expect_identical(r$msgs, ctl$msgs, label = lab)
      expect_false(any(grepl("NaN/Inf", r$msgs)), label = lab)
    }
  }
  # Polygon features, joined by their largest overlap.
  feats <- sf::st_buffer(pts[1:10, ], 50)
  for (largest in c(TRUE, FALSE)) {
    ctl <- .fun_run(assign_features_to_polygons(feats, .fun_squares("empty"),
                                                largest = largest, quiet = TRUE))
    for (nm in names(.fun_bad)) {
      r <- .fun_run(assign_features_to_polygons(feats, .fun_squares(.fun_bad[[nm]]),
                                                largest = largest, quiet = TRUE))
      expect_identical(r$value, ctl$value, label = paste(largest, nm))
      expect_identical(r$msgs, ctl$msgs, label = paste(largest, nm))
    }
  }
})
