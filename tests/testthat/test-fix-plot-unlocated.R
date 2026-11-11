# tests/testthat/test-fix-plot-unlocated.R
# ---------------------------------------------------------------------------
# BUG-243.  plot_tessellation_map() and plot_folds() drew rows without a
# location: a vertex at Inf gave the map an infinite x range in silence, one
# at 1e200 a map 1e200 wide with ggplot2's bare "Position guide is
# perpendicular to the intended axis", and plot_folds() drew the point
# make_folds() had dropped.  Such rows are now left out with one warning.
# BUG-244.  plot_tessellation_map(labels = TRUE) handed geom_sf_text() a
# label for a cell with an NA ID or an empty geometry, and ggplot2 warned at
# print "Removed 1 row containing missing values".  Such cells are not
# labelled.
# BUG-245.  .transform_or_stamp() stamped a CRS in longitude and latitude on
# CRS-less coordinates outside that range, and sf warned "bounding box has
# potentially an invalid value range for longlat data".  It now refuses.
# BUG-266.  plot_tessellation_map() moved a layer with a CRS with plain
# sf::st_transform(), so a row the plot CRS cannot place was drawn emptied
# or short of a vertex in silence.  It now warns (spatialkit_rows_unplaced).
#
# Layers are built whole with sf::st_sfc(): sf's `[[<-` on a geometry column
# keeps the old bounding box, which hides a far row from st_bbox().
# ---------------------------------------------------------------------------

.pu_sq <- function(x0, y0, s = 100)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                            c(x0, y0 + s), c(x0, y0))))
.pu_good <- function()
  list(.pu_sq(5e5, 4e6), .pu_sq(5e5 + 100, 4e6), .pu_sq(5e5 + 200, 4e6),
       .pu_sq(5e5 + 300, 4e6))
.pu_cells <- function(geoms, crs = 32617)
  sf::st_sf(cell_id = seq_along(geoms), value = seq_along(geoms),
            geometry = sf::st_sfc(geoms, crs = crs))
# Cell 3 with a vertex at x = `at`.
.pu_bad_cell <- function(at)
  sf::st_polygon(list(rbind(c(5e5 + 200, 4e6), c(at, 4e6), c(5e5 + 300, 4e6 + 100),
                            c(5e5 + 200, 4e6 + 100), c(5e5 + 200, 4e6))))
.pu_pts <- function(x2 = 5e5 + 150)
  sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_point(c(5e5 + 50, 4e6 + 50)), sf::st_point(c(x2, 4e6 + 50)),
    sf::st_point(c(5e5 + 250, 4e6 + 50)), crs = 32617))
.pu_bnd <- function(at = NULL) {
  ring <- rbind(c(5e5 - 50, 4e6 - 50), c(5e5 + 450, 4e6 - 50), c(5e5 + 450, 4e6 + 150),
                c(5e5 - 50, 4e6 + 150), c(5e5 - 50, 4e6 - 50))
  if (!is.null(at)) ring[2, 1] <- at
  sf::st_sf(b = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(ring)), .pu_sq(5e5 + 100, 4e6 + 20, 50), crs = 32617))
}

# Built and drawn on a null device, as print() does: ggplot2 warns about
# rows it drops at the draw step.
.pu_build <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  b <- ggplot2::ggplot_build(p)
  ggplot2::ggplot_gtable(b)
  b
}
# Every warning an expression raises, and its value; nothing is muffled
# from the result.
.pu_catch <- function(expr) {
  ws <- list()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws, text = vapply(ws, conditionMessage, ""))
}
.pu_range <- function(b) {
  pp <- b$layout$panel_params[[1]]
  c(pp$x_range, pp$y_range)
}
# The parts of a built plot that are drawn: the layers' data and the
# panel's ranges.
.pu_drawn <- function(b) list(data = b$data, range = .pu_range(b))

test_that("plot_tessellation_map() leaves out a row of any layer at Inf or 1e200, with one warning", {
  skip_if_not_installed("ggplot2")
  good <- .pu_good()
  for (at in c(Inf, 1e200)) {
    lbl <- format(at)
    # tessellation_sf: the map of the other three cells.
    bad <- good; bad[[3]] <- .pu_bad_cell(at)
    r <- .pu_catch(.pu_build(plot_tessellation_map(.pu_cells(bad), fill_col = "value",
                                                   labels = TRUE)))
    expect_identical(r$text, paste0(
      "plot_tessellation_map(): not drawing 1 row(s) of `tessellation_sf` without a ",
      "location (a coordinate that is not finite or lies beyond 1e150)."), info = lbl)
    ref <- .pu_build(plot_tessellation_map(.pu_cells(good)[-3, ], fill_col = "value",
                                           labels = TRUE))
    expect_true(all(is.finite(.pu_range(r$value))), info = lbl)
    expect_identical(.pu_range(r$value), .pu_range(ref), info = lbl)
    expect_identical(nrow(r$value$data[[1]]), 3L, info = lbl)
    expect_identical(r$value$data[[1]]$fill, ref$data[[1]]$fill, info = lbl)
    expect_identical(r$value$data[[2]]$label, c(1L, 2L, 4L), info = lbl)

    # features_sf, seeds_sf and boundary: the map of the layer without it.
    for (arg in c("features_sf", "seeds_sf", "boundary")) {
      lay <- if (arg == "boundary") .pu_bnd(at) else .pu_pts(at)
      a <- list(.pu_cells(good), fill_col = "value"); a[[arg]] <- lay
      r <- .pu_catch(.pu_build(do.call(plot_tessellation_map, a)))
      expect_length(r$warnings, 1L)
      expect_match(r$text, sprintf("not drawing 1 row(s) of `%s` without a location", arg),
                   fixed = TRUE, all = TRUE, info = paste(lbl, arg))
      a[[arg]] <- lay[-if (arg == "boundary") 1L else 2L, ]
      ref <- .pu_build(do.call(plot_tessellation_map, a))
      expect_true(all(is.finite(.pu_range(r$value))), info = paste(lbl, arg))
      expect_identical(.pu_drawn(r$value), .pu_drawn(ref), info = paste(lbl, arg))
    }
  }
})

test_that("rows of several layers are counted in one warning, and a move does not count them again", {
  skip_if_not_installed("ggplot2")
  good <- .pu_good()
  bad <- good; bad[[3]] <- .pu_bad_cell(1e200)
  r <- .pu_catch(.pu_build(plot_tessellation_map(.pu_cells(bad), features_sf = .pu_pts(Inf),
                                                 seeds_sf = .pu_pts(-Inf))))
  expect_identical(r$text, paste0(
    "plot_tessellation_map(): not drawing 1 row(s) of `tessellation_sf`, 1 row(s) of ",
    "`seeds_sf` and 1 row(s) of `features_sf` without a location (a coordinate that is ",
    "not finite or lies beyond 1e150)."))
  # Moved into another CRS (target_crs), the row is left out the same way,
  # and the move raises nothing more.
  r <- .pu_catch(.pu_build(plot_tessellation_map(.pu_cells(bad), target_crs = 32616)))
  expect_length(r$warnings, 1L)
  expect_identical(nrow(r$value$data[[1]]), 3L)
  ref <- .pu_build(plot_tessellation_map(.pu_cells(good)[-3, ], target_crs = 32616))
  expect_identical(.pu_drawn(r$value), .pu_drawn(ref))
  # Every row of a layer without a location: that layer is not drawn.
  r <- .pu_catch(.pu_build(plot_tessellation_map(
    .pu_cells(good), features_sf = sf::st_sf(id = 1, geometry = sf::st_sfc(
      sf::st_point(c(1e200, 4e6)), crs = 32617)))))
  expect_length(r$warnings, 1L)
  expect_length(r$value$data, 1L)
})

test_that("a row a move cannot place is named by its row in the caller's layer after a row is left out", {
  skip_if_not_installed("ggplot2")
  # Row 1 at x = 1e200 has no location; row 2, (10, 0.1), has one that UTM
  # zone 46 cannot place.
  c46 <- .pu_cells(list(.pu_sq(3e5, 0, 2e5), .pu_sq(5e5, 0, 2e5)), crs = 32646)
  fl <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_point(c(1e200, 0.1)), sf::st_point(c(10, 0.1)), sf::st_point(c(94, 0.1)),
    crs = 4326))
  r <- .pu_catch(.pu_build(plot_tessellation_map(c46, features_sf = fl)))
  expect_length(r$warnings, 2L)
  expect_match(r$text[1], "not drawing 1 row(s) of `features_sf` without a location",
               fixed = TRUE)
  w <- r$warnings[[2]]
  expect_s3_class(w, "spatialkit_rows_unplaced")
  expect_identical(w$rows, 2L)
  expect_identical(w$n, 3L)
  # Rows 2 (emptied by the move) and 3 are drawn, as the move leaves them.
  expect_identical(nrow(r$value$data[[2]]), 2L)
})

test_that("plot_folds() leaves out the points make_folds() left out, with one warning", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  xy <- lapply(1:40, function(i) c(5e5 + runif(1, 0, 2000), 4e6 + runif(1, 0, 2000)))
  mk <- function(p7) {
    g <- lapply(xy, sf::st_point); g[[7]] <- p7
    sf::st_sf(id = 1:40, geometry = sf::st_sfc(g, crs = 32617))
  }
  cases <- list(far = sf::st_point(c(1e200, 4e6)), inf = sf::st_point(c(Inf, 4e6)),
                empty = sf::st_point())
  for (nm in names(cases)) {
    pts <- mk(cases[[nm]])
    f <- suppressWarnings(make_folds(pts, k = 4, method = "block_kfold", seed = 1))
    expect_false(7L %in% f$assignment$row_id, info = nm)
    r <- .pu_catch(.pu_build(plot_folds(f, pts)))
    expect_identical(r$text, paste0(
      "plot_folds(): not drawing 1 row(s) of `points_sf` without a location (a coordinate ",
      "that is not finite or lies beyond 1e150",
      if (nm == "empty") paste0(", or a point make_folds() left out of every fold as ",
                                "empty or not placed in the CRS of the folds"),
      ")."), info = nm)
    L <- r$value$plot$layers[[length(r$value$plot$layers)]]$data
    expect_identical(nrow(L), 39L, info = nm)
    expect_false(anyNA(L$.fold), info = nm)
    expect_identical(L$..row_id, setdiff(1:40, 7L), info = nm)
    cc <- do.call(rbind, xy[-7])
    pr <- .pu_range(r$value)
    expect_true(all(is.finite(pr)), info = nm)
    expect_true(pr[1] <= min(cc[, 1]) && pr[2] >= max(cc[, 1]) &&
                pr[2] - pr[1] < 1.2 * diff(range(cc[, 1])), info = nm)
  }
  # A boundary row with a vertex at Inf is left out too, in the same warning.
  pts <- mk(sf::st_point(c(1e200, 4e6)))
  f <- suppressWarnings(make_folds(pts, k = 4, method = "block_kfold", seed = 1))
  bnd <- sf::st_sf(b = 1:2, geometry = sf::st_sfc(
    .pu_sq(5e5 - 100, 4e6 - 100, 2200),
    sf::st_polygon(list(rbind(c(5e5, 4e6), c(Inf, 4e6), c(5e5, 4e6 + 10), c(5e5, 4e6)))),
    crs = 32617))
  r <- .pu_catch(.pu_build(plot_folds(f, pts, boundary = bnd)))
  expect_length(r$warnings, 1L)
  expect_match(r$text, "1 row(s) of `points_sf` and 1 row(s) of `boundary`", fixed = TRUE)
  expect_identical(nrow(r$value$plot$layers[[1]]$data), 1L)
  expect_true(all(is.finite(.pu_range(r$value))))
  # Only a boundary row left out: the words say what that row has, not
  # what a point make_folds() left out has.
  pts <- mk(sf::st_point(c(5e5 + 1, 4e6 + 1)))
  f <- make_folds(pts, k = 4, method = "block_kfold", seed = 1)
  r <- .pu_catch(.pu_build(plot_folds(f, pts, boundary = bnd)))
  expect_identical(r$text, paste0(
    "plot_folds(): not drawing 1 row(s) of `boundary` without a location (a coordinate ",
    "that is not finite or lies beyond 1e150)."))
  # Every point located: nothing is said and every row is drawn.
  pts <- mk(sf::st_point(c(5e5 + 1, 4e6 + 1)))
  f <- make_folds(pts, k = 4, method = "block_kfold", seed = 1)
  r <- .pu_catch(.pu_build(plot_folds(f, pts)))
  expect_length(r$warnings, 0L)
  expect_identical(nrow(r$value$plot$layers[[length(r$value$plot$layers)]]$data), 40L)
})

test_that("plot_folds() leaves out a point make_folds() dropped once projected (latitude 95)", {
  skip_if_not_installed("ggplot2")
  # Point 7 at latitude 95 has a location as given in EPSG:4326, and none
  # in the UTM zone make_folds() builds the blocks in: make_folds() drops
  # it "once projected".  plot_folds() drew it, with blocks emptied in
  # silence at build time, and without blocks at latitude 95 (y range to
  # 98) with no word.
  set.seed(2)
  g <- lapply(1:40, function(i) sf::st_point(c(-80 + runif(1, 0, 0.5), 36 + runif(1, 0, 0.5))))
  g[[7]] <- sf::st_point(c(-80.2, 95))
  pll <- sf::st_sf(id = 1:40, geometry = sf::st_sfc(g, crs = 4326))
  f <- suppressWarnings(make_folds(pll, k = 4, method = "block_kfold", seed = 1))
  expect_false(7L %in% f$assignment$row_id)
  for (blk in c(TRUE, FALSE)) {
    r <- .pu_catch(.pu_build(plot_folds(f, pll, blocks = blk)))
    expect_identical(r$text, paste0(
      "plot_folds(): not drawing 1 row(s) of `points_sf` without a location (a coordinate ",
      "that is not finite or lies beyond 1e150, or a point make_folds() left out of every ",
      "fold as empty or not placed in the CRS of the folds)."), info = paste("blocks", blk))
    L <- r$value$plot$layers[[length(r$value$plot$layers)]]$data
    expect_identical(L$id, setdiff(1:40, 7L), info = paste("blocks", blk))
    expect_false(anyNA(L$.fold), info = paste("blocks", blk))
  }
  # Without blocks the map is drawn in degrees, and ends short of the pole.
  expect_lt(.pu_range(r$value)[4], 37)
  # A located point with no fold (a layer with a row the folds were not
  # built from) is drawn in the colour for NA, as before.
  extra <- rbind(pll[-7, ], sf::st_sf(id = 41L, geometry = sf::st_sfc(
    sf::st_point(c(-80.1, 36.1)), crs = 4326)))
  f2 <- suppressWarnings(make_folds(pll[-7, ], k = 4, method = "block_kfold", seed = 1))
  r <- .pu_catch(.pu_build(plot_folds(f2, extra)))
  expect_length(r$warnings, 0L)
  L <- r$value$plot$layers[[length(r$value$plot$layers)]]$data
  expect_identical(nrow(L), 40L)
  expect_identical(sum(is.na(L$.fold)), 1L)
})

test_that("plot_tessellation_map(labels = TRUE) labels the cells with a label and a location only", {
  skip_if_not_installed("ggplot2")
  good <- .pu_good()
  lab_of <- function(b) {
    l <- Filter(function(d) "label" %in% names(d), b$data)
    if (length(l)) l[[1]] else NULL
  }
  # One NA ID: no ggplot2 warning at print; cells 1, 3 and 4 labelled.
  c1 <- .pu_cells(good); c1$cell_id[2] <- NA
  r <- .pu_catch(.pu_build(plot_tessellation_map(c1, fill_col = "value", labels = TRUE)))
  expect_length(r$warnings, 0L)
  ref <- .pu_build(plot_tessellation_map(.pu_cells(good), fill_col = "value", labels = TRUE))
  expect_identical(lab_of(r$value)$label, c(1L, 3L, 4L))
  expect_identical(lab_of(r$value)[, c("x", "y")], lab_of(ref)[-2, c("x", "y")],
                   ignore_attr = TRUE)
  # Cell 2 empty: the same.
  ge <- good; ge[[2]] <- sf::st_polygon()
  r <- .pu_catch(.pu_build(plot_tessellation_map(.pu_cells(ge), fill_col = "value",
                                                 labels = TRUE)))
  expect_length(r$warnings, 0L)
  expect_identical(lab_of(r$value)$label, c(1L, 3L, 4L))
  # Every ID NA: no labels, the package's log line, no ggplot2 warning.
  c0 <- .pu_cells(good); c0$cell_id <- NA_integer_
  lines <- capture_spatialkit_log(
    r <- .pu_catch(.pu_build(plot_tessellation_map(c0, fill_col = "value", labels = TRUE))),
    level = logger::WARN)
  expect_length(r$warnings, 0L)
  expect_null(lab_of(r$value))
  expect_true(log_has(lines, "no cell has both a label in 'cell_id' and a location; skipping labels"))
  # Every cell with an ID and a location: every cell labelled, as before.
  r <- .pu_catch(.pu_build(plot_tessellation_map(.pu_cells(good), fill_col = "value",
                                                 labels = TRUE)))
  expect_length(r$warnings, 0L)
  expect_identical(lab_of(r$value)$label, 1:4)
})

test_that("a CRS in longitude and latitude is not stamped on CRS-less coordinates outside its range", {
  skip_if_not_installed("ggplot2")
  nc <- .pu_cells(.pu_good()); sf::st_crs(nc) <- NA
  r <- .pu_catch(tryCatch(plot_tessellation_map(nc, target_crs = 4326),
                          error = function(e) e))
  expect_s3_class(r$value, "error")
  expect_length(r$warnings, 0L)
  msg <- conditionMessage(r$value)
  expect_match(msg, paste0(
    "^plot_tessellation_map\\(\\): `tessellation_sf` has no CRS and its coordinates lie ",
    "outside longitude -180 to 360 and latitude -90 to 90 \\(xmin=500000, xmax=500400, ",
    "ymin=4e\\+06, ymax=4\\.0001e\\+06\\), so they cannot be in the target CRS \\('EPSG:4326'\\)"))
  expect_match(msg, "Set the CRS of `tessellation_sf` with sf::st_crs(), or pass a projected `target_crs`.",
               fixed = TRUE)
  # A projected target is stamped with the warning, as before.
  r <- .pu_catch(.pu_build(plot_tessellation_map(nc, target_crs = 32617)))
  expect_length(r$warnings, 1L)
  expect_match(r$text, "stamping the target CRS ('EPSG:32617') WITHOUT reprojection", fixed = TRUE)
  # CRS-less degrees are reprojected, as before.
  ll <- sf::st_transform(.pu_cells(.pu_good()), 4326); sf::st_crs(ll) <- NA
  r <- .pu_catch(.pu_build(plot_tessellation_map(ll, target_crs = 32617)))
  expect_length(r$warnings, 1L)
  expect_match(r$text, "look like lon/lat", fixed = TRUE)
  # Inside the range (longitudes 200 to 201, which sf reads as lon/lat),
  # the stamp is made with its warning, as before.
  e <- sf::st_sfc(.pu_sq(200, 10, 1))
  r <- .pu_catch(.transform_or_stamp(e, 4326, what = "x", caller = "f"))
  expect_length(r$warnings, 1L)
  expect_identical(sf::st_crs(r$value), sf::st_crs(4326))
  # A y beyond 90 north or south is refused even where x is a longitude:
  # latitude 100 to 120, 170 to 175, -120 to -100, and 4e6.  sf adds
  # nothing.
  refused <- "^f\\(\\): `x` has no CRS and its coordinates lie outside longitude"
  for (e in list(.pu_sq(20, 89, 4), .pu_sq(10, 100, 20), .pu_sq(10, 170, 5),
                 .pu_sq(10, -120, 20), .pu_sq(20, 4e6, 100))) {
    r <- .pu_catch(tryCatch(.transform_or_stamp(sf::st_sfc(e), 4326, what = "x", caller = "f"),
                            error = function(e) e))
    expect_s3_class(r$value, "error")
    expect_match(if (inherits(r$value, "error")) conditionMessage(r$value) else "<returned>",
                 refused)
    expect_length(r$warnings, 0L)
  }
  # And an x beyond either end: west of -180 (x = -5e5, y = 10), east of
  # 360; and a y south of -90 alone (x = 10, y = -4e6).
  for (e in list(.pu_sq(-5e5, 10, 1), .pu_sq(400, 10, 1), .pu_sq(10, -4e6, 1)))
    expect_error(.transform_or_stamp(sf::st_sfc(e), 4326, what = "x", caller = "f"), refused)
  # At the bounds, the coordinates are not refused: x from -180 and from
  # 359 to 360, y to 90 and from -90.
  for (e in list(.pu_sq(-180, 10, 1), .pu_sq(359, 10, 1), .pu_sq(200, 89, 1),
                 .pu_sq(200, -90, 1))) {
    r <- .pu_catch(.transform_or_stamp(sf::st_sfc(e), 4326, what = "x", caller = "f"))
    expect_identical(sf::st_crs(r$value), sf::st_crs(4326))
  }
  # A caller that names a row past the pole itself (lat_max = 180) has a y
  # between 90 and 180 stamped, with the warning (test-fix-boundary-degrees.R
  # pins its refusal), and a y beyond 180 refused.
  r <- .pu_catch(.transform_or_stamp(sf::st_sfc(.pu_sq(20, 89, 4)), 4326, what = "x",
                                     caller = "f", lat_max = 180))
  expect_length(r$warnings, 1L)
  expect_identical(sf::st_crs(r$value), sf::st_crs(4326))
  expect_error(.transform_or_stamp(sf::st_sfc(.pu_sq(20, 175, 10)), 4326, what = "x",
                                   caller = "f", lat_max = 180), refused)
  expect_error(.transform_or_stamp(sf::st_sfc(.pu_sq(20, -195, 10)), 4326, what = "x",
                                   caller = "f", lat_max = 180), refused)
  # plot_tessellation_map(target_crs = 4326) on CRS-less cells at latitude
  # 100 to 120, and at 170 to 175: refused, and sf adds nothing.
  for (y0 in c(100, 170)) {
    t1 <- sf::st_sf(cell_id = 1:2, geometry = sf::st_sfc(.pu_sq(10, y0, 5), .pu_sq(30, y0, 5)))
    r <- .pu_catch(tryCatch(plot_tessellation_map(t1, target_crs = 4326),
                            error = function(e) e))
    expect_s3_class(r$value, "error")
    expect_match(if (inherits(r$value, "error")) conditionMessage(r$value) else "<returned>",
                 "^plot_tessellation_map\\(\\): `tessellation_sf` has no CRS")
    expect_length(r$warnings, 0L)
  }
  # Another caller: plot_folds() with points in EPSG:4326 and a CRS-less
  # boundary in metres.
  set.seed(1)
  pll <- sf::st_sf(id = 1:30, geometry = sf::st_sfc(lapply(1:30, function(i)
    sf::st_point(c(-80 + runif(1, 0, 0.02), 36 + runif(1, 0, 0.02)))), crs = 4326))
  f <- suppressWarnings(make_folds(pll, k = 3, method = "random_kfold", seed = 1))
  bnd <- sf::st_sfc(.pu_sq(5e5, 4e6, 3000))
  expect_error(plot_folds(f, pll, boundary = bnd),
               "^plot_folds\\(\\): `boundary` has no CRS and its coordinates lie outside")
})

test_that("plot_tessellation_map() warns about an overlay row the plot CRS cannot place, and draws what it drew", {
  skip_if_not_installed("ggplot2")
  # Two 200 km cells in UTM zone 46; overlays in EPSG:4326 whose row 2
  # (10, 0.1) the zone cannot place, and a ring with a vertex at latitude 95.
  c46 <- .pu_cells(list(.pu_sq(3e5, 0, 2e5), .pu_sq(5e5, 0, 2e5)), crs = 32646)
  fl <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_point(c(92, 0.1)), sf::st_point(c(10, 0.1)), sf::st_point(c(94, 0.1)), crs = 4326))
  ring <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(92, 0), c(95, 0), c(95, 95), c(94, 1.5), c(92, 1.5), c(92, 0)))), crs = 4326))
  # The overlays moved with sf::st_transform(), as the function moved them
  # before: drawn in the plot CRS, they are not moved again.
  pre <- function(x) suppressWarnings(sf::st_transform(x, 32646))
  for (arg in c("features_sf", "seeds_sf", "boundary")) {
    lay <- if (arg == "boundary") ring else fl
    a <- list(c46); a[[arg]] <- lay
    r <- .pu_catch(.pu_build(do.call(plot_tessellation_map, a)))
    expect_length(r$warnings, 1L)
    w <- r$warnings[[1]]
    expect_s3_class(w, "spatialkit_rows_unplaced")
    expect_identical(w$rows, if (arg == "boundary") 1L else 2L, info = arg)
    expect_match(conditionMessage(w), sprintf(
      "^plot_tessellation_map\\(\\): the move of `%s` to EPSG:32646", arg), info = arg)
    a[[arg]] <- pre(lay)
    ref <- .pu_build(do.call(plot_tessellation_map, a))
    expect_identical(.pu_drawn(r$value), .pu_drawn(ref), info = arg)
  }
  # tessellation_sf moved by target_crs.
  cll <- .pu_cells(list(.pu_sq(92, 0, 1), .pu_sq(10, 0, 1)), crs = 4326)
  r <- .pu_catch(.pu_build(plot_tessellation_map(cll, target_crs = 32646)))
  expect_length(r$warnings, 1L)
  expect_s3_class(r$warnings[[1]], "spatialkit_rows_unplaced")
  expect_match(r$text, "the move of `tessellation_sf` to EPSG:32646 cannot place every vertex of row 2 of its 2",
               fixed = TRUE)
  # Overlays the plot CRS places in full, and one already in the plot CRS:
  # no warning, and the same drawing as before.
  for (lay in list(fl[-2, ], pre(fl[-2, ]))) {
    r <- .pu_catch(.pu_build(plot_tessellation_map(c46, features_sf = lay)))
    expect_length(r$warnings, 0L)
    ref <- .pu_build(plot_tessellation_map(c46, features_sf = pre(fl[-2, ])))
    expect_identical(.pu_drawn(r$value), .pu_drawn(ref))
  }
})
