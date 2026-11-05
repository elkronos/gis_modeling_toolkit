# tests/testthat/test-fix-mixed-zm-coerce.R
# ---------------------------------------------------------------------------
# coerce_to_points() and the callers that reduce the caller's layer to points
# before anything else reads it.  sf hands GEOS and s2 every geometry of a
# layer in the dimension of its first, so a polygon or line layer that mixes
# XY and XYZ (or XYM) geometries stopped GEOS ("ParseException: Unexpected
# EOF parsing WKB", "Points of LinearRing do not form a closed linestring"),
# stopped st_transform() of its lines ("OGR error"), or was read as rings and
# lines of other points.  rbind() of an XY and an XYZ layer builds such a
# layer without a word.  coerce_to_points() now reads it in x and y and
# returns XY points with a warning of class "spatialkit_zm_dropped"; the
# callers that hand back no copy of the layer flatten it before and after the
# coerce, so they raise no such warning.  Each mixed layer is held to the
# same rows in XY, to the bit.
#
# Also: in mode "auto" a MULTILINESTRING whose longest parts are of equal
# length is given the midpoint lowest in x, then y, whatever the order of
# its parts.
# ---------------------------------------------------------------------------

.mzc_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value, conds = got,
       classes = vapply(got, function(w) class(w)[1L], character(1)),
       texts = vapply(got, conditionMessage, character(1)))
}

# Row `i` of `l` given a third coordinate at its own x and y: Z, or M.
.mzc_row <- function(l, i, dim = "XYZ") {
  g  <- sf::st_geometry(l)
  gi <- sf::st_zm(g[[i]], drop = TRUE, what = "ZM")
  add <- function(m) cbind(m, 3)
  co <- unclass(gi)
  gi <- switch(class(gi)[2L],
    POINT = sf::st_point(c(co, 3), dim = dim),
    LINESTRING = sf::st_linestring(add(co), dim = dim),
    POLYGON = sf::st_polygon(lapply(co, add), dim = dim),
    MULTILINESTRING = sf::st_multilinestring(lapply(co, add), dim = dim),
    MULTIPOLYGON = sf::st_multipolygon(lapply(co, function(p) lapply(p, add)), dim = dim))
  lst <- as.list(g); lst[[i]] <- gi
  sf::st_set_geometry(l, sf::st_sfc(lst, crs = sf::st_crs(l), recompute_bbox = TRUE))
}
.mzc_allz <- function(l) { for (i in seq_len(nrow(l))) l <- .mzc_row(l, i); l }

.mzc_store <- new.env(parent = emptyenv())
.mzc_fixtures <- function() {
  if (is.null(.mzc_store$f)) {
    d <- local({
      withr::local_seed(6901, .rng_kind = "Mersenne-Twister",
                        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
      n <- 150
      x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
      z <- sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 0.3)
      data.frame(x = 5e5 + x, y = 5e6 + y, z = z, a1 = stats::rnorm(n),
                 a2 = 0.5 * z + stats::rnorm(n, 0, 0.3))
    })
    pts   <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
    polys <- sf::st_set_geometry(pts, sf::st_buffer(sf::st_geometry(pts), 10,
                                                    nQuadSegs = 2))
    lines <- sf::st_set_geometry(pts, sf::st_sfc(lapply(seq_len(nrow(d)), function(i)
      sf::st_linestring(rbind(c(d$x[i], d$y[i]), c(d$x[i] + 500, d$y[i] + 800),
                              c(d$x[i] + 900, d$y[i] + 300)))), crs = 32617))
    bb <- sf::st_bbox(pts)
    squares <- sf::st_transform(sf::st_sf(id = 1:16, geometry =
      sf::st_make_grid(sf::st_as_sfc(bb), n = c(4, 4))), 4326)
    .mzc_store$f <- list(polys = polys, lines = lines, squares = squares)
  }
  .mzc_store$f
}

test_that("the premise: GEOS reads a mixed polygon layer in the dimension of its first geometry", {
  F <- .mzc_fixtures()
  p <- F$polys[1:10, ]
  # sf raises GEOS's error from an on.exit() handler, which expect_error()
  # does not see; tryCatch() does.
  geos <- function(l) tryCatch(sf::st_point_on_surface(sf::st_geometry(l)),
                               error = function(e) e)
  expect_s3_class(geos(.mzc_row(p, 1)), "error")
  expect_s3_class(geos(.mzc_row(p, 5)), "error")
  expect_true(.mixes_zm(.mzc_row(p, 5)))
})

test_that("coerce_to_points() reads mixed polygons as the same polygons in XY, with a warning", {
  F <- .mzc_fixtures()
  p <- F$polys[1:40, ]
  for (m in c("auto", "centroid", "point_on_surface")) {
    ref <- .mzc_catch(coerce_to_points(p, m))
    expect_length(ref$classes, 0L)
    for (dim in c("XYZ", "XYM")) for (i in c(1L, 5L)) {
      r <- .mzc_catch(coerce_to_points(.mzc_row(p, i, dim), m))
      info <- sprintf("%s, row %d %s", m, i, dim)
      expect_identical(r$value, ref$value, info = info)
      expect_identical(r$classes, "spatialkit_zm_dropped", info = info)
      expect_identical(r$conds[[1L]]$rows, i, info = info)
      expect_match(r$texts, "^coerce_to_points\\(\\): `x` mixes geometries", info = info)
      expect_match(r$texts, "the points are returned in XY", fixed = TRUE, info = info)
    }
  }
})

test_that("coerce_to_points() reads mixed squares in degrees as the same squares in XY", {
  F <- .mzc_fixtures()
  for (m in c("auto", "centroid")) {
    ref <- .mzc_catch(coerce_to_points(F$squares, m))
    for (i in c(1L, 5L)) {
      r <- .mzc_catch(coerce_to_points(.mzc_row(F$squares, i), m))
      # The centroid of row 5 lay 8,666 km from its place before.
      expect_identical(r$value, ref$value, info = paste(m, i))
      expect_identical(r$classes, c("spatialkit_zm_dropped", ref$classes), info = paste(m, i))
    }
  }
})

test_that("coerce_to_points() reads mixed lines as the same lines in XY", {
  F <- .mzc_fixtures()
  l <- F$lines[1:20, ]
  for (m in c("auto", "line_midpoint", "centroid")) {
    ref <- .mzc_catch(coerce_to_points(l, m))
    for (i in c(1L, 5L)) {
      r <- .mzc_catch(coerce_to_points(.mzc_row(l, i), m))
      expect_identical(r$value, ref$value, info = paste(m, i))
      expect_identical(r$classes, c("spatialkit_zm_dropped", ref$classes), info = paste(m, i))
    }
  }
  # A MULTILINESTRING among them.
  ml <- l[1:6, ]
  g <- sf::st_geometry(ml)
  g2 <- sf::st_sfc(c(list(sf::st_multilinestring(list(unclass(g[[1]]),
                                                      unclass(g[[2]]) + 2000))),
                     as.list(g[-1])), crs = 32617)
  ml <- sf::st_set_geometry(ml, g2)
  ref <- .mzc_catch(coerce_to_points(ml, "auto"))
  r <- .mzc_catch(coerce_to_points(.mzc_row(ml, 1), "auto"))
  expect_identical(r$value, ref$value)
  expect_identical(r$classes, "spatialkit_zm_dropped")
})

test_that("bbox_center reads a mixed layer row by row and gives XY points with no warning", {
  F <- .mzc_fixtures()
  p <- F$polys[1:10, ]
  ref <- .mzc_catch(coerce_to_points(p, "bbox_center"))
  r <- .mzc_catch(coerce_to_points(.mzc_row(p, 1), "bbox_center"))
  expect_identical(r$value, ref$value)
  expect_length(r$classes, 0L)
})

test_that("a layer whose geometry types each keep one dimension comes back in XY, with a true warning", {
  # XY squares beside a POINT Z, or a LINESTRING Z.  Each type is read on
  # its own, so the points were right, but row 4 came back with Z beside XY
  # rows.  The points are now XY throughout, and the warning does not say
  # the layer was read wrongly.
  sq <- function(x0, y0) sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 10, y0), c(x0 + 10, y0 + 10),
                                                   c(x0, y0 + 10), c(x0, y0))))
  g <- lapply(1:30, function(i) sq(5e5 + 100 * i, 5e6 + 37 * i))
  rows4 <- list(sf::st_point(c(5e5 + 400, 5e6 + 500, 12)),
                sf::st_linestring(rbind(c(5e5, 5e6, 1), c(5e5 + 300, 5e6 + 400, 2)), dim = "XYZ"))
  for (r4 in rows4) {
    g[[4]] <- r4
    x <- sf::st_sf(id = 1:30, geometry = sf::st_sfc(g, crs = 32617))
    ref <- coerce_to_points(sf::st_zm(x), "auto")
    r <- .mzc_catch(coerce_to_points(x, "auto"))
    info <- class(r4)[2L]
    expect_identical(r$value, ref, info = info)
    expect_true(all(.geom_dims(r$value) == "XY"), info = info)
    expect_identical(r$classes, "spatialkit_zm_dropped", info = info)
    expect_identical(r$conds[[1L]]$rows, 4L, info = info)
    expect_match(r$texts, "so it is read in x and y: Z and M are dropped", fixed = TRUE, info = info)
    expect_false(grepl("wrongly", r$texts), info = info)
  }
})

test_that("layers in one dimension take the path they took", {
  F <- .mzc_fixtures()
  p <- F$polys[1:20, ]; l <- F$lines[1:20, ]
  pz <- .mzc_allz(p); lz <- .mzc_allz(l)
  for (m in c("auto", "centroid", "point_on_surface")) {
    r <- .mzc_catch(coerce_to_points(pz, m))
    expect_length(r$classes, 0L)
    expect_identical(r$value, coerce_to_points(p, m))
  }
  # The midpoints of lines XYZ throughout keep their Z.
  r <- .mzc_catch(coerce_to_points(lz, "auto"))
  expect_length(r$classes, 0L)
  expect_identical(sf::st_zm(r$value), coerce_to_points(l, "auto"))
  expect_true(all(.geom_dims(r$value) == "XYZ"))
  # An empty polygon of another dimension after the first: no coordinate is
  # misread, so the layer is read as it is, with no warning.
  sq <- function(x0, z = NULL) {
    m <- rbind(c(x0, 0), c(x0 + 10, 0), c(x0 + 10, 10), c(x0, 10), c(x0, 0))
    if (is.null(z)) sf::st_polygon(list(m)) else sf::st_polygon(list(cbind(m, z)), dim = "XYZ")
  }
  a <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(sq(0, 1), sf::st_polygon(), sq(20, 1), crs = 32617))
  b <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(sq(0, 1), sf::st_polygon(dim = "XYZ"), sq(20, 1),
                                                 crs = 32617))
  for (m in c("auto", "centroid", "point_on_surface")) {
    r <- .mzc_catch(coerce_to_points(a, m))
    expect_length(r$classes, 0L)
    expect_identical(r$value$geometry, coerce_to_points(b, m)$geometry)
  }
})

test_that("make_folds() and estimate_sac_range() on mixed polygons run as on the same polygons in XY", {
  F <- .mzc_fixtures()
  p <- F$polys
  ref <- .mzc_catch(make_folds(p, k = 5, method = "random_kfold", seed = 1))
  for (i in c(1L, 5L)) for (dim in c("XYZ", "XYM")) {
    r <- .mzc_catch(make_folds(.mzc_row(p, i, dim), k = 5, method = "random_kfold", seed = 1))
    expect_identical(r$value, ref$value, info = paste(i, dim))
    expect_identical(r$classes, ref$classes, info = paste(i, dim))
  }
  ref <- .mzc_catch(estimate_sac_range(p, "z", seed = 1))
  for (i in c(1L, 5L)) {
    r <- .mzc_catch(estimate_sac_range(.mzc_row(p, i), "z", seed = 1))
    expect_identical(r$value, ref$value, info = i)
    expect_identical(r$classes, ref$classes, info = i)
  }
  # Lines too.
  ref <- .mzc_catch(make_folds(F$lines, k = 5, method = "random_kfold", seed = 1))
  r <- .mzc_catch(make_folds(.mzc_row(F$lines, 1), k = 5, method = "random_kfold", seed = 1))
  expect_identical(r$value, ref$value)
  expect_identical(r$classes, ref$classes)
})

test_that("estimate_sac_range() on mixed polygons identifies the range the same polygons in XY give", {
  # A surface whose range is identified, so the test does not compare NA
  # with NA.
  d <- local({
    withr::local_seed(3, .rng_kind = "Mersenne-Twister",
                      .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    n <- 300
    x <- 5e5 + stats::runif(n, 0, 1e5); y <- 5e6 + stats::runif(n, 0, 1e5)
    data.frame(x = x, y = y, z = sin(x / 6e3) + cos(y / 6e3) + stats::rnorm(n, 0, 0.2))
  })
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  polys <- sf::st_set_geometry(pts, sf::st_buffer(sf::st_geometry(pts), 10, nQuadSegs = 2))
  ref <- .mzc_catch(estimate_sac_range(polys, "z", seed = 1))
  expect_true(is.finite(as.numeric(ref$value)))
  for (i in c(1L, 5L)) {
    r <- .mzc_catch(estimate_sac_range(.mzc_row(polys, i), "z", seed = 1))
    expect_identical(r$value, ref$value, info = i)
    expect_identical(r$classes, ref$classes, info = i)
  }
})

test_that("the callers read a layer XYZ throughout that mixes points with polygons and lines as in XY", {
  # Row 1 a POINT Z, every seventh row a LINESTRING Z, the rest POLYGON Z.
  # The points taken from the polygons and lines are XY beside the POINT Z
  # kept, so make_folds() stopped with "OGR error", and
  # select_features_forward(split) and compare_models_cv() with it.
  F <- .mzc_fixtures()
  p <- F$polys
  co <- sf::st_coordinates(coerce_to_points(p, "centroid"))
  gz <- lapply(seq_len(nrow(p)), function(i) {
    if (i == 1L) return(sf::st_point(c(co[1, 1], co[1, 2], 3)))
    if (i %% 7L == 0L)
      return(sf::st_linestring(rbind(c(co[i, 1], co[i, 2], 3),
                                     c(co[i, 1] + 50, co[i, 2] + 80, 3)), dim = "XYZ"))
    sf::st_polygon(lapply(unclass(sf::st_geometry(p)[[i]]), function(m) cbind(m, 3)), dim = "XYZ")
  })
  lz <- sf::st_set_geometry(p, sf::st_sfc(gz, crs = 32617))
  lx <- sf::st_zm(lz)
  expect_false(.mixes_zm(lz))
  mk <- function(l) make_folds(l, k = 5, method = "random_kfold", seed = 1)
  ref <- .mzc_catch(mk(lx)); r <- .mzc_catch(mk(lz))
  expect_identical(r$value, ref$value)
  expect_identical(r$classes, ref$classes)
  skip_if_not_installed("ranger")
  sel <- function(l) select_features_forward(
    l, "z", c("a1", "a2"), fit_fn = function(tr, v) lm_spatial_fit(tr, "z", v),
    k = 3, quiet = TRUE, select_on = "split")
  ref <- .mzc_catch(sel(lx)); r <- .mzc_catch(sel(lz))
  expect_identical(r$value$split, ref$value$split)
  expect_identical(r$value$selected, ref$value$selected)
  expect_identical(r$classes, ref$classes)
  cmp <- function(l) compare_models_cv(l, "z", c("a1", "a2"), models = "RF", k = 3,
                                       seed = 11, rf_args = list(num_trees = 20),
                                       quiet = TRUE)
  ref <- .mzc_catch(cmp(lx)); r <- .mzc_catch(cmp(lz))
  expect_identical(r$value$cv_results$rf_cv$folds, ref$value$cv_results$rf_cv$folds)
  expect_identical(r$value$comparison, ref$value$comparison)
  expect_identical(r$classes, ref$classes)
})

test_that("make_folds() reads a point layer with an XYZ point first among XY points as in XY", {
  # It stopped with sf's "OGR error" at its Z/M test.
  F <- .mzc_fixtures()
  pts <- coerce_to_points(F$polys, "centroid")
  ref <- .mzc_catch(make_folds(pts, k = 5, method = "random_kfold", seed = 1))
  for (i in c(1L, 5L)) {
    r <- .mzc_catch(make_folds(.mzc_row(pts, i), k = 5, method = "random_kfold", seed = 1))
    expect_identical(r$value, ref$value, info = i)
    expect_identical(r$classes, ref$classes, info = i)
  }
})

test_that("prep_model_data() and the fits on mixed polygons run as on the same polygons in XY", {
  F <- .mzc_fixtures()
  p <- F$polys
  for (pz in c("auto", "centroid")) {
    ref <- .mzc_catch(prep_model_data(p, "z", c("a1", "a2"), pointize = pz))
    for (i in c(1L, 5L)) {
      r <- .mzc_catch(prep_model_data(.mzc_row(p, i), "z", c("a1", "a2"), pointize = pz))
      expect_identical(r$value, ref$value, info = paste(pz, i))
      expect_identical(r$classes, ref$classes, info = paste(pz, i))
    }
  }
  # The lines, with an XYM row.
  ref <- .mzc_catch(prep_model_data(F$lines, "z", c("a1", "a2")))
  r <- .mzc_catch(prep_model_data(.mzc_row(F$lines, 1, "XYM"), "z", c("a1", "a2")))
  expect_identical(r$value, ref$value)
  expect_identical(r$classes, ref$classes)
  # The squares in degrees.
  sq <- F$squares
  sq$z <- seq_len(nrow(sq)) / 4; sq$a1 <- cos(seq_len(nrow(sq)))
  ref <- .mzc_catch(prep_model_data(sq, "z", "a1"))
  for (i in c(1L, 5L)) {
    r <- .mzc_catch(prep_model_data(.mzc_row(sq, i), "z", "a1"))
    expect_identical(r$value, ref$value, info = i)
    expect_identical(r$classes, ref$classes, info = i)
  }
  skip_if_not_installed("ranger")
  fit <- function(l) fit_rf_model(l, "z", c("a1", "a2"), num_trees = 20,
                                  include_coords = TRUE, seed = 1)$model$predictions
  ref <- .mzc_catch(fit(p))
  r <- .mzc_catch(fit(.mzc_row(p, 1)))
  expect_identical(r$value, ref$value)
  expect_false("spatialkit_zm_dropped" %in% r$classes)
})

test_that("compare_models_cv() and select_features_forward(select_on = 'split') on mixed polygons run as in XY", {
  skip_if_not_installed("ranger")
  F <- .mzc_fixtures()
  p <- F$polys
  cmp <- function(l) compare_models_cv(l, "z", c("a1", "a2"), models = "RF", k = 3,
                                       seed = 11, rf_args = list(num_trees = 20),
                                       quiet = TRUE)
  ref <- .mzc_catch(cmp(p))
  for (i in c(1L, 5L)) {
    r <- .mzc_catch(cmp(.mzc_row(p, i)))
    expect_identical(r$value$cv_results$rf_cv$folds, ref$value$cv_results$rf_cv$folds, info = i)
    expect_identical(r$value$comparison, ref$value$comparison, info = i)
    expect_false("spatialkit_zm_dropped" %in% r$classes, info = i)
  }
  sel <- function(l) select_features_forward(
    l, "z", c("a1", "a2"), fit_fn = function(tr, v) lm_spatial_fit(tr, "z", v),
    k = 3, quiet = TRUE, select_on = "split")
  ref <- .mzc_catch(sel(p))
  for (i in c(1L, 5L)) {
    r <- .mzc_catch(sel(.mzc_row(p, i)))
    expect_identical(r$value$split, ref$value$split, info = i)
    expect_identical(r$value$selected, ref$value$selected, info = i)
    expect_identical(r$classes, ref$classes, info = i)
  }
})

test_that("a MULTILINESTRING with equal longest parts gets one point whatever the order of its parts", {
  l1 <- rbind(c(0, 0), c(100, 0)); l2 <- rbind(c(0, 1000), c(100, 1000))
  l3 <- rbind(c(0, 2000), c(150, 2000))
  mk <- function(..., crs = 32617)
    sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_multilinestring(list(...)), crs = crs))
  at <- function(x) unname(sf::st_coordinates(coerce_to_points(x, "auto"))[1, 1:2])
  # Before, (50, 0) with the parts as given and (50, 1000) reversed.
  expect_identical(at(mk(l1, l2)), c(50, 0))
  expect_identical(at(mk(l2, l1)), c(50, 0))
  # Lowest in x first, then in y.
  l4 <- rbind(c(-50, 500), c(50, 500))
  expect_identical(at(mk(l1, l4, l2)), c(0, 500))
  expect_identical(at(mk(l2, l4, l1)), c(0, 500))
  # One longest part: its midpoint, as before.
  expect_identical(at(mk(l1, l3)), c(75, 2000))
  expect_identical(at(mk(l3, l1)), c(75, 2000))
  # Parts of no length.
  d1 <- rbind(c(5, 5), c(5, 5)); d2 <- rbind(c(0, 0), c(0, 0))
  expect_identical(at(mk(d1, d2)), at(mk(d2, d1)))
  # In degrees, through the temporary projection: parts that mirror about
  # the equator have one length there and one midpoint x, so y decides.
  g1 <- rbind(c(-80, 10), c(-79.9, 10)); g2 <- rbind(c(-79.9, -10), c(-80, -10))
  a <- suppressWarnings(at(mk(g1, g2, crs = 4326)))
  b <- suppressWarnings(at(mk(g2, g1, crs = 4326)))
  expect_identical(a, b)
  # The midpoint of the part at 10 S, to the rounding of the projection.
  expect_lt(max(abs(a - c(-79.95, -10))), 1e-4)
})

test_that("congruent longest parts are equal whatever the origin, so the point moves with it", {
  # Two diagonal parts of one length, the second 0.71 to the left and 1000.37
  # up.  At some origins the rounding of the coordinates made the first
  # longer in its last bits and at others the second, so the point taken
  # jumped 1 km with a shift of the origin.
  a <- rbind(c(0.13, 0.29), c(60.17, 80.31))
  b <- sweep(a, 2, c(-0.71, 1000.37), "+")
  withr::local_seed(3, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  dy <- vapply(1:40, function(k) {
    o <- c(stats::runif(1, 1e5, 9e5) + stats::runif(1),
           stats::runif(1, 1e6, 9e6) + stats::runif(1))
    x <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_multilinestring(
      list(sweep(a, 2, o, "+"), sweep(b, 2, o, "+"))), crs = 32617))
    sf::st_coordinates(coerce_to_points(x, "auto"))[1, 2] - o[2]
  }, numeric(1))
  # The second part's midpoint is lower in x: it is taken at every origin.
  expect_true(all(abs(dy - 1040.67) < 1e-3))
})

test_that("?coerce_to_points states the reading of mixed layers and the rule for equal parts", {
  # The source page from the source tree, the installed page under R CMD
  # check, which runs the tests where there is no man/.
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  txt <- page("coerce_to_points.Rd")
  expect_match(txt, paste0("When several parts share the greatest length (to one part in ten ",
                           "million, measured in the projection used), their midpoints are ",
                           "taken and the one lowest in x, then in y, is kept"), fixed = TRUE)
  expect_match(txt, "spatialkit_zm_dropped", fixed = TRUE)
  expect_match(txt, "is read in x and y", fixed = TRUE)
})
