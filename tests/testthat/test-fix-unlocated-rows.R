# tests/testthat/test-fix-unlocated-rows.R
# ---------------------------------------------------------------------------
# A row whose geometry is empty, or with a coordinate that is not finite or
# lies beyond 1e150, has no location.  prep_model_data() dropped only empty
# and non-finite rows, so a point at a northing of 1e200 was fitted, and the
# same point in EPSG:32616 beside a boundary in EPSG:32617 was moved to a
# finite northing and fitted there.  area_of_applicability() gave such a
# prediction row a location index and refused training data with one as
# "every training row is at the same location".  gp_lengthscale_bounds()
# stopped on an NA row and returned an upper bound of Inf for an Inf or
# 1e200 row.  The bound is now read on the coordinates as given, before any
# move, and again after it, and such rows are read as the same rows empty.
# ---------------------------------------------------------------------------

.ful_layer <- function(n = 120L, seed = 77L) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
  z <- sin(x / 6000) + cos(y / 5000) + stats::rnorm(n, 0, 0.3)
  sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, z = z, wq = stats::rnorm(n)),
               coords = c("x", "y"), crs = 32617)
}

# Row i of `layer`, moved into `epsg`, put at `at` (an empty point for NULL).
.ful_put <- function(layer, epsg, at, i = 1L) {
  l <- sf::st_transform(layer, epsg)
  g <- sf::st_geometry(l)
  g[i] <- sf::st_sfc(if (is.null(at)) sf::st_point() else sf::st_point(at),
                     crs = sf::st_crs(g))
  sf::st_set_geometry(l, g)
}

# Value, warning messages and classes, without letting them through.
.ful_run <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  list(value = v, msgs = vapply(w, conditionMessage, ""),
       classes = unlist(lapply(w, function(c) class(c)[1L])))
}

.ful_far <- list(c(5e5, 1e200), c(1e200, 5e6), c(5e5, -1e200), c(-1e151, 5e6))

# --- BUG-143: prep_model_data() ---------------------------------------------

test_that("prep_model_data() drops a row beyond 1e150 as given, as the same row emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  ctl <- .ful_run(prep_model_data(.ful_put(L, 32617, NULL), "z", "wq"))
  expect_identical(nrow(ctl$value), 119L)
  for (at in .ful_far) {
    r <- .ful_run(prep_model_data(.ful_put(L, 32617, at), "z", "wq"))
    lab <- paste(at, collapse = ", ")
    expect_identical(r$value, ctl$value, label = lab)
    expect_identical(r$msgs, ctl$msgs, label = lab)
    expect_identical(r$classes, ctl$classes, label = lab)
  }
  expect_identical(ctl$classes, "spatialkit_rows_dropped")
  expect_match(ctl$msgs, paste0("0 with a non-finite value and 1 without a location (an ",
                                "empty point, or a coordinate that is not finite or lies ",
                                "beyond 1e150)"), fixed = TRUE)
  dr <- attr(ctl$value, "dropped")
  expect_identical(dr$n_geometry, 1L)
  expect_identical(dr$which, 1L)
  expect_identical(dr$reason, "geometry")
})

test_that("prep_model_data() reads the bound before the move into the boundary's CRS", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(
    c(xmin = 5e5, ymin = 5e6, xmax = 6e5, ymax = 5.1e6), crs = 32617)))
  ctl <- .ful_run(prep_model_data(.ful_put(L, 32616, NULL), "z", "wq", boundary = bnd))
  expect_identical(nrow(ctl$value), 119L)
  expect_identical(sf::st_crs(ctl$value), sf::st_crs(32617))
  for (at in .ful_far[1:2]) {
    r <- .ful_run(prep_model_data(.ful_put(L, 32616, at), "z", "wq", boundary = bnd))
    lab <- paste(at, collapse = ", ")
    expect_identical(r$value, ctl$value, label = lab)
    expect_identical(r$msgs, ctl$msgs, label = lab)
    expect_true(max(abs(sf::st_coordinates(r$value))) < 1e7, label = lab)
  }
})

test_that("prep_model_data() drops a polygon with a vertex beyond 1e150, as the same row emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  sq <- function(x0, y0, far = FALSE) {
    m <- rbind(c(x0, y0), c(x0 + 100, y0), c(x0 + 100, y0 + 100), c(x0, y0 + 100), c(x0, y0))
    if (far) m[3L, 2L] <- 1e200
    sf::st_polygon(list(m))
  }
  xy <- expand.grid(x = 5e5 + 1000 * (0:4), y = 5e6 + 1000 * (0:4))
  set.seed(3)
  P <- sf::st_sf(z = stats::rnorm(25), wq = stats::rnorm(25),
                 geometry = sf::st_sfc(Map(sq, xy$x, xy$y), crs = 32617))
  far <- P; sf::st_geometry(far)[7] <- sf::st_sfc(sq(xy$x[7], xy$y[7], TRUE), crs = 32617)
  emp <- P; sf::st_geometry(emp)[7] <- sf::st_sfc(sf::st_polygon(), crs = 32617)
  a <- .ful_run(prep_model_data(far, "z", "wq"))
  b <- .ful_run(prep_model_data(emp, "z", "wq"))
  expect_identical(a$value, b$value)
  expect_identical(a$msgs, b$msgs)
  expect_identical(attr(a$value, "dropped")$which, 7L)
  # A layer whose rows all have a location loses none.
  ok <- .ful_run(prep_model_data(P, "z", "wq"))
  expect_identical(nrow(ok$value), 25L)
  expect_length(ok$msgs, 0L)
})

test_that("prep_model_data() drops a lon/lat polygon with a vertex beyond 1e150 before its centre is taken", {
  # s2 refuses the polygon with the far vertex, so the row has to be emptied
  # before coerce_to_points() takes its centre: it stopped the call with a
  # refusal naming the row, where the same row emptied is dropped.
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  sq <- function(x0, y0, far = FALSE) {
    m <- rbind(c(x0, y0), c(x0 + .01, y0), c(x0 + .01, y0 + .01), c(x0, y0 + .01), c(x0, y0))
    if (far) m[3L, 2L] <- 1e200
    sf::st_polygon(list(m))
  }
  xy <- expand.grid(x = -80 + 0.05 * (0:4), y = 45 + 0.05 * (0:4))
  set.seed(3)
  P <- sf::st_sf(z = stats::rnorm(25), wq = stats::rnorm(25),
                 geometry = sf::st_sfc(Map(sq, xy$x, xy$y), crs = 4326))
  far <- P; sf::st_geometry(far)[7] <- sf::st_sfc(sq(xy$x[7], xy$y[7], TRUE), crs = 4326)
  emp <- P; sf::st_geometry(emp)[7] <- sf::st_sfc(sf::st_polygon(), crs = 4326)
  for (pz in c("auto", "centroid", "point_on_surface")) {
    a <- .ful_run(prep_model_data(far, "z", "wq", pointize = pz))
    b <- .ful_run(prep_model_data(emp, "z", "wq", pointize = pz))
    expect_identical(a$value, b$value, label = pz)
    expect_identical(a$msgs, b$msgs, label = pz)
    expect_identical(attr(a$value, "dropped")$which, 7L, label = pz)
  }
})

test_that("prep_model_data() on a layer whose rows all have a location keeps every point as given", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  # Coordinates at the bound itself are located.
  g <- sf::st_geometry(L)
  g[2] <- sf::st_sfc(sf::st_point(c(5e5, 1e150)), crs = 32617)
  g[3] <- sf::st_sfc(sf::st_point(c(-1e150, 5e6)), crs = 32617)
  L2 <- sf::st_set_geometry(L, g)
  for (lay in list(L, L2)) {
    r <- .ful_run(prep_model_data(lay, "z", "wq"))
    expect_length(r$msgs, 0L)
    expect_identical(sf::st_coordinates(r$value), sf::st_coordinates(lay))
    expect_identical(attr(r$value, "dropped")$n, 0L)
  }
})

test_that("prep_model_data(require_response = FALSE) records a row beyond 1e150 as a geometry drop", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  a <- prep_model_data(.ful_put(L, 32617, c(5e5, 1e200)), "z", "wq", require_response = FALSE)
  b <- prep_model_data(.ful_put(L, 32617, NULL), "z", "wq", require_response = FALSE)
  expect_identical(a, b)
  expect_identical(attr(a, "dropped")$n_geometry, 1L)
})

test_that("the fits read a row beyond 1e150 as the same row emptied", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("ranger")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  far <- .ful_put(L, 32617, c(5e5, 1e200)); emp <- .ful_put(L, 32617, NULL)
  ga <- suppressWarnings(fit_gwr_model(far, "z", "wq"))
  gb <- suppressWarnings(fit_gwr_model(emp, "z", "wq"))
  expect_identical(nrow(ga$data_sf), 119L)
  expect_identical(ga$data_sf, gb$data_sf)
  expect_identical(ga$info$n_dropped, gb$info$n_dropped)
  expect_identical(stats::fitted(ga), stats::fitted(gb))
  ra <- suppressWarnings(fit_rf_model(far, "z", "wq", num_trees = 20, importance = "none",
                                      include_coords = TRUE, seed = 1))
  rb <- suppressWarnings(fit_rf_model(emp, "z", "wq", num_trees = 20, importance = "none",
                                      include_coords = TRUE, seed = 1))
  expect_identical(ra$data_sf, rb$data_sf)
  expect_identical(stats::predict(ra, L[2:20, ]), stats::predict(rb, L[2:20, ]))
})

test_that("the predict() methods read a newdata row beyond 1e150 before the move into the fit's CRS", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("ranger")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(); tr <- L[1:90, ]; nd <- L[91:120, ]
  fits <- list(
    rf_xy = suppressWarnings(fit_rf_model(tr, "z", "wq", num_trees = 20, importance = "none",
                                          include_coords = TRUE, seed = 1)),
    rf    = suppressWarnings(fit_rf_model(tr, "z", "wq", num_trees = 20, importance = "none",
                                          seed = 1)),
    gwr   = suppressWarnings(fit_gwr_model(tr, "z", "wq")))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    for (epsg in c(32616, 32617)) {
      ctl <- .ful_run(stats::predict(f, .ful_put(nd, epsg, NULL)))
      expect_true(is.na(ctl$value[1]), label = paste(nm, epsg))
      expect_true(all(is.finite(ctl$value[-1])), label = paste(nm, epsg))
      for (at in .ful_far[1:2]) {
        r <- .ful_run(stats::predict(f, .ful_put(nd, epsg, at)))
        lab <- paste(nm, epsg, paste(at, collapse = ", "))
        expect_identical(r$value, ctl$value, label = lab)
        expect_identical(r$msgs, ctl$msgs, label = lab)
      }
    }
  }
})

test_that("predict.bayesian_fit() reads a newdata row beyond 1e150 before the move into the fit's CRS", {
  # No Stan: the engine is an empty brmsfit (no compilation, no draws), and
  # posterior_epred() is replaced by a stand-in whose value is the row's
  # scaled coordinates, so a row predicted where the move put it shows.
  skip_if_not_installed("brms")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(); tr <- L[1:90, ]; nd <- L[91:120, ]
  xy <- sf::st_coordinates(tr)
  cs <- list(x_center = mean(xy[, 1]), x_scale = stats::sd(xy[, 1]),
             y_center = mean(xy[, 2]), y_scale = stats::sd(xy[, 2]))
  df <- sf::st_drop_geometry(tr)
  df$`..x` <- (xy[, 1] - cs$x_center) / cs$x_scale
  df$`..y` <- (xy[, 2] - cs$y_center) / cs$y_scale
  eng <- suppressMessages(brms::brm(
    z ~ wq + gp(..x, ..y, k = 5, c = 1.5, scale = FALSE), data = df, empty = TRUE))
  fit <- new_spatial_fit("bayesian_fit", engine = eng, formula = z ~ wq,
                         response_var = "z", predictor_vars = "wq", data_sf = tr,
                         info = list(coord_scaling = cs))
  local_mocked_bindings(
    posterior_epred = function(object, newdata, ...)
      matrix(rep(newdata$..x + newdata$..y, each = 4L), nrow = 4L),
    .package = "brms")
  ctl <- .ful_run(suppressMessages(stats::predict(fit, .ful_put(nd, 32616, NULL))))
  expect_true(is.na(ctl$value[1]))
  expect_true(all(is.finite(ctl$value[-1])))
  for (at in .ful_far[1:2]) {
    r <- .ful_run(suppressMessages(stats::predict(fit, .ful_put(nd, 32616, at))))
    lab <- paste(at, collapse = ", ")
    expect_identical(r$value, ctl$value, label = lab)
    expect_identical(r$msgs, ctl$msgs, label = lab)
  }
})

test_that("predict_surface() reads a grid or covariates row beyond 1e150 before the move into the fit's CRS", {
  skip_if_not_installed("ranger")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(); tr <- L[1:90, ]; nd <- L[91:120, ]
  rf <- suppressWarnings(fit_rf_model(tr, "z", "wq", num_trees = 20, importance = "none",
                                      include_coords = TRUE, seed = 1))
  # A grid in EPSG:32616, moved into the fit's EPSG:32617: the row comes
  # back empty with an NA prediction, as the same row emptied.
  ctl <- suppressWarnings(predict_surface(rf, grid = .ful_put(nd, 32616, NULL)))
  expect_true(is.na(ctl$.pred[1]) && sf::st_is_empty(ctl)[1])
  for (at in .ful_far[1:2]) {
    s <- suppressWarnings(predict_surface(rf, grid = .ful_put(nd, 32616, at)))
    expect_identical(s, ctl, label = paste(at, collapse = ", "))
  }
  # A grid already in the fit's CRS keeps the row as given, with NA.
  s17 <- suppressWarnings(predict_surface(rf, grid = .ful_put(nd, 32617, c(5e5, 1e200))))
  expect_true(is.na(s17$.pred[1]))
  expect_identical(unname(sf::st_coordinates(s17)[1, ]), c(5e5, 1e200))
  expect_identical(s17$.pred[-1], ctl$.pred[-1])
  # A grid of lon/lat polygons, one with a vertex at latitude 1e200.
  sq <- function(x0, y0, far = FALSE) {
    m <- rbind(c(x0, y0), c(x0 + .01, y0), c(x0 + .01, y0 + .01), c(x0, y0 + .01), c(x0, y0))
    if (far) m[3L, 2L] <- 1e200
    sf::st_polygon(list(m))
  }
  ll <- sf::st_coordinates(sf::st_transform(nd[1:5, ], 4326))
  P <- sf::st_sf(wq = nd$wq[1:5], geometry = sf::st_sfc(Map(sq, ll[, 1], ll[, 2]), crs = 4326))
  Pf <- P; sf::st_geometry(Pf)[2] <- sf::st_sfc(sq(ll[2, 1], ll[2, 2], TRUE), crs = 4326)
  Pe <- P; sf::st_geometry(Pe)[2] <- sf::st_sfc(sf::st_polygon(), crs = 4326)
  a <- suppressWarnings(predict_surface(rf, grid = Pf))
  b <- suppressWarnings(predict_surface(rf, grid = Pe))
  expect_identical(a, b)
  expect_true(is.na(a$.pred[2]))
  # A row of `covariates` moved into the fit's CRS is not the nearest feature
  # of a grid point near where the move put it.
  rf0 <- suppressWarnings(fit_rf_model(tr, "z", "wq", num_trees = 20, importance = "none",
                                       seed = 1))
  moved <- sf::st_coordinates(sf::st_transform(sf::st_sfc(sf::st_point(c(5e5, 1e200)),
                                                          crs = 32616), 32617))
  grid <- sf::st_as_sf(data.frame(x = c(moved[1, 1] + 10, 5.5e5), y = c(moved[1, 2], 5.05e6)),
                       coords = c("x", "y"), crs = 32617)
  cov_far <- .ful_put(L, 32616, c(5e5, 1e200)); cov_far$wq[1] <- 99
  cov_emp <- .ful_put(L, 32616, NULL); cov_emp$wq[1] <- 99
  a <- suppressWarnings(predict_surface(rf0, grid = grid, covariates = cov_far))
  b <- suppressWarnings(predict_surface(rf0, grid = grid, covariates = cov_emp))
  expect_identical(a, b)
  expect_false(any(a$wq == 99))
})

# --- BUG-144: area_of_applicability() --------------------------------------

# Everything the call returns, with the geometry of `aoa` set aside.
.ful_aoa <- function(a) {
  a$aoa <- sf::st_drop_geometry(a$aoa)
  a
}

test_that("area_of_applicability() gives a newdata row beyond 1e150 no location index", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(); tr <- L[1:90, ]; nd <- L[91:120, ]
  ctl <- area_of_applicability(.ful_put(nd, 32617, NULL), train_sf = tr,
                               predictor_vars = "wq", coords = "yes")
  expect_true(is.na(ctl$aoa$DI_location[1]))
  expect_true(is.na(ctl$aoa$AOA[1]))
  for (at in c(.ful_far, list(c(Inf, 5e6)))) {
    a <- area_of_applicability(.ful_put(nd, 32617, at), train_sf = tr,
                               predictor_vars = "wq", coords = "yes")
    expect_identical(.ful_aoa(a), .ful_aoa(ctl), label = paste(at, collapse = ", "))
  }
  # Moved from EPSG:32616 into the training CRS.
  ctl16 <- area_of_applicability(.ful_put(nd, 32616, NULL), train_sf = tr,
                                 predictor_vars = "wq", coords = "yes")
  for (at in .ful_far[1:2]) {
    a <- area_of_applicability(.ful_put(nd, 32616, at), train_sf = tr,
                               predictor_vars = "wq", coords = "yes")
    expect_identical(.ful_aoa(a), .ful_aoa(ctl16), label = paste(at, collapse = ", "))
    expect_true(is.na(a$aoa$DI_location[1]))
  }
})

test_that("area_of_applicability() leaves a training row beyond 1e150 out of the location distances", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(); tr <- L[1:90, ]; nd <- L[91:120, ]
  for (epsg in c(32617, 32616)) {
    ctl <- area_of_applicability(nd, train_sf = .ful_put(tr, epsg, NULL),
                                 predictor_vars = "wq", coords = "yes")
    for (at in .ful_far[1:2]) {
      a <- area_of_applicability(nd, train_sf = .ful_put(tr, epsg, at),
                                 predictor_vars = "wq", coords = "yes")
      expect_identical(.ful_aoa(a), .ful_aoa(ctl),
                       label = paste(epsg, paste(at, collapse = ", ")))
    }
    expect_length(ctl$location$train_DI, 89L)
  }
})

test_that("area_of_applicability() with a forest on the coordinates reads such rows as empty", {
  skip_if_not_installed("ranger")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(); tr <- L[1:90, ]; nd <- L[91:120, ]
  rf <- suppressWarnings(fit_rf_model(tr, "z", "wq", num_trees = 20, importance = "none",
                                      include_coords = TRUE, seed = 1))
  ctl <- area_of_applicability(.ful_put(nd, 32617, NULL), model = rf)
  expect_true(is.na(ctl$aoa$..x[1]) && is.na(ctl$aoa$DI[1]) && is.na(ctl$aoa$AOA[1]))
  for (at in c(.ful_far, list(c(Inf, 5e6)))) {
    a <- area_of_applicability(.ful_put(nd, 32617, at), model = rf)
    expect_identical(.ful_aoa(a), .ful_aoa(ctl), label = paste(at, collapse = ", "))
  }
  # Moved into the training CRS: the row comes back empty, as the emptied row.
  ctl16 <- area_of_applicability(.ful_put(nd, 32616, NULL), model = rf)
  for (at in .ful_far[1:2]) {
    a <- area_of_applicability(.ful_put(nd, 32616, at), model = rf)
    expect_identical(.ful_aoa(a), .ful_aoa(ctl16), label = paste(at, collapse = ", "))
    expect_true(sf::st_is_empty(a$aoa)[1])
    expect_true(sf::st_is_empty(ctl16$aoa)[1])
    expect_identical(sf::st_geometry(a$aoa)[-1], sf::st_geometry(ctl16$aoa)[-1])
  }
  # A training row without a location, beside the forest.
  ctl_tr <- area_of_applicability(nd, model = rf, train_sf = .ful_put(tr, 32617, NULL))
  a <- area_of_applicability(nd, model = rf, train_sf = .ful_put(tr, 32617, c(5e5, 1e200)))
  expect_identical(.ful_aoa(a), .ful_aoa(ctl_tr))
})

test_that("area_of_applicability() names a missing location in its refusals for too few training rows", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  tr <- .ful_put(L[1:2, ], 32617, c(5e5, 1e200))
  expect_error(area_of_applicability(L[3:10, ], train_sf = tr, predictor_vars = "wq",
                                     coords = "yes"),
               "fewer than two training rows with finite values for every predictor and a location.",
               fixed = TRUE)
  # Without a location index no location is named.
  tr2 <- L[1:2, ]; tr2$wq[1] <- NA
  expect_error(area_of_applicability(L[3:10, ], train_sf = tr2, predictor_vars = "wq",
                                     coords = "no"),
               "fewer than two training rows with finite values for every predictor.",
               fixed = TRUE)
})

test_that("area_of_applicability() names a missing location when `folds` meet a training row without one", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(); tr <- L[1:90, ]; nd <- L[91:120, ]
  folds <- rep_len(1:5, 90)
  expect_error(area_of_applicability(nd, train_sf = .ful_put(tr, 32617, c(5e5, 1e200)),
                                     predictor_vars = "wq", coords = "yes", folds = folds),
               "has rows with missing or non-finite predictor values or no location, so `folds`",
               fixed = TRUE)
  tr2 <- tr; tr2$wq[1] <- Inf
  expect_error(area_of_applicability(nd, train_sf = tr2, predictor_vars = "wq",
                                     coords = "no", folds = folds),
               "has rows with missing or non-finite predictor values, so `folds`",
               fixed = TRUE)
})

# --- BUG-146: gp_lengthscale_bounds() --------------------------------------

test_that("gp_lengthscale_bounds() leaves rows without a location out, with a warning that counts them", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  set.seed(1)
  xy <- cbind(stats::runif(100, 0, 1e4), stats::runif(100, 0, 1e4))
  expect_no_warning(ctl <- gp_lengthscale_bounds(xy))
  for (extra in list(c(NA, NA), c(NA, 5), c(Inf, 5), c(5, -Inf), c(NaN, 5),
                     c(5, 1e200), c(-1e151, 5))) {
    lab <- paste(extra, collapse = ", ")
    r <- .ful_run(gp_lengthscale_bounds(rbind(xy, extra)))
    expect_identical(r$value, ctl, label = lab)
    expect_identical(r$classes, "spatialkit_rows_dropped", label = lab)
    expect_identical(r$msgs, paste0(
      "gp_lengthscale_bounds(): 1 of 101 row(s) of `coords_xy` have no location ",
      "(an empty point, or a coordinate that is not finite or lies beyond 1e150); ",
      "the bounds are taken from the other 100."), label = lab)
    # A data frame, and the row first.
    r2 <- .ful_run(gp_lengthscale_bounds(as.data.frame(rbind(extra, xy))))
    expect_identical(r2$value, ctl, label = lab)
  }
  # Coordinates at the bound are located: the row is kept, with no warning,
  # and its distances widen the upper bound.
  for (edge in list(c(1e150, 5), c(5, -1e150))) {
    expect_no_warning(b <- gp_lengthscale_bounds(rbind(xy, edge)))
    expect_true(all(is.finite(b)), label = paste(edge, collapse = ", "))
    expect_gt(b[["upper"]], 1e149)
  }
  # Fewer than two located rows fall to the default.
  r <- .ful_run(gp_lengthscale_bounds(rbind(c(1, 2), c(Inf, 5), c(NA, NA))))
  expect_identical(r$value, c(lower = 0.001, upper = 1))
  expect_match(r$msgs, "2 of 3 row(s)", fixed = TRUE)
  w <- tryCatch(gp_lengthscale_bounds(rbind(xy, c(NA, NA), c(Inf, 1))),
                spatialkit_rows_dropped = function(c) c)
  expect_identical(w$n, 2L)
  expect_identical(w$which, c(101L, 102L))
})

# --- BUG-138, BUG-139, BUG-145: the variogram, the folds, the profile -------
# estimate_sac_range(), make_folds(), resolution_profile(),
# determine_optimal_levels() and fold_separation() dropped only empty or
# non-finite rows, read after the projection or the move into another CRS:
# a row at Inf or beyond 1e150 gave a range of NA for a false reason, a block
# of its own, a place in a test fold, a WSS of NA at every level, or a
# prediction point the move brought back finite.  The bound is now read as
# given and again after the move.

.ful_unloc <- "without a location (an empty point, or a coordinate that is not finite or lies beyond 1e150)"

# The messages of a call on a layer whose row 1 is put at `at`, against those
# of the same row emptied (`ctl`).  A latitude beyond 90 has a location as
# given that the move to the working CRS cannot place, and ensure_projected()
# says so before the caller drops the row (BUG-193).
.ful_msgs_at <- function(ctl, at, n, crs = "EPSG:32617")
  if (identical(at, c(-80, 95)))
    c(sprintf(paste0("ensure_projected(): the move of `x` to %s cannot place every ",
                     "vertex of row 1 of its %d, which has a location in longitude ",
                     "and latitude; it comes back empty or without the vertices it ",
                     "cannot place. Pass a target_crs that holds every location of ",
                     "the layer."), crs, n), ctl) else ctl

# A make_folds() result without the probe of the layer as supplied, which
# records the point a row has as given, so a latitude beyond 90 differs from
# the same row emptied.
.ful_folds <- function(f) { f$params$row_probe <- NULL; f }

# Squares in EPSG:4326 around the points of `layer`; row `i` gets a vertex at
# `bad`, or is emptied.
.ful_squares <- function(layer, bad = NULL, empty = FALSE, i = 7L, d = 0.001) {
  ll <- sf::st_transform(layer, 4326)
  xy <- sf::st_coordinates(ll)
  sq <- function(x, y, v = NULL) {
    m <- rbind(c(x - d, y - d), c(x + d, y - d), c(x + d, y + d), c(x - d, y + d),
               c(x - d, y - d))
    if (!is.null(v)) m[3, ] <- v
    sf::st_polygon(list(m))
  }
  g <- lapply(seq_len(nrow(xy)), function(j) sq(xy[j, 1], xy[j, 2]))
  if (!is.null(bad)) g[[i]] <- sq(xy[i, 1], xy[i, 2], bad)
  if (empty) g[[i]] <- sf::st_polygon()
  sf::st_sf(sf::st_drop_geometry(ll), geometry = sf::st_sfc(g, crs = 4326))
}

test_that("estimate_sac_range() drops a row at Inf or beyond 1e150 as given, as the same row emptied", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  for (epsg in c(32617, 32616, 4326)) {
    ctl <- estimate_sac_range(.ful_put(L, epsg, NULL), "z", seed = 1)
    expect_true(is.finite(ctl), label = paste("control", epsg))
    ats <- if (epsg == 4326) list(c(-80, 1e200), c(1e200, 45), c(Inf, 45))
           else c(list(c(Inf, 5e6)), .ful_far)
    for (at in ats) {
      lab <- paste(epsg, paste(at, collapse = ", "))
      lines <- capture_spatialkit_log(r <- estimate_sac_range(.ful_put(L, epsg, at), "z", seed = 1))
      expect_identical(r, ctl, label = lab)
      expect_true(any(grepl(paste0("estimate_sac_range(): dropping 1 point(s) ", .ful_unloc, "."),
                            lines, fixed = TRUE)), label = lab)
    }
  }
})

test_that("estimate_sac_range() and make_folds() drop a row that loses its location in the projection", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  ctl_l <- .ful_put(L, 4326, NULL)
  bad_l <- suppressWarnings(.ful_put(L, 4326, c(-80, 95)))   # a latitude beyond 90
  ctl <- estimate_sac_range(ctl_l, "z", seed = 1)
  lines <- capture_spatialkit_log(r <- suppressWarnings(estimate_sac_range(bad_l, "z", seed = 1)))
  expect_identical(r, ctl)
  expect_true(any(grepl("estimate_sac_range(): dropping 1 point(s) without a location once projected (an empty point",
                        lines, fixed = TRUE)))
  grid <- sf::st_transform(.ful_layer(40L, seed = 5L), 4326)
  cases <- list(block_kfold  = list(k = 4),
                buffered_loo = list(buffer = 5000),
                nndm         = list(prediction_points = grid),
                knndm        = list(k = 4, prediction_points = grid))
  for (m in names(cases)) {
    args <- c(list(method = m, seed = 1), cases[[m]])
    a <- do.call(make_folds, c(list(points_sf = ctl_l), args))
    lines <- capture_spatialkit_log(
      b <- suppressWarnings(do.call(make_folds, c(list(points_sf = bad_l), args))))
    expect_identical(.ful_folds(b), .ful_folds(a), label = m)
    expect_true(any(grepl(sprintf("make_folds(%s): dropping 1 point(s) without a location once projected (an empty point",
                                  m), lines, fixed = TRUE)), label = m)
  }
})

test_that("make_folds() drops a row at Inf or beyond 1e150 as given, under every method", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  grid <- .ful_layer(40L, seed = 5L)
  cases <- list(random_kfold = list(k = 5),
                block_kfold  = list(k = 5),
                buffered_loo = list(buffer = 5000),
                nndm         = list(prediction_points = grid),
                knndm        = list(k = 4, prediction_points = grid))
  for (epsg in c(32617, 32616, 4326)) {
    ats <- if (epsg == 4326) list(c(-80, 1e200), c(Inf, 45))
           else c(list(c(Inf, 5e6)), .ful_far)
    for (m in names(cases)) {
      args <- c(list(method = m, seed = 1), cases[[m]])
      ctl <- do.call(make_folds, c(list(points_sf = .ful_put(L, epsg, NULL)), args))
      expect_identical(nrow(ctl$assignment), 119L, label = m)
      expect_false(1L %in% ctl$assignment$row_id, label = m)
      for (at in ats) {
        lab <- paste(m, epsg, paste(at, collapse = ", "))
        lines <- capture_spatialkit_log(
          r <- do.call(make_folds, c(list(points_sf = .ful_put(L, epsg, at)), args)))
        expect_identical(r, ctl, label = lab)
        expect_true(any(grepl(paste0("make_folds(): dropping 1 point(s) ", .ful_unloc, "."),
                              lines, fixed = TRUE)), label = lab)
      }
    }
  }
  # Coordinates at the bound have a location: the row is kept.
  for (edge in list(c(5e5, 1e150), c(-1e150, 5e6))) {
    f <- make_folds(.ful_put(L, 32617, edge), k = 5, method = "random_kfold", seed = 1)
    expect_identical(nrow(f$assignment), 120L, label = paste(edge, collapse = ", "))
  }
})

test_that("make_folds() and estimate_sac_range() read a lon/lat polygon with a vertex beyond 1e150 as the same row emptied", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  ctl_l <- .ful_squares(L, empty = TRUE)
  ctl_r <- suppressWarnings(estimate_sac_range(ctl_l, "z", seed = 1))
  for (bad in list(c(-80, 1e200), c(Inf, 45))) {
    lab <- paste(bad, collapse = ", ")
    bad_l <- .ful_squares(L, bad = bad)
    expect_identical(suppressWarnings(estimate_sac_range(bad_l, "z", seed = 1)), ctl_r,
                     label = lab)
    for (m in c("random_kfold", "block_kfold")) {
      a <- suppressWarnings(make_folds(ctl_l, k = 5, method = m, seed = 1))
      b <- suppressWarnings(make_folds(bad_l, k = 5, method = m, seed = 1))
      expect_identical(b, a, label = paste(m, lab))
    }
  }
})

test_that("make_folds(nndm, knndm) reads a prediction point beyond 1e150 as given, before the move", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  tr <- .ful_layer()
  pp <- .ful_layer(60L, seed = 5L)
  for (m in c("nndm", "knndm")) {
    for (epsg in c(32617, 32616)) {
      ctl <- make_folds(tr, k = 5, method = m, seed = 1,
                        prediction_points = .ful_put(pp, epsg, NULL))
      # The emptied point takes no part in the target and is not counted.
      expect_identical(ctl$params$n_prediction_points, 59L, label = m)
      expect_length(ctl$params$target_distances, 59L)
      for (at in c(list(c(Inf, 5e6)), .ful_far)) {
        lab <- paste(m, epsg, paste(at, collapse = ", "))
        r <- .ful_run(make_folds(tr, k = 5, method = m, seed = 1,
                                 prediction_points = .ful_put(pp, epsg, at)))
        expect_identical(r$value, ctl, label = lab)
        expect_length(r$msgs, 0L)
      }
    }
  }
  # Prediction points that all have a location are all counted.
  f <- make_folds(tr, k = 5, method = "nndm", seed = 1, prediction_points = pp)
  expect_identical(f$params$n_prediction_points, 60L)
})

test_that("fold_separation() reads a prediction point beyond 1e150 as given, and does not count it", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  tr <- .ful_layer()
  pp <- .ful_layer(60L, seed = 5L)
  f <- make_folds(tr, k = 5, method = "random_kfold", seed = 1)
  for (epsg in c(32617, 32616, 4326)) {
    ctl <- fold_separation(f, tr, sac = 20000, prediction_points = .ful_put(pp, epsg, NULL))
    expect_identical(attr(ctl, "target")$n_prediction_points, 59L)
    ats <- if (epsg == 4326) list(c(-80, 1e200), c(-80, 95), c(Inf, 45))
           else c(list(c(Inf, 5e6)), .ful_far)
    for (at in ats) {
      lab <- paste(epsg, paste(at, collapse = ", "))
      r <- suppressWarnings(fold_separation(f, tr, sac = 20000,
                                            prediction_points = .ful_put(pp, epsg, at)))
      expect_identical(r, ctl, label = lab)
      expect_true(any(grepl("prediction target: 59 point(s)", capture.output(print(r)),
                            fixed = TRUE)), label = lab)
    }
  }
})

test_that("fold_separation() reads a data row beyond 1e150 as given, before the move into the range's CRS", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  sac <- estimate_sac_range(L, "z", seed = 1)
  # Bare splits that name row 1, as folds built on the layer with the row in
  # place name it.
  splits <- make_folds(L, k = 5, method = "random_kfold", seed = 1)$folds
  for (epsg in c(32617, 32616)) {
    ctl <- fold_separation(splits, .ful_put(L, epsg, NULL), sac = sac)
    for (at in c(list(c(Inf, 5e6)), .ful_far)) {
      lab <- paste(epsg, paste(at, collapse = ", "))
      r <- fold_separation(splits, .ful_put(L, epsg, at), sac = sac)
      expect_identical(r, ctl, label = lab)
    }
  }
})

test_that("resolution_profile() drops a row beyond 1e150 as given, before the projection", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(150L)
  for (epsg in c(32617, 4326)) {
    ctl <- .ful_run(resolution_profile(.ful_put(L, epsg, NULL), "z", sac = NULL, n_levels = 3))
    expect_true(all(is.finite(ctl$value$wss)))
    ats <- if (epsg == 4326) list(c(-80, 1e200), c(-80, 95), c(Inf, 45))
           else c(list(c(Inf, 5e6)), .ful_far)
    for (at in ats) {
      lab <- paste(epsg, paste(at, collapse = ", "))
      r <- .ful_run(withCallingHandlers(
        resolution_profile(.ful_put(L, epsg, at), "z", sac = NULL, n_levels = 3),
        warning = function(w)
          if (grepl("bounding box has potentially an invalid value range",
                    conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")))
      expect_identical(r$value, ctl$value, label = lab)
      expect_identical(r$msgs, .ful_msgs_at(ctl$msgs, at, 150L), label = lab)
    }
    expect_true(any(grepl(paste0("resolution_profile(): dropping 1 point(s) ", .ful_unloc, "."),
                          ctl$msgs, fixed = TRUE)))
  }
  # A lon/lat polygon with a vertex beyond 1e150 or at Inf is the same row
  # emptied.
  ctl <- .ful_run(resolution_profile(.ful_squares(L, empty = TRUE), "z", n_levels = 3))
  for (bad in list(c(-80, 1e200), c(Inf, 45))) {
    r <- .ful_run(resolution_profile(.ful_squares(L, bad = bad), "z", n_levels = 3))
    expect_identical(r$value, ctl$value, label = paste(bad, collapse = ", "))
    expect_identical(r$msgs, ctl$msgs, label = paste(bad, collapse = ", "))
  }
})

test_that("resolution_profile() reads the bound before and after the move into the CRS of `sac`", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(150L)
  sac <- estimate_sac_range(L, "z", seed = 1)
  expect_identical(attr(sac, "crs"), sf::st_crs(32617))
  # A row beyond 1e150 as given, in another projected CRS: dropped before the
  # move, which would bring it back finite.
  for (epsg in c(32616, 3857)) {
    ctl <- .ful_run(resolution_profile(.ful_put(L, epsg, NULL), "z", sac = sac, n_levels = 3))
    r <- .ful_run(resolution_profile(.ful_put(L, epsg, c(5e5, 1e200)), "z", sac = sac,
                                     n_levels = 3))
    expect_identical(r$value, ctl$value, label = paste(epsg))
    expect_identical(r$msgs, ctl$msgs, label = paste(epsg))
  }
  # A row with a location in EPSG:5070 that has none in EPSG:32617 (the move
  # gives NaN): dropped after the move, with the drop warning, where it
  # stopped the call.  The move says first that it cannot place the row
  # (BUG-203, class spatialkit_rows_unplaced, in the words of
  # resolution_profile() and naming EPSG:5070).
  ctl <- .ful_run(resolution_profile(.ful_put(L, 5070, NULL), "z", sac = sac, n_levels = 3))
  r <- .ful_run(resolution_profile(.ful_put(L, 5070, c(3e7, 3e7)), "z", sac = sac,
                                   n_levels = 3))
  expect_identical(r$value, ctl$value)
  expect_identical(r$msgs[-(1:2)], ctl$msgs[-1L])
  expect_identical(r$classes[1L], "spatialkit_rows_unplaced")
  expect_identical(r$msgs[1L], paste0("resolution_profile(): the move of `data_sf` to ",
                                      "EPSG:32617 cannot place every vertex of row 1 of its ",
                                      "150, which has a location in EPSG:5070; it comes back ",
                                      "empty or without the vertices it cannot place. ",
                                      "Correct or drop those rows."))
  expect_identical(r$msgs[2L], paste0("resolution_profile(): dropping 1 point(s) without a ",
                                      "location in the CRS of `sac` (an empty point, or a ",
                                      "coordinate that is not finite or lies beyond 1e150)."))
  # With the split, a row the move leaves without a location is in neither
  # half and the halves index the layer as passed.  The move is made to
  # leave row 5 at 1e200 (a far row as given would stop the split itself),
  # with row 1 empty as given, so both mappings are read.
  ctl <- suppressWarnings(resolution_profile(.ful_put(L, 32616, NULL), "z", sac = sac,
                                             n_levels = 3, select_on = "split", seed = 1))
  real <- spatialkit:::.transform_or_stamp
  local_mocked_bindings(.transform_or_stamp = function(x, crs, what = NULL, caller = NULL, ...) {
    out <- real(x, crs, what = what, caller = caller, ...)
    if (identical(what, "data_sf") && identical(caller, "resolution_profile")) {
      g <- sf::st_geometry(out)
      g[4] <- sf::st_sfc(sf::st_point(c(5e5, 1e200)), crs = sf::st_crs(g))
      out <- sf::st_set_geometry(out, g)
    }
    out
  })
  r <- .ful_run(resolution_profile(.ful_put(L, 32616, NULL), "z", sac = sac,
                                   n_levels = 3, select_on = "split", seed = 1))
  sp0 <- attr(ctl, "split")
  sp <- attr(r$value, "split")
  # Row 5 of the layer is row 4 after the empty row 1 is dropped.
  expect_true(any(grepl("dropping 1 point(s) without a location in the CRS of `sac`",
                        r$msgs, fixed = TRUE)))
  expect_false(any(c(1L, 5L) %in% c(sp$selection, sp$estimation)))
  expect_setequal(c(sp$selection, sp$estimation), setdiff(1:150, c(1L, 5L)))
  expect_length(intersect(sp$selection, sp$estimation), 0L)
  expect_identical(sp$selection, setdiff(sp0$selection, 5L))
  expect_identical(sp$estimation, setdiff(sp0$estimation, 5L))
})

test_that("resolution_profile() boxes the halves of the split again when the move into the CRS of `sac` drops a row", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(150L)
  sac <- estimate_sac_range(L, "z", seed = 1)
  given <- .ful_put(L, 32616, NULL)
  # The row with the largest x: its half's box is wider with it than without.
  p <- 1L + which.max(sf::st_coordinates(given)[-1L, 1L])
  real <- spatialkit:::.transform_or_stamp
  local_mocked_bindings(.transform_or_stamp = function(x, crs, what = NULL, caller = NULL, ...) {
    out <- real(x, crs, what = what, caller = caller, ...)
    if (identical(what, "data_sf") && identical(caller, "resolution_profile")) {
      g <- sf::st_geometry(out)
      g[p - 1L] <- sf::st_sfc(sf::st_point(c(5e5, 1e200)), crs = sf::st_crs(g))
      out <- sf::st_set_geometry(out, g)
    }
    out
  })
  r <- .ful_run(resolution_profile(given, "z", sac = sac, n_levels = 3,
                                   select_on = "split", seed = 1))
  sp <- attr(r$value, "split")
  expect_false(p %in% c(sp$selection, sp$estimation))
  # The boxes are of the rows each half holds, in data_sf's CRS, and the
  # balance is of their counts: they were of the halves before the drop.
  g <- sf::st_geometry(given)
  expect_identical(sp$extent$selection, sf::st_bbox(g[sp$selection]))
  expect_identical(sp$extent$estimation, sf::st_bbox(g[sp$estimation]))
  n <- c(length(sp$selection), length(sp$estimation))
  expect_identical(sp$balance, max(n) / min(n))
})

test_that("the fold probe leaves a row beyond 1e150 out", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  # Without a CRS, a row at 1e200 set the tolerance of the location check to
  # 1e191, so folds built on other data of the same size were accepted.  Read
  # first: it fails without waiting on PROJ.
  crsless <- function(seed) {
    d <- .ful_layer(seed = seed)
    sf::st_crs(d) <- NA
    g <- sf::st_geometry(d)
    g[1] <- sf::st_sfc(sf::st_point(c(5e5, 1e200)))
    sf::st_set_geometry(d, g)
  }
  f_other <- suppressWarnings(make_folds(crsless(78L), k = 5, method = "random_kfold", seed = 1))
  expect_error(suppressWarnings(fold_separation(f_other, crsless(77L))),
               "built from different data")
  f_same <- suppressWarnings(make_folds(crsless(77L), k = 5, method = "random_kfold", seed = 1))
  expect_s3_class(suppressWarnings(fold_separation(f_same, crsless(77L))), "fold_separation")
  # In EPSG:3857 an x of 1e200 never came back from PROJ's move to EPSG:4326,
  # which the probe of the layer as supplied takes: make_folds() did not
  # return.  The probe is first taken with a move that refuses such a row, so
  # a probe that hands it to PROJ fails here and does not hang.
  far <- .ful_put(L, 3857, c(1e200, 5e6))
  far$..row_id <- seq_len(nrow(far))
  real_tf <- sf::st_transform
  pr <- local({
    local_mocked_bindings(st_transform = function(x, crs, ...) {
      xy <- sf::st_coordinates(x)
      if (any(is.finite(xy[, 1:2]) & abs(xy[, 1:2]) > 1e150))
        stop("the probe moved a row beyond 1e150")
      real_tf(x, crs, ...)
    }, .package = "sf")
    .fold_row_probe(far, catch = FALSE)
  })
  expect_identical(pr$row_id[1L], 1L)
  expect_true(is.na(pr$x[1L]) && is.na(pr$y[1L]))
  ctl <- make_folds(.ful_put(L, 3857, NULL), k = 5, method = "block_kfold", seed = 1)
  r <- make_folds(.ful_put(L, 3857, c(1e200, 5e6)), k = 5, method = "block_kfold", seed = 1)
  expect_identical(r, ctl)
})

test_that("determine_optimal_levels() drops a row beyond 1e150 as given, as the same row emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(150L)
  for (epsg in c(32617, 32616, 4326)) {
    ctl <- .ful_run(determine_optimal_levels(.ful_put(L, epsg, NULL), max_levels = 5,
                                             set_seed = 1))
    ats <- if (epsg == 4326) list(c(-80, 1e200), c(-80, 95), c(Inf, 45))
           else c(list(c(Inf, 5e6)), .ful_far)
    for (at in ats) {
      lab <- paste(epsg, paste(at, collapse = ", "))
      r <- .ful_run(withCallingHandlers(
        determine_optimal_levels(.ful_put(L, epsg, at), max_levels = 5, set_seed = 1),
        warning = function(w)
          if (grepl("bounding box has potentially an invalid value range",
                    conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")))
      expect_identical(r$value, ctl$value, label = lab)
      expect_identical(r$msgs, .ful_msgs_at(ctl$msgs, at, 150L), label = lab)
    }
    expect_true(any(grepl(paste0("determine_optimal_levels(): dropping 1 point(s) ",
                                 .ful_unloc, "."), ctl$msgs, fixed = TRUE)))
  }
  ctl <- .ful_run(determine_optimal_levels(.ful_squares(L, empty = TRUE), max_levels = 5,
                                           set_seed = 1))
  for (bad in list(c(-80, 1e200), c(Inf, 45))) {
    r <- .ful_run(determine_optimal_levels(.ful_squares(L, bad = bad), max_levels = 5,
                                           set_seed = 1))
    expect_identical(r$value, ctl$value, label = paste(bad, collapse = ", "))
    expect_identical(r$msgs, ctl$msgs, label = paste(bad, collapse = ", "))
  }
})

test_that("cv_block_size_sweep() runs on a layer with a row at Inf or beyond 1e150, as on the row emptied", {
  skip_if_not_installed("ranger")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(150L)
  rf_fn <- function(tr) fit_rf_model(tr, "z", "wq", num_trees = 30, importance = "none")
  sweep <- function(d) suppressWarnings(cv_block_size_sweep(
    d, "z", "wq", fit_fn = rf_fn, k = 3, n_sizes = 2, include_random = FALSE,
    quiet = TRUE, seed = 1))
  ctl <- sweep(.ful_put(L, 32617, NULL))
  for (at in list(c(Inf, 5e6), c(5e5, 1e200)))
    expect_equal(sweep(.ful_put(L, 32617, at)), ctl, label = paste(at, collapse = ", "))
})

# --- BUG-140, BUG-141, BUG-142, BUG-147, BUG-181: seeds, cells, joins, moves --
# The seeders clustered a sample point at a northing of 1e200, or the place a
# move into the boundary's CRS gave it; the Voronoi cells and the clip target
# stopped with GEOS's message on a seed at Inf or 1e200, or built an infinite
# target; assign_features_to_polygons() joined such a feature where the move
# put it and blamed the polygons' cover for one at Inf; the variogram design
# effect and kriging_adequacy() correlated and kriged it; and the readers
# that move a layer (harmonize_crs(), ensure_projected(target_crs =) on a
# projected layer, and the internal moves of a boundary or of cells)
# returned a polygon with such a vertex as a smaller finite shape.  Such rows
# are now read as given, before any move, and read as the same rows empty.

.ful_box <- function(crs = 32617)
  sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(
    c(xmin = 5e5, ymin = 5e6, xmax = 6e5, ymax = 5.1e6), crs = crs)))

.ful_bad <- list(c(Inf, 5e6), c(5e5, 1e200), c(1e200, 5e6), c(-1e151, 5e6))

# Each variant of row `i` of `layer` in `epsg` against the same row emptied:
# value, warning messages and classes identical.
.ful_same_as_emptied <- function(layer, epsg, fn, ats = .ful_bad, i = 1L) {
  ctl <- .ful_run(fn(.ful_put(layer, epsg, NULL, i)))
  for (at in ats) {
    r <- .ful_run(fn(.ful_put(layer, epsg, at, i)))
    lab <- paste0("EPSG:", epsg, " at (", paste(at, collapse = ", "), ")")
    expect_identical(r$value, ctl$value, label = lab)
    expect_identical(r$msgs, ctl$msgs, label = lab)
    expect_identical(r$classes, ctl$classes, label = lab)
  }
  ctl
}

test_that("voronoi_seeds_kmeans() drops a point at Inf or beyond 1e150 as given, as the same point emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  ctl <- .ful_same_as_emptied(L, 32617, function(d) voronoi_seeds_kmeans(d, k = 8))
  expect_identical(nrow(ctl$value), 8L)
  expect_identical(ctl$msgs, paste0("voronoi_seeds_kmeans(): dropping 1 point(s) ",
                                    .ful_unloc, "."))
  # In lon/lat, a latitude of 1e200 takes no part in the projection either.
  .ful_same_as_emptied(L, 4326, function(d) voronoi_seeds_kmeans(d, k = 8),
                       ats = list(c(-80, 1e200), c(Inf, 45)))
})

test_that("get_voronoi_seeds(kmeans) reads a sample point beyond 1e150 as given, before the move into the boundary's CRS", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  seeds <- function(d) get_voronoi_seeds(boundary = .ful_box(), method = "kmeans", n = 8,
                                         sample_points = d, set_seed = 9)
  for (epsg in c(32617, 32616)) {
    ctl <- .ful_same_as_emptied(L, epsg, seeds)
    expect_identical(nrow(ctl$value), 8L)
    expect_identical(sf::st_crs(ctl$value), sf::st_crs(32617))
    expect_identical(attr(ctl$value, "kmeans")$rows, 2:120)
    expect_identical(ctl$msgs, paste0("get_voronoi_seeds(kmeans): dropping 1 point(s) ",
                                      .ful_unloc, "."))
  }
})

test_that("the record of scored seeds leaves a cloud point beyond 1e150 out, as the same point emptied", {
  L <- .ful_layer()
  seeds <- sf::st_as_sf(data.frame(x = c(5.2e5, 5.8e5), y = c(5.05e6, 5.05e6)),
                        coords = c("x", "y"), crs = 32617)
  for (epsg in c(32617, 32616)) {
    ctl <- .scored_seeds_record(.ful_put(L, epsg, NULL), seeds, 3L)
    expect_identical(ctl$rows, 2:120)
    expect_true(is.finite(ctl$tot_withinss))
    for (at in .ful_bad)
      expect_identical(.scored_seeds_record(.ful_put(L, epsg, at), seeds, 3L), ctl,
                       label = paste(epsg, paste(at, collapse = ", ")))
  }
})

test_that("create_voronoi_polygons() and build_tessellation() index a seed beyond 1e150 NA, as the same seed emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  S <- suppressWarnings(voronoi_seeds_kmeans(L, k = 12))
  cv <- function(d) create_voronoi_polygons(d, boundary = .ful_box(), quiet = TRUE)
  for (epsg in c(32617, 32616)) {
    ctl <- .ful_same_as_emptied(S, epsg, cv)
    expect_identical(nrow(ctl$value$cells), 11L)
    expect_true(is.na(ctl$value$index[1L]))
    expect_identical(ctl$classes, "spatialkit_unassigned")
    expect_match(ctl$msgs, paste0("They have no location (an empty point, or a coordinate ",
                                  "that is not finite or lies beyond 1e150), so there is ",
                                  "nothing to index."), fixed = TRUE)
  }
  # No boundary: the cells are clipped to the hull of the other seeds.
  ctl <- .ful_same_as_emptied(S, 32617, function(d)
    build_tessellation(d, method = "voronoi", quiet = TRUE))
  expect_identical(nrow(ctl$value$cells), 11L)
})

test_that("build_tessellation() counts an observation beyond 1e150 with the empty ones, under every method", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  for (m in c("hex", "square", "voronoi")) {
    bt <- function(d) build_tessellation(d, boundary = .ful_box(), method = m,
                                         approx_n_cells = if (m != "voronoi") 20,
                                         quiet = TRUE)
    ctl <- .ful_same_as_emptied(L, 32617, bt, ats = .ful_bad[1:2])
    expect_true(is.na(ctl$value$index[1L]), label = m)
    expect_match(ctl$msgs, "They have no location (an empty point", fixed = TRUE,
                 label = m)
  }
})

test_that("clip_target_for() leaves a point beyond 1e150 out of the box, as the same point emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  for (epsg in c(32617, 32616)) {
    ctl <- .ful_same_as_emptied(L, epsg, function(d) clip_target_for(d, expand = 0.05,
                                                                     quiet = TRUE))
    expect_true(is.finite(sum(as.numeric(sf::st_area(ctl$value)))))
  }
})

test_that("assign_features_to_polygons() matches a feature beyond 1e150 to no polygon, whatever the move", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  cells <- suppressWarnings(build_tessellation(L, boundary = .ful_box(), method = "square",
                                               approx_n_cells = 16, quiet = TRUE))$cells
  asg <- function(d) assign_features_to_polygons(d, cells, polygon_id_col = "cell_id",
                                                 quiet = TRUE)
  ctl <- .ful_same_as_emptied(L, 32617, asg)
  expect_identical(nrow(ctl$value), 119L)
  expect_identical(attr(ctl$value, "unassigned")$which, 1L)
  expect_match(ctl$msgs, paste0("They have no location (an empty point, or a coordinate ",
                                "that is not finite or lies beyond 1e150), so there is ",
                                "nothing to assign"), fixed = TRUE)
  expect_no_match(ctl$msgs, "cover the features", fixed = TRUE)
  # In EPSG:32616 beside cells in EPSG:32617: the move put a northing of 1e200
  # at (87304, 5766268), and a polygon there took it.
  big <- sf::st_sf(cell_id = 1:2, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(4e5, 4.9e6), c(7e5, 4.9e6), c(7e5, 5.2e6), c(4e5, 5.2e6),
                              c(4e5, 4.9e6)))),
    sf::st_polygon(list(rbind(c(0, 5.6e6), c(2e5, 5.6e6), c(2e5, 5.9e6), c(0, 5.9e6),
                              c(0, 5.6e6)))), crs = 32617))
  # The rows come back with the geometry they were given.
  ctl <- .ful_same_as_emptied(L, 32616, function(d) {
    a <- assign_features_to_polygons(d, big, polygon_id_col = "cell_id",
                                     keep_unassigned = TRUE, quiet = TRUE)
    list(a = sf::st_drop_geometry(a), unassigned = attr(a, "unassigned"))
  })
  expect_true(is.na(ctl$value$a$cell_id[1L]))
  expect_false(any(ctl$value$a$cell_id %in% 2L))
  expect_identical(ctl$value$unassigned$which, 1L)
  # A feature without a location beside one outside the polygons: counted
  # apart, and the coverage advice kept for the other.
  d <- .ful_put(.ful_put(L, 32617, c(5e5, 1e200)), 32617, c(9e5, 9e6), i = 2L)
  r <- .ful_run(asg(d))
  expect_match(r$msgs, paste0("1 of them has no location (an empty point, or a coordinate ",
                              "that is not finite or lies beyond 1e150); for the rest, ",
                              "check that the polygons cover"), fixed = TRUE)
})

test_that("assign_features_to_polygons() matches a polygon feature with a vertex beyond 1e150 to no polygon", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(40L)
  cells <- suppressWarnings(build_tessellation(sf::st_transform(L, 4326),
                                               boundary = sf::st_transform(.ful_box(), 4326),
                                               method = "square", approx_n_cells = 9,
                                               quiet = TRUE))$cells
  ctl <- .ful_run(assign_features_to_polygons(.ful_squares(L, empty = TRUE), cells,
                                              polygon_id_col = "cell_id", quiet = TRUE))
  expect_identical(attr(ctl$value, "unassigned")$which, 7L)
  for (bad in list(c(1e200, 45), c(-80, 1e200), c(Inf, 45))) {
    r <- .ful_run(assign_features_to_polygons(.ful_squares(L, bad = bad), cells,
                                              polygon_id_col = "cell_id", quiet = TRUE))
    lab <- paste(bad, collapse = ", ")
    expect_identical(sf::st_drop_geometry(r$value), sf::st_drop_geometry(ctl$value), label = lab)
    expect_identical(r$msgs, ctl$msgs, label = lab)
  }
})

test_that("summarize_by_cell(deff = 'variogram') leaves a point beyond 1e150 out of the correlation, as the same point emptied", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(200L)
  sac <- suppressWarnings(estimate_sac_range(L, "z", seed = 1))
  tess <- suppressWarnings(build_tessellation(L, boundary = .ful_box(), method = "square",
                                              approx_n_cells = 16, quiet = TRUE))
  A <- assign_features_to_polygons(L, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  sm <- function(d) summarize_by_cell(d, "z", id_col = "cell_id", cells_sf = tess$cells,
                                      deff = "variogram", sac = sac, estimate = "mean")
  for (epsg in c(32617, 32616)) {
    ctl <- .ful_same_as_emptied(A, epsg, sm)
    expect_match(ctl$msgs, paste0("1 point(s) have no location (an empty point, or a ",
                                  "coordinate that is not finite or lies beyond 1e150); ",
                                  "their values are summarised"), fixed = TRUE, all = FALSE)
  }
  # The correlation statistics read the bound too.
  cf <- function(h) exp(-h / 3e4)
  xy <- sf::st_coordinates(A)[1:30, 1:2]
  ref <- .cor_stats_from_coords(xy[-1, ], cf)
  for (at in .ful_bad) {
    xy1 <- xy; xy1[1, ] <- at
    expect_identical(.cor_stats_from_coords(xy1, cf), ref, label = paste(at, collapse = ", "))
  }
})

test_that("kriging_adequacy() leaves a point beyond 1e150 out before gstat and the folds, as the same point emptied", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(200L)
  sac <- suppressWarnings(estimate_sac_range(L, "z", seed = 1))
  tess <- suppressWarnings(build_tessellation(L, boundary = .ful_box(), method = "square",
                                              approx_n_cells = 16, quiet = TRUE))
  A <- assign_features_to_polygons(L, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  ka <- function(d) kriging_adequacy(d, "z", tess$cells, id_col = "cell_id", sac = sac,
                                     k = 4, seed = 1)
  for (epsg in c(32617, 32616)) {
    ctl <- .ful_same_as_emptied(A, epsg, ka)
    expect_identical(nrow(ctl$value), nrow(tess$cells))
  }
})

# A half-degree square near (-93.5, 45.5) and, in row 2, another with vertex
# 3 at `v` (an empty polygon for NULL), as polygons or lines.
.ful_pair <- function(v, crs = 4326, line = FALSE, x0 = c(-93.2, -93.9),
                      y0 = c(45.1, 45.3), d = 0.5) {
  sq <- function(x, y) rbind(c(x, y), c(x + d, y), c(x + d, y + d), c(x, y + d), c(x, y))
  b <- sq(x0[2], y0[2]); if (!is.null(v)) b[3L, ] <- v
  g2 <- if (line) { if (is.null(v)) sf::st_linestring() else sf::st_linestring(b[1:4, ]) }
        else if (is.null(v)) sf::st_polygon() else sf::st_polygon(list(b))
  g1 <- if (line) sf::st_linestring(sq(x0[1], y0[1])[1:4, ])
        else sf::st_polygon(list(sq(x0[1], y0[1])))
  sf::st_sf(v = 1:2, geometry = sf::st_sfc(g1, g2, crs = crs))
}

test_that("the readers that move a layer return a polygon or line with a vertex beyond 1e150 empty", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  tgt <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(5e5, 5e6)), crs = 32615))
  readers <- list(
    harmonize_b  = function(a) harmonize_crs(a, tgt, prefer = "b")$a,
    harmonize_t  = function(a) harmonize_crs(a, tgt, target_crs = 32615)$a,
    harmonize_a  = function(a) harmonize_crs(tgt, a, prefer = "a")$b,
    stamp        = function(a) .transform_or_stamp(a, sf::st_crs(32615), "x", "p"),
    align        = function(a) .align_crs(a, tgt),
    ensure       = function(a) ensure_projected(a, target_crs = 32615))
  cases <- list(
    list(crs = 4326, bad = list(c(1e200, 45.8), c(-93.4, 1e200), c(-1e200, 45.8),
                                c(Inf, 45.8))),
    list(crs = NA, bad = list(c(1e200, 45.8), c(-93.4, 1e200), c(Inf, 45.8))),
    list(crs = 32616, x0 = c(4e5, 3e5), y0 = c(5e6, 5e6), d = 4e4,
         bad = list(c(1e200, 5.05e6), c(3.2e5, 1e200), c(Inf, 5.05e6))))
  for (cs in cases) for (line in c(FALSE, TRUE)) {
    mk <- function(v) {
      a <- if (is.null(cs$x0)) .ful_pair(v, line = line)
           else .ful_pair(v, crs = cs$crs, line = line, x0 = cs$x0, y0 = cs$y0, d = cs$d)
      if (is.na(cs$crs)) sf::st_set_crs(a, NA) else a
    }
    for (nm in names(readers)) {
      if (is.na(cs$crs) && nm %in% c("align", "ensure")) next
      ctl <- .ful_run(readers[[nm]](mk(NULL)))
      for (v in cs$bad) {
        lab <- sprintf("%s, CRS %s, %s, vertex (%s)", nm, cs$crs,
                       if (line) "line" else "polygon", paste(v, collapse = ", "))
        r <- .ful_run(readers[[nm]](mk(v)))
        expect_identical(r$value, ctl$value, label = lab)
        # Without a CRS, the lon/lat test reads the other vertices of the
        # row (BUG-31), so the box its warning quotes holds them.
        expect_identical(gsub("\\(xmin=[^)]*\\)", "(box)", r$msgs),
                         gsub("\\(xmin=[^)]*\\)", "(box)", ctl$msgs), label = lab)
        expect_true(sf::st_is_empty(sf::st_geometry(r$value)[[2]]), label = lab)
        expect_s3_class(sf::st_geometry(r$value)[[2]],
                        if (line) "LINESTRING" else "POLYGON")
      }
    }
  }
  # A point row comes back an empty point.
  p <- sf::st_sf(v = 1:2, geometry = sf::st_sfc(sf::st_point(c(4e5, 5e6)),
                                                 sf::st_point(c(5e5, 1e200)), crs = 32616))
  for (nm in c("harmonize_b", "stamp", "align", "ensure")) {
    r <- readers[[nm]](p)
    expect_true(sf::st_is_empty(sf::st_geometry(r)[[2]]), label = nm)
    expect_s3_class(sf::st_geometry(r)[[2]], "POINT")
  }
  # A layer whose rows all have a location is moved as sf::st_transform()
  # moves it.
  a <- .ful_pair(c(-93.4, 45.8))
  for (nm in setdiff(names(readers), "harmonize_a"))
    expect_identical(readers[[nm]](a), sf::st_transform(a, 32615), label = nm)
  b <- .ful_pair(c(3.2e5, 5.04e6), crs = 32616, x0 = c(4e5, 3e5), y0 = c(5e6, 5e6), d = 4e4)
  expect_identical(ensure_projected(b, target_crs = 32615), sf::st_transform(b, 32615))
  expect_identical(harmonize_crs(b, tgt, prefer = "b")$a, sf::st_transform(b, 32615))
})

test_that("a polygon layer with a vertex beyond 1e150 joins and grids as the same row emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  pts <- sf::st_transform(sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_point(c(-93.85, 45.35)), sf::st_point(c(-93.45, 45.75)),
    sf::st_point(c(-93.0, 45.3)), crs = 4326)), 32615)
  P0 <- .ful_pair(NULL); P0$poly_id <- c("A", "B")
  ctl <- .ful_run(assign_features_to_polygons(pts, P0, keep_unassigned = TRUE, quiet = TRUE))
  expect_identical(ctl$value$poly_id, c(NA, NA, "A"))
  g0 <- tryCatch(create_grid_polygons(P0[2, ], cellsize = 5000, crs = 32615, quiet = TRUE),
                 error = conditionMessage)
  for (v in list(c(1e200, 45.8), c(-93.4, 1e200), c(Inf, 45.8))) {
    P <- .ful_pair(v); P$poly_id <- c("A", "B")
    lab <- paste(v, collapse = ", ")
    r <- .ful_run(assign_features_to_polygons(pts, P, keep_unassigned = TRUE, quiet = TRUE))
    expect_identical(r$value$poly_id, ctl$value$poly_id, label = lab)
    expect_identical(r$msgs, ctl$msgs, label = lab)
    g <- tryCatch(create_grid_polygons(P[2, ], cellsize = 5000, crs = 32615, quiet = TRUE),
                  error = conditionMessage)
    expect_identical(g, g0, label = lab)
  }
})

# --- BUG-141, BUG-181: the builders' cover ----------------------------------
# A point of a point `cover` at Inf or 1e200 set the whole cover aside
# ("gives no polygon"), and a polygon cover with such a vertex was moved into
# the working CRS as a smaller shape, or one far larger.

test_that("the builders read a point of `cover` beyond 1e150 as given, as the same point emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  S <- suppressWarnings(voronoi_seeds_kmeans(L, k = 10))
  attr(S, "cover") <- NULL
  carried <- function(d) { s <- S; attr(s, "cover") <- d; s }
  for (epsg in c(32617, 32616)) {
    ctl <- .ful_same_as_emptied(L, epsg, function(d)
      build_tessellation(S, method = "voronoi", cover = d, quiet = TRUE))
    expect_identical(ctl$value$params$boundary_source, "cover")
    expect_identical(ctl$msgs, character(0))
    ctl <- .ful_same_as_emptied(L, epsg, function(d)
      create_voronoi_polygons(S, cover = d, quiet = TRUE))
    expect_identical(ctl$value$params$boundary_source, "cover")
    # The cover the seeds carry.
    ctl <- .ful_same_as_emptied(L, epsg, function(d)
      create_voronoi_polygons(carried(d), quiet = TRUE))
    expect_identical(ctl$value$params$boundary_source, "seeds' cover")
    ctl <- .ful_same_as_emptied(L, epsg, function(d)
      build_tessellation(carried(d), method = "voronoi", quiet = TRUE))
    expect_identical(ctl$value$params$boundary_source, "seeds' cover")
    # With a boundary the cover is only checked against the cells.
    for (m in c("voronoi", "hex")) {
      ctl <- .ful_same_as_emptied(L, epsg, function(d)
        build_tessellation(S, boundary = .ful_box(), method = m, cover = d,
                           approx_n_cells = if (m != "voronoi") 20, quiet = TRUE))
      expect_identical(ctl$value$params$cover_outside, 0L, label = m)
    }
  }
})

test_that("a polygon `cover` with a vertex beyond 1e150 is moved into the working CRS as the same row emptied", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  S <- suppressWarnings(voronoi_seeds_kmeans(L, k = 10))
  attr(S, "cover") <- NULL
  # The box, in EPSG:32616, in two halves; vertex 3 of the second at `v`.
  b <- sf::st_coordinates(sf::st_transform(.ful_box(), 32616))[, 1:2]
  m1 <- (b[1, ] + b[2, ]) / 2; m2 <- (b[3, ] + b[4, ]) / 2
  h1 <- rbind(b[1, ], m1, m2, b[4, ], b[1, ])
  h2 <- rbind(m1, b[2, ], b[3, ], m2, m1)
  mk <- function(v) {
    g2 <- if (is.null(v)) sf::st_polygon() else { q <- h2; q[3, ] <- v; sf::st_polygon(list(q)) }
    sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(unname(h1))), g2, crs = 32616))
  }
  builders <- list(
    voronoi = function(cv) create_voronoi_polygons(S, cover = cv, quiet = TRUE),
    build   = function(cv) build_tessellation(S, method = "voronoi", cover = cv, quiet = TRUE),
    carried = function(cv) { s <- S; attr(s, "cover") <- sf::st_geometry(cv)
                             create_voronoi_polygons(s, quiet = TRUE) })
  for (nm in names(builders)) {
    ctl <- .ful_run(builders[[nm]](mk(NULL)))
    expect_identical(ctl$msgs, character(0), label = nm)
    for (v in list(c(3.2e5, 1e200), c(1e200, 5.05e6), c(Inf, 5.05e6))) {
      lab <- sprintf("%s, vertex (%s)", nm, paste(v, collapse = ", "))
      r <- .ful_run(builders[[nm]](mk(v)))
      expect_identical(r$value, ctl$value, label = lab)
      expect_identical(r$msgs, ctl$msgs, label = lab)
    }
  }
  # Moved alone, the cover comes back with the row empty, in its own type.
  g <- .cover_follow(sf::st_geometry(mk(c(3.2e5, 1e200))), S, S)
  expect_true(sf::st_is_empty(g[[2]]))
  expect_s3_class(g[[2]], "POLYGON")
  expect_identical(g[1], sf::st_transform(sf::st_geometry(mk(NULL))[1], 32617))
})

# --- BUG-142: polygon features in the variogram design effect ----------------

test_that("summarize_by_cell(deff = 'variogram') reads a polygon feature with a vertex beyond 1e150 as the same row emptied", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(200L)
  sac <- suppressWarnings(estimate_sac_range(L, "z", seed = 1))
  tess <- suppressWarnings(build_tessellation(L, boundary = .ful_box(), method = "square",
                                              approx_n_cells = 16, quiet = TRUE))
  xy <- sf::st_coordinates(L)
  sq <- function(x, y, d = 500)
    sf::st_polygon(list(rbind(c(x, y), c(x + d, y), c(x + d, y + d), c(x, y + d), c(x, y))))
  PL <- sf::st_sf(z = L$z, geometry = sf::st_sfc(
    lapply(seq_len(nrow(xy)), function(i) sq(xy[i, 1], xy[i, 2])), crs = 32617))
  A <- suppressWarnings(assign_features_to_polygons(PL, tess$cells, polygon_id_col = "cell_id",
                                                    quiet = TRUE))
  # Row 1 emptied (NULL), or vertex 3 put at `v` (an NA keeps that
  # coordinate).  Reduced to a point before any move, the square with a
  # vertex at (Inf, 5.05e6) gives a finite point inside it.
  bad1 <- function(A, v) {
    g <- sf::st_geometry(A)
    if (is.null(v)) {
      g[1] <- sf::st_sfc(sf::st_polygon(), crs = sf::st_crs(g))
    } else {
      b <- sf::st_coordinates(g[[1]])[, 1:2]
      b[3, !is.na(v)] <- v[!is.na(v)]
      g[1] <- sf::st_sfc(sf::st_polygon(list(unname(b))), crs = sf::st_crs(g))
    }
    sf::st_set_geometry(A, g)
  }
  sm <- function(d) summarize_by_cell(d, "z", id_col = "cell_id", cells_sf = tess$cells,
                                      deff = "variogram", sac = sac, estimate = "mean")
  for (epsg in c(32617, 32616)) {
    Ae <- sf::st_transform(A, epsg)
    ctl <- .ful_run(sm(bad1(Ae, NULL)))
    expect_match(ctl$msgs, "1 point(s) have no location", fixed = TRUE, all = FALSE)
    for (v in list(c(Inf, 5.05e6), c(Inf, NA), c(NA, 1e200), c(1e200, 5.05e6))) {
      lab <- sprintf("EPSG:%s, vertex 3 at (%s)", epsg, paste(v, collapse = ", "))
      r <- .ful_run(sm(bad1(Ae, v)))
      expect_identical(r$value, ctl$value, label = lab)
      expect_identical(r$msgs, ctl$msgs, label = lab)
    }
  }
})

# --- BUG-181: the other moves of a caller's polygons -------------------------

test_that("assign_features_to_polygons() moves lon/lat cells into the scored record's CRS with a row beyond 1e150 empty", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  pts <- sf::st_transform(sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_point(c(-93.85, 45.35)), sf::st_point(c(-93.45, 45.75)),
    sf::st_point(c(-93.0, 45.3)), crs = 4326)), 32615)
  # A record that says the cells are the ones a profile scored in EPSG:32615.
  cen <- sf::st_coordinates(sf::st_transform(sf::st_sfc(
    sf::st_point(c(-92.95, 45.35)), sf::st_point(c(-93.65, 45.55)), crs = 4326), 32615))[, 1:2]
  rec <- structure(list(version = 1L, rule = "size", estimate = "mean", cells = 2L,
                        scored = list(centres = unname(cen)), crs = sf::st_crs(32615),
                        scored_cells = TRUE),
                   class = "spatialkit_count_rule")
  cells <- function(v) {
    P <- .ful_pair(v); P$poly_id <- c("A", "B")
    attr(P, "count_rule") <- rec
    P
  }
  expect_true(.scored_crs(.count_rule_get(cells(NULL))) == sf::st_crs(32615))
  asg <- function(v) .ful_run(assign_features_to_polygons(pts, cells(v), keep_unassigned = TRUE,
                                                          quiet = TRUE))
  ctl <- asg(NULL)
  expect_identical(ctl$value$poly_id, c(NA, NA, "A"))
  for (v in list(c(1e200, 45.8), c(-93.4, 1e200), c(Inf, 45.8))) {
    lab <- paste(v, collapse = ", ")
    r <- asg(v)
    expect_identical(r$value$poly_id, ctl$value$poly_id, label = lab)
    expect_identical(r$msgs, ctl$msgs, label = lab)
  }
})

test_that("cell_values_at()'s free cells read a polygon with a vertex beyond 1e150 as the same row emptied", {
  P0 <- .ful_pair(NULL); P0$cell_id <- 1:2
  seed_xy <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(-93.85, 45.35)), crs = 4326), 32615))[, 1:2, drop = FALSE]
  free <- function(P) .values_at_free_cells(NA_integer_, seed_xy, P, "cell_id",
                                            sf::st_crs(32615), sf::NA_crs_)
  ctl <- free(P0)
  expect_identical(ctl$id, "1")
  for (v in list(c(1e200, 45.8), c(-93.4, 1e200), c(Inf, 45.8))) {
    P <- .ful_pair(v); P$cell_id <- 1:2
    expect_identical(free(P), ctl, label = paste(v, collapse = ", "))
  }
})

test_that("build_tessellation() on lon/lat points reads a latitude beyond 1e150 as given, under every method", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  bnd <- sf::st_transform(.ful_box(), 4326)
  for (m in c("hex", "square", "voronoi")) {
    ctl <- .ful_same_as_emptied(L, 4326, function(d)
      build_tessellation(d, boundary = bnd, method = m,
                         approx_n_cells = if (m != "voronoi") 20, quiet = TRUE),
      ats = list(c(-80.5, 1e200), c(Inf, 45.3), c(1e200, 45.3)))
    expect_true(is.na(ctl$value$index[1L]), label = m)
  }
})

# --- The paths a reader's emptying does not reach ----------------------------

test_that("get_voronoi_seeds(kmeans) without a boundary attaches the hull of the located sample points", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer()
  seeds <- function(d) get_voronoi_seeds(method = "kmeans", n = 8, sample_points = d,
                                         set_seed = 9)
  for (epsg in c(32617, 32616)) {
    ctl <- .ful_same_as_emptied(L, epsg, seeds)
    expect_true(all(abs(sf::st_bbox(attr(ctl$value, "cover"))) < 1e7))
  }
})

test_that("build_tessellation() says a point that has no place once projected has no location", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_put(.ful_layer(), 4326, c(-80.5, 95))
  r <- .ful_run(build_tessellation(L, boundary = sf::st_transform(.ful_box(), 4326),
                                   method = "hex", approx_n_cells = 20, quiet = TRUE))
  expect_true(is.na(r$value$index[1L]))
  expect_match(r$msgs, paste0("fall in no cell and are NA in `index`; `params$n_outside` ",
                              "holds the count. They have no location"),
               fixed = TRUE, all = FALSE)
})

test_that("assign_features_to_polygons() moves a feature beyond 1e150 into the scored record's CRS empty", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  # From EPSG:32616 a northing of 1e200 lands at a finite place in
  # EPSG:32615; the lon/lat cells are drawn round that place.
  bad <- c(3.2e5, 1e200)
  land <- sf::st_transform(sf::st_sfc(sf::st_point(bad), crs = 32616), 32615)
  expect_true(all(is.finite(sf::st_coordinates(land))))
  ll <- sf::st_coordinates(sf::st_transform(land, 4326))[1, 1:2]
  sq <- function(x, y, d = 0.2)
    sf::st_polygon(list(rbind(c(x, y), c(x + d, y), c(x + d, y + d), c(x, y + d), c(x, y))))
  cells <- sf::st_sf(poly_id = c("A", "B"), geometry = sf::st_sfc(
    sq(ll[1] - 0.1, ll[2] - 0.1), sq(ll[1] + 0.1, ll[2] - 0.1), crs = 4326))
  cen <- sf::st_coordinates(sf::st_transform(sf::st_centroid(sf::st_geometry(cells)),
                                             32615))[, 1:2]
  attr(cells, "count_rule") <- structure(
    list(version = 1L, rule = "size", estimate = "mean", cells = 2L,
         scored = list(centres = unname(cen)), crs = sf::st_crs(32615), scored_cells = TRUE),
    class = "spatialkit_count_rule")
  pts <- sf::st_transform(sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_point(ll), sf::st_point(ll + c(0.05, 0)), sf::st_point(ll + c(0.15, 0)),
    crs = 4326)), 32616)
  asg <- function(at) .ful_run(assign_features_to_polygons(
    .ful_put(pts, 32616, at), cells, keep_unassigned = TRUE, quiet = TRUE))
  ctl <- asg(NULL)
  expect_identical(ctl$value$poly_id, c(NA, "A", "B"))
  r <- asg(bad)
  expect_identical(sf::st_drop_geometry(r$value), sf::st_drop_geometry(ctl$value))
  expect_identical(r$msgs, ctl$msgs)
})

test_that("kriging_adequacy() reads a multipoint with a member beyond 1e150 as the same row emptied", {
  skip_if_not_installed("gstat")
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  L <- .ful_layer(200L)
  sac <- suppressWarnings(estimate_sac_range(L, "z", seed = 1))
  tess <- suppressWarnings(build_tessellation(L, boundary = .ful_box(), method = "square",
                                              approx_n_cells = 16, quiet = TRUE))
  A <- assign_features_to_polygons(L, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  M <- sf::st_cast(A, "MULTIPOINT")
  row1 <- function(g1) { g <- sf::st_geometry(M); g[1] <- sf::st_sfc(g1, crs = 32617)
                         sf::st_set_geometry(M, g) }
  ka <- function(d) kriging_adequacy(d, "z", tess$cells, id_col = "cell_id", sac = sac,
                                     k = 4, seed = 1)
  ctl <- .ful_run(ka(row1(sf::st_multipoint())))
  p1 <- sf::st_coordinates(sf::st_geometry(A)[[1]])[1:2]
  for (m in list(c(p1[1], 1e200), c(Inf, p1[2]), c(1e200, p1[2]))) {
    lab <- paste(m, collapse = ", ")
    r <- .ful_run(ka(row1(sf::st_multipoint(rbind(p1, m)))))
    expect_identical(r$value, ctl$value, label = lab)
    expect_identical(r$msgs, ctl$msgs, label = lab)
  }
})
