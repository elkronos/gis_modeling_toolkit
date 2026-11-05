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
