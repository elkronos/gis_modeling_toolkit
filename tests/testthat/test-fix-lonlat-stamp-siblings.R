# BUG-272: harmonize_crs(), ensure_projected(target_crs =) and
# assign_features_to_polygons() stamped a CRS in longitude and latitude on
# a layer without a CRS whose coordinates lie outside longitude -180 to 360
# and latitude -90 to 90, which .transform_or_stamp() refuses since BUG-245.
# They now refuse it with the same words, naming their own argument.

.ls_pts <- function(x0, y0, w, crs = NA, n = 10L) {
  # Ten points on a fixed lattice in a square of side w at (x0, y0): the
  # box is known exactly, and nothing depends on the random stream.
  i <- seq_len(n) - 1L
  sf::st_as_sf(data.frame(x = x0 + w * (i %% 5L) / 4, y = y0 + w * (i %/% 5L)),
               coords = c("x", "y"), crs = crs)
}
.ls_sq <- function(x0, y0, w, crs = NA) {
  sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w), c(x0, y0 + w), c(x0, y0)))),
    crs = sf::st_crs(crs)))
}
# Every condition the call raises, and its value or its error.
.ls_catch <- function(expr) {
  ws <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  list(value = v, warnings = ws)
}
.ls_refused <- function(r, caller, what, remedy = FALSE) {
  expect_s3_class(r$value, "error")
  msg <- if (inherits(r$value, "error")) conditionMessage(r$value) else "<returned>"
  expect_match(msg, paste0(
    "^", caller, "\\(\\): `", what, "` has no CRS and its coordinates lie outside ",
    "longitude -180 to 360 and latitude -90 to 90 \\(xmin=[^)]+\\), so they cannot ",
    "be in the target CRS \\('EPSG:4326'\\), which is in longitude and latitude\\. ",
    "Set the CRS of `", what, "` with sf::st_crs\\(\\)",
    if (remedy) ", or pass a projected `target_crs`\\.$" else "\\.$"))
  # Neither the stamp's warning nor sf's bare bounding-box warning.
  expect_length(r$warnings, 0L)
  invisible(msg)
}

metres <- .ls_pts(5e5, 4e6, 2000)          # x 5e5 to 502000, y 4e6 to 4002000
y_over <- .ls_pts(10, 100, 20)             # a y beyond 90 alone
ll     <- .ls_pts(-80, 36, 0.5, crs = 4326)
utm    <- .ls_pts(5e5, 4e6, 2000, crs = 32617)

test_that("harmonize_crs() refuses a CRS in degrees on CRS-less coordinates outside the envelope", {
  for (m in list(metres, y_over)) {
    .ls_refused(.ls_catch(harmonize_crs(m, ll)), "harmonize_crs", "a")
    .ls_refused(.ls_catch(harmonize_crs(ll, m)), "harmonize_crs", "b")
    .ls_refused(.ls_catch(harmonize_crs(m, ll, prefer = "b")), "harmonize_crs", "a")
    .ls_refused(.ls_catch(harmonize_crs(m, utm, target_crs = 4326)),
                "harmonize_crs", "a", remedy = TRUE)
    .ls_refused(.ls_catch(harmonize_crs(utm, m, target_crs = 4326)),
                "harmonize_crs", "b", remedy = TRUE)
  }
  # The box is the one of the coordinates, as .transform_or_stamp() gives it.
  msg <- .ls_refused(.ls_catch(harmonize_crs(metres, ll)), "harmonize_crs", "a")
  expect_match(msg, "(xmin=500000, xmax=502000, ymin=4e+06, ymax=4.002e+06)", fixed = TRUE)
})

test_that("ensure_projected(target_crs =) refuses a CRS in degrees on CRS-less coordinates outside the envelope", {
  for (m in list(metres, y_over))
    .ls_refused(.ls_catch(ensure_projected(m, target_crs = 4326)),
                "ensure_projected", "x", remedy = TRUE)
})

test_that("assign_features_to_polygons() refuses either layer without a CRS beside one in degrees", {
  # Features in metres beside a cell in degrees: they fell in no polygon,
  # with sf's bounding-box warning three times.
  .ls_refused(.ls_catch(assign_features_to_polygons(metres, .ls_sq(-81, 35, 2, crs = 4326))),
              "assign_features_to_polygons", "features_sf")
  # Polygons in metres beside points in degrees: the call stopped on s2's
  # refusal, blamed on reading `polygons_sf` as longitude and latitude.
  r <- .ls_catch(assign_features_to_polygons(ll, .ls_sq(4.9e5, 3.9e6, 2e5)))
  msg <- .ls_refused(r, "assign_features_to_polygons", "polygons_sf")
  expect_no_match(msg, "s2", fixed = TRUE)
  .ls_refused(.ls_catch(assign_features_to_polygons(ll, .ls_sq(10, 100, 20))),
              "assign_features_to_polygons", "polygons_sf")
})

test_that("each bound alone refuses the stamp", {
  # x below -180, x above 360, y below -90, y above 90; the other
  # coordinate inside the envelope.
  for (m in list(.ls_pts(-5e5, 10, 1), .ls_pts(400, 10, 1),
                 .ls_pts(10, -4e6, 1), .ls_pts(10, 95, 1))) {
    .ls_refused(.ls_catch(harmonize_crs(m, ll)), "harmonize_crs", "a")
    .ls_refused(.ls_catch(ensure_projected(m, target_crs = 4326)),
                "ensure_projected", "x", remedy = TRUE)
  }
})

# Five points: the four corners of the box and its centre.
.ls_box <- function(xmin, xmax, ymin, ymax) {
  sf::st_as_sf(data.frame(x = c(xmin, xmax, xmin, xmax, (xmin + xmax) / 2),
                          y = c(ymin, ymin, ymax, ymax, (ymin + ymax) / 2)),
               coords = c("x", "y"))
}

test_that("each bound refuses half a unit past it and stamps at it", {
  # Half a degree past each bound, the other coordinate inside: refused.
  past <- list(.ls_box(-180.5, -179.5, 10, 11), .ls_box(359.5, 360.5, 10, 11),
               .ls_box(10, 11, -90.5, -89.5), .ls_box(10, 11, 89.5, 90.5))
  for (m in past) {
    msg <- .ls_refused(.ls_catch(harmonize_crs(m, ll)), "harmonize_crs", "a")
    .ls_refused(.ls_catch(ensure_projected(m, target_crs = 4326)),
                "ensure_projected", "x", remedy = TRUE)
  }
  expect_match(msg, "(xmin=10, xmax=11, ymin=89.5, ymax=90.5)", fixed = TRUE)
  # At each bound: inside the envelope, so stamped as before.  The x beyond
  # 180 (or the width of 480) keeps .looks_like_lonlat() from reading the
  # box as degrees, so each reaches the stamp.
  at <- list(.ls_box(-180, 300, 10, 11), .ls_box(200, 360, 10, 11),
             .ls_box(200, 201, -90, -89), .ls_box(200, 201, 89, 90))
  for (m in at) {
    r <- .ls_catch(harmonize_crs(m, ll))
    expect_false(inherits(r$value, "error"))
    expect_length(r$warnings, 1L)
    expect_match(r$warnings, "stamping 'EPSG:4326' onto it WITHOUT reprojection",
                 fixed = TRUE)
    expect_identical(sf::st_crs(r$value$a), sf::st_crs(4326))
    expect_equal(sf::st_coordinates(r$value$a), sf::st_coordinates(m))
    lines <- capture_spatialkit_log(r <- .ls_catch(ensure_projected(m, target_crs = 4326)))
    expect_false(inherits(r$value, "error"))
    expect_true(log_has(lines, "stamping the supplied target CRS \\('EPSG:4326'\\) WITHOUT reprojection"))
    expect_identical(sf::st_crs(r$value), sf::st_crs(4326))
  }
})

test_that("a projected target, degrees, a box inside the envelope and no location are as before", {
  # A projected target is stamped with the warning (harmonize_crs()) or the
  # log line (ensure_projected()), as before.
  r <- .ls_catch(harmonize_crs(metres, utm))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "stamping 'EPSG:32617' onto it WITHOUT reprojection", fixed = TRUE)
  expect_identical(sf::st_crs(r$value$a), sf::st_crs(32617))
  expect_equal(sf::st_coordinates(r$value$a), sf::st_coordinates(metres))
  lines <- capture_spatialkit_log(r <- .ls_catch(ensure_projected(metres, target_crs = 32617)))
  expect_length(r$warnings, 0L)
  expect_true(log_has(lines, "stamping the supplied target CRS \\('EPSG:32617'\\) WITHOUT reprojection"))
  expect_identical(sf::st_crs(r$value), sf::st_crs(32617))
  r <- .ls_catch(assign_features_to_polygons(metres, .ls_sq(4.9e5, 3.9e6, 2e5, crs = 32617),
                                             quiet = TRUE))
  expect_false(inherits(r$value, "error"))
  expect_identical(r$value$poly_id, rep(1L, 10L))
  # CRS-less degrees are read as EPSG:4326 and moved, as before.
  deg <- ll; sf::st_crs(deg) <- NA
  r <- .ls_catch(harmonize_crs(deg, utm))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "look like lon/lat", fixed = TRUE)
  expect_equal(sf::st_coordinates(r$value$a),
               sf::st_coordinates(sf::st_transform(ll, 32617)), tolerance = 1e-9)
  # x 200 to 250, y 10 to 60 lies inside the envelope: stamped with the
  # warning, as .transform_or_stamp() does.
  inside <- .ls_pts(200, 10, 50)
  r <- .ls_catch(harmonize_crs(inside, ll))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "stamping 'EPSG:4326' onto it WITHOUT reprojection", fixed = TRUE)
  expect_identical(sf::st_crs(r$value$a), sf::st_crs(4326))
  lines <- capture_spatialkit_log(r <- .ls_catch(ensure_projected(inside, target_crs = 4326)))
  expect_identical(sf::st_crs(r$value), sf::st_crs(4326))
  # A layer with no located coordinate (every coordinate Inf) has no box
  # to test and is stamped as before.
  inf <- sf::st_as_sf(data.frame(x = c(Inf, Inf), y = c(Inf, Inf)),
                      coords = c("x", "y"), na.fail = FALSE)
  r <- .ls_catch(harmonize_crs(inf, ll))
  expect_false(inherits(r$value, "error"))
  expect_identical(sf::st_crs(r$value$a), sf::st_crs(4326))
  lines <- capture_spatialkit_log(r <- .ls_catch(ensure_projected(inf, target_crs = 4326)))
  expect_identical(sf::st_crs(r$value), sf::st_crs(4326))
})

# The replay of a fit's lon/lat assumption on CRS-less new data
# (.replay_crs_assumption()) and cell_values_at()'s reading of CRS-less
# `newdata` and `seeds` beside CRS-less `points` in degrees stamped
# EPSG:4326 in the same way.  Predictions came back NA with
# ensure_projected()'s warning that the rows have a location in longitude
# and latitude, and cell_values_at() returned NA rows with no word about
# `newdata`.

# CRS-less degrees on a 9 x 9 lattice, x -80 to -78, y 36 to 38, with a
# predictor and a response that do not depend on the random stream.
.ls_train <- function() {
  g <- expand.grid(i = 0:8, j = 0:8)
  d <- data.frame(x = -80 + g$i / 4, y = 36 + g$j / 4)
  d$z <- sin(g$i) + cos(g$j)
  d$v <- 2 * d$z + (g$i %% 3) / 10
  sf::st_as_sf(d, coords = c("x", "y"))
}
metres_z <- .ls_pts(5e5, 4e6, 2000)
metres_z$z <- seq(-1, 1, length.out = nrow(metres_z))

test_that("the replay of a fit's lon/lat assumption refuses CRS-less coordinates outside the envelope", {
  tr <- .ls_train()
  ep <- suppressWarnings(ensure_projected(tr))
  expect_identical(attr(ep, "crs_assumed"), "EPSG:4326")
  for (m in list(metres_z, y_over))
    .ls_refused(.ls_catch(.replay_crs_assumption(m, ep, "predict.rf_fit")),
                "predict.rf_fit", "newdata")
  .ls_refused(.ls_catch(.replay_crs_assumption(metres_z, ep, "predict_surface", "grid")),
              "predict_surface", "grid")
  # Degrees inside the envelope are still given the assumption, with its
  # warning, and a fit that assumed nothing leaves the layer alone.
  r <- .ls_catch(.replay_crs_assumption(tr[1:5, ], ep, "predict.rf_fit"))
  expect_identical(sf::st_crs(r$value), sf::st_crs(4326))
  expect_match(r$warnings, "interpreting it as EPSG:4326", fixed = TRUE)
  inside <- .ls_pts(200, 10, 50)
  r <- .ls_catch(.replay_crs_assumption(inside, ep, "predict.rf_fit"))
  expect_identical(sf::st_crs(r$value), sf::st_crs(4326))
  plain <- ep; attr(plain, "crs_assumed") <- NULL
  r <- .ls_catch(.replay_crs_assumption(metres_z, plain, "predict.rf_fit"))
  expect_true(is.na(sf::st_crs(r$value)))
  expect_length(r$warnings, 0L)
})

test_that("predict() methods and predict_surface() refuse CRS-less new data in metres after a fit on degrees", {
  tr <- .ls_train()
  ep <- suppressWarnings(ensure_projected(tr))
  # predict.gwr_fit() and predict.bayesian_fit() replay the assumption
  # before they read the model, so a fit object with the fit's data is
  # enough to reach it.
  fake <- function(cls) structure(
    list(engine = structure(list(), class = "brmsfit"), data_sf = ep,
         response_var = "v", predictor_vars = "z"),
    class = c(cls, "spatial_fit"))
  if (requireNamespace("sp", quietly = TRUE))
    .ls_refused(.ls_catch(predict(fake("gwr_fit"), newdata = metres_z)),
                "predict.gwr_fit", "newdata")
  if (requireNamespace("brms", quietly = TRUE))
    .ls_refused(.ls_catch(predict(fake("bayesian_fit"), newdata = metres_z)),
                "predict.bayesian_fit", "newdata")
  skip_if_not_installed("ranger")
  fit <- suppressWarnings(fit_rf_model(tr, "v", "z", num_trees = 20, seed = 1))
  expect_identical(attr(fit$data_sf, "crs_assumed"), "EPSG:4326")
  .ls_refused(.ls_catch(predict(fit, newdata = metres_z)), "predict.rf_fit", "newdata")
  .ls_refused(.ls_catch(predict_surface(fit, grid = metres_z)), "predict_surface", "grid")
  .ls_refused(.ls_catch(predict_surface(fit, n_cells = 25, covariates = metres_z)),
              "predict_surface", "covariates")
  .ls_refused(.ls_catch(predict_surface(fit, grid = fit$data_sf[, "z"],
                                        boundary = .ls_sq(4.9e5, 3.9e6, 2e5))),
              "predict_surface", "boundary")
  # The training rows without a CRS are still placed where the fit placed
  # them.
  r <- .ls_catch(predict(fit, newdata = tr))
  expect_match(r$warnings, "interpreting it as EPSG:4326", fixed = TRUE)
  expect_equal(r$value, predict(fit, newdata = fit$data_sf), tolerance = 1e-12)
})

test_that("cell_values_at() refuses CRS-less newdata and seeds in metres beside CRS-less points in degrees", {
  tr <- .ls_train()
  tr$cell_id <- 1L + (sf::st_coordinates(tr)[, 1] > -79) +
    2L * (sf::st_coordinates(tr)[, 2] > 37)
  # Four cells of one degree, without a CRS, numbered as cell_id is.
  cells <- do.call(rbind, lapply(1:4, function(k) {
    s <- .ls_sq(-80 + ((k - 1L) %% 2L), 36 + ((k - 1L) %/% 2L), 1)
    s$cell_id <- k; s[, "cell_id"]
  }))
  vals <- suppressWarnings(summarize_by_cell(tr, "v", id_col = "cell_id",
                                             estimate = "mean", cells_sf = cells))
  seeds <- .ls_pts(5e5, 4e6, 2000, n = 4L)
  seeds$cell_id <- 1:4
  for (m in list(metres, y_over)) {
    r <- .ls_catch(cell_values_at(m, vals, tr, "v", id_col = "cell_id"))
    # The one warning is ensure_projected()'s, that `points` is read as
    # degrees; none is about `newdata`.
    expect_length(r$warnings, 1L)
    expect_match(r$warnings, "^ensure_projected\\(\\): CRS is missing; assuming EPSG:4326")
    r$warnings <- character()
    .ls_refused(r, "cell_values_at", "newdata")
  }
  # The refusal of `seeds` comes before anything is read off them; with
  # `newdata` in degrees it is the call's only error.
  r <- .ls_catch(cell_values_at(tr[1:3, ], vals, tr, "v", id_col = "cell_id",
                                seeds = seeds))
  expect_s3_class(r$value, "error")
  expect_match(conditionMessage(r$value), paste0(
    "^cell_values_at\\(\\): `seeds` has no CRS and its coordinates lie outside ",
    "longitude -180 to 360 and latitude -90 to 90 \\(xmin=500000, xmax=501500, "))
  # CRS-less degrees beside them are read as `points` is, as before.
  # (-79.75, 36.25) and (-78.25, 37.75), inside cells 1 and 4.
  r <- .ls_catch(cell_values_at(tr[c(11, 71), ], vals, tr, "v", id_col = "cell_id"))
  expect_false(inherits(r$value, "error"))
  expect_identical(r$value$cell_id, c(1L, 4L))
})
