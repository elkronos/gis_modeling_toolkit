# A refusal or warning that sat on one path holds on its sibling (BUG-148,
# BUG-149, BUG-150, BUG-175, BUG-213, BUG-214, BUG-215): a layer with no
# complete predictor row warns as a split half without one does; a `values`
# of empty polygons is refused as one without polygons is; points none of
# which has a location are refused with or without a boundary, and zero
# cells still give the count of points in none; every argument a seeding
# method ignores is reported; the refusal of zero cells holds after the clip;
# the snap of cells in degrees takes its tolerance and distances on the
# sphere, as its margin does.

# A call with its warnings kept and muffled, and the log lines it wrote.
.gs_catch <- function(expr) {
  got <- list()
  value <- NULL
  logged <- capture_spatialkit_log(
    value <- withCallingHandlers(expr, warning = function(w) {
      got[[length(got) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }))
  list(value = value, warnings = got, logged = logged)
}
.gs_of_class <- function(r, cls) Filter(function(w) inherits(w, cls), r$warnings)
.gs_texts <- function(r) vapply(r$warnings, conditionMessage, character(1))
.gs_error <- function(expr) tryCatch({ expr; NULL }, error = function(e) e)
.gs_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

.gs_pts <- function(n = 200, seed = 1, crs = 32632) {
  withr::local_seed(seed)
  d <- sf::st_as_sf(data.frame(x = stats::runif(n, 0, 1e4), y = stats::runif(n, 0, 1e4)),
                    coords = c("x", "y"), crs = crs)
  d$z <- sin(sf::st_coordinates(d)[, 1] / 2000) + stats::rnorm(n, sd = 0.3)
  d$w <- stats::rnorm(n)
  d$wNA <- NA_real_
  d
}
.gs_square <- function(x0, y0, w, crs = 32632)
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w), c(x0, y0 + w), c(x0, y0)))), crs = crs))

# ---- BUG-148: resolution_profile() with no complete predictor row ----------

.gs_profile <- function(d, pred, ...)
  .gs_catch(resolution_profile(d, "z", predictor_vars = pred, n_levels = 6,
                               estimate = NULL, partition = "kmeans", ...))

test_that("a layer with a response and no complete predictor row warns spatialkit_no_scored_rows, with or without the split", {
  d <- .gs_pts()
  for (split in c(FALSE, TRUE)) {
    r <- if (split) .gs_profile(d, "wNA", select_on = "split") else .gs_profile(d, "wNA")
    info <- paste("split", split)
    nsr <- .gs_of_class(r, "spatialkit_no_scored_rows")
    expect_length(nsr, 1L)
    msg <- conditionMessage(nsr[[1L]])
    expect_match(msg, "no row of the layer has a finite value of 'z' and of every predictor",
                 fixed = TRUE, info = info)
    expect_match(msg, "`rss`, `cp`, `moran_i` and `moran_z` are NA at every level",
                 fixed = TRUE, info = info)
    expect_match(msg, "The 200 rows with a finite value of 'z' each have a predictor",
                 fixed = TRUE, info = info)
    expect_false(grepl("selection half", msg), info = info)
    expect_identical(nsr[[1L]]$reason, "no row with a response and every predictor")
    expect_identical(nsr[[1L]]$n_response, 200L)
    expect_identical(nsr[[1L]]$n, 0L)
    expect_identical(nsr[[1L]]$selection,
                     if (split) length(attr(r$value, "split")$selection) else NA_integer_)
    # No fit on zero rows: neither base R's bare warning nor the log line
    # that said the raw response was scored.
    expect_false(any(grepl("not a multiple of vector length", .gs_texts(r))), info = info)
    expect_false(log_has(r$logged, "scoring the raw response"), info = info)
    p <- r$value
    expect_true(all(is.na(p$rss)) && all(is.na(p$cp)) && all(is.na(p$moran_z)))
    b <- attr(p, "bounds")
    expect_identical(b$n_response, 0L)
    expect_identical(b$n_response_only, 200L)
    # summary() and select_resolution() name the predictors as the cause.
    e <- .gs_error(suppressWarnings(summary(p)))
    expect_match(conditionMessage(e), paste0(
      "No row of the layer has a response and every predictor (the 200 rows with a ",
      "response each have a predictor that is missing or not finite)"), fixed = TRUE)
    expect_false(grepl("cp needs a response", conditionMessage(e), fixed = TRUE))
    for (crit in c("cp", "moran_z")) {
      e <- .gs_error(suppressWarnings(select_resolution(p, criterion = crit, diagnostic = TRUE)))
      expect_match(conditionMessage(e), sprintf(paste0(
        "`%s` is NA at every level because no row of the layer has a response and every ",
        "predictor: the 200 rows with a response each have a predictor"), crit),
        fixed = TRUE, info = paste(info, crit))
    }
  }
})

test_that("a split half with no complete row keeps its warning and text, without base R's warning", {
  d <- .gs_pts()
  p0 <- .gs_quiet(resolution_profile(d, "z", predictor_vars = "w", n_levels = 6,
                                     estimate = NULL, partition = "kmeans",
                                     select_on = "split"))
  d$wS <- d$w
  d$wS[attr(p0, "split")$selection] <- NA
  r <- .gs_profile(d, "wS", select_on = "split")
  nsr <- .gs_of_class(r, "spatialkit_no_scored_rows")
  expect_length(nsr, 1L)
  expect_match(conditionMessage(nsr[[1L]]), paste0(
    "no row of the selection half has a finite value of 'z' and of every predictor, so ",
    "`rss`, `cp`, `moran_i` and `moran_z`, which read the selection half, are NA at every ",
    "level."), fixed = TRUE)
  expect_identical(nsr[[1L]]$reason, "no response in the selection half")
  expect_false(any(grepl("not a multiple of vector length", .gs_texts(r))))
  expect_false(log_has(r$logged, "scoring the raw response"))
  # The layer has complete rows: the predictors are not the cause.
  b <- attr(r$value, "bounds")
  expect_null(b$n_response_only)
  expect_gt(b$n_response, 0L)
  expect_false(.no_complete_row(r$value))
  e <- .gs_error(suppressWarnings(select_resolution(r$value, criterion = "cp")))
  expect_match(conditionMessage(e), "no row of the selection half has a response and every predictor",
               fixed = TRUE)
})

test_that("a layer with complete rows, or without a response, is not told of the predictors", {
  d <- .gs_pts()
  r <- .gs_profile(d, "w")
  expect_length(.gs_of_class(r, "spatialkit_no_scored_rows"), 0L)
  expect_identical(attr(r$value, "variable"), "residuals")
  expect_null(attr(r$value, "bounds")$n_response_only)
  expect_true(all(is.finite(r$value$cp)))
  # No row with a response: the predictors are not the cause, and nothing
  # is fitted on zero rows.
  d$z <- NA_real_
  r <- .gs_profile(d, "wNA")
  expect_length(.gs_of_class(r, "spatialkit_no_scored_rows"), 0L)
  expect_false(any(grepl("not a multiple of vector length", .gs_texts(r))))
  expect_null(attr(r$value, "bounds")$n_response_only)
  e <- .gs_error(suppressWarnings(select_resolution(r$value, criterion = "cp")))
  expect_match(conditionMessage(e), "(it needs a response,", fixed = TRUE)
})

# ---- BUG-149: cell_values_at() on a `values` of empty polygons -------------

.gs_values <- function() {
  d <- .gs_pts(240, seed = 7)
  rows <- d[1:200, ]; held <- d[201:240, ]
  tess <- .gs_quiet(build_tessellation(rows, boundary = .gs_square(0, 0, 1e4),
                                       method = "hex", approx_n_cells = 20, quiet = TRUE))
  asg <- .gs_quiet(assign_features_to_polygons(rows, tess$cells, polygon_id_col = "cell_id",
                                               quiet = TRUE))
  values <- .gs_quiet(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                        estimate = "pooled"))
  list(values = values, asg = asg, held = held)
}

test_that("cell_values_at() without seeds refuses a `values` none of whose polygons has a location", {
  V <- .gs_values()
  refusal <- paste0("cell_values_at(): `values` carries no cell polygons, so a location ",
                    "cannot be placed in a cell.")
  for (empty in list(sf::st_polygon(), sf::st_multipolygon())) {
    ve <- V$values
    sf::st_geometry(ve) <- sf::st_sfc(lapply(seq_len(nrow(ve)), function(i) empty),
                                      crs = sf::st_crs(ve))
    e <- .gs_error(cell_values_at(V$held, ve, V$asg, "z", id_col = "cell_id"))
    expect_s3_class(e, "simpleError")
    expect_match(conditionMessage(e), refusal, fixed = TRUE)
  }
  # The sibling refusal keeps its text.
  vp <- V$values
  sf::st_geometry(vp) <- suppressWarnings(sf::st_centroid(sf::st_geometry(vp)))
  e <- .gs_error(cell_values_at(V$held, vp, V$asg, "z", id_col = "cell_id"))
  expect_match(conditionMessage(e), refusal, fixed = TRUE)
})

test_that("cell_values_at() runs on a `values` with some empty polygons, as before", {
  V <- .gs_values()
  whole <- .gs_quiet(cell_values_at(V$held, V$values, V$asg, "z", id_col = "cell_id"))
  v1 <- V$values
  gone <- unique(whole$cell_id[!is.na(whole$cell_id)])[1:2]
  at <- which(v1$cell_id %in% gone)
  sf::st_geometry(v1)[at] <- sf::st_sfc(lapply(at, function(i) sf::st_polygon()),
                                        crs = sf::st_crs(v1))
  some <- .gs_quiet(cell_values_at(V$held, v1, V$asg, "z", id_col = "cell_id"))
  expect_identical(nrow(some), nrow(V$held))
  moved <- whole$cell_id %in% gone
  expect_true(any(moved))
  expect_true(all(some$filled[moved]))
  expect_identical(some$cell_id[!moved], whole$cell_id[!moved])
  expect_identical(some$value[!moved], whole$value[!moved])
})

test_that("cell_values_at() runs on a `values` of multipolygons whose first part is empty", {
  # sf calls such a multipolygon empty, though its second part is the cell.
  V <- .gs_values()
  whole <- .gs_quiet(cell_values_at(V$held, V$values, V$asg, "z", id_col = "cell_id"))
  vm <- V$values
  g <- sf::st_geometry(vm)
  sf::st_geometry(vm) <- sf::st_sfc(lapply(seq_along(g), function(i)
    sf::st_multipolygon(list(list(), unclass(g[[i]])))), crs = sf::st_crs(vm))
  expect_true(all(sf::st_is_empty(vm)))
  got <- .gs_quiet(cell_values_at(V$held, vm, V$asg, "z", id_col = "cell_id"))
  expect_identical(got$cell_id, whole$cell_id)
  expect_identical(got$value, whole$value)
  expect_identical(got$filled, whole$filled)
})

# ---- BUG-150: Voronoi cells of points none of which has a location ---------

test_that("create_voronoi_polygons() and build_tessellation(voronoi) refuse points none of which has a location, with or without a boundary", {
  d <- .gs_pts(20)
  b <- .gs_square(0, 0, 1e4)
  de <- d
  sf::st_geometry(de) <- sf::st_sfc(lapply(1:20, function(i) sf::st_point()), crs = 32632)
  far <- d
  sf::st_geometry(far) <- sf::st_sfc(lapply(1:20, function(i) sf::st_point(c(1e200, i))),
                                     crs = 32632)
  refusal <- paste0("create_voronoi_polygons(): no point of `points_sf` has a location ",
                    "(each is an empty point, or has a coordinate that is not finite or ",
                    "lies beyond 1e150), so no cell can be built.")
  for (pts in list(de, far)) {
    calls <- list(
      function() build_tessellation(pts, boundary = b, method = "voronoi", quiet = TRUE),
      function() build_tessellation(pts, method = "voronoi", quiet = TRUE),
      function() create_voronoi_polygons(pts, boundary = b, quiet = TRUE),
      function() create_voronoi_polygons(pts, quiet = TRUE))
    for (k in seq_along(calls)) {
      r <- tryCatch(.gs_catch(calls[[k]]()), error = function(e) e)
      expect_s3_class(r, "error")
      expect_identical(conditionMessage(r), refusal, info = paste("call", k))
    }
  }
  # Seeds passed on by get_voronoi_seeds(method = "provided").
  se <- suppressWarnings(get_voronoi_seeds(method = "provided", seeds = de[1:5, ]))
  e <- .gs_error(.gs_quiet(build_tessellation(se, method = "voronoi", quiet = TRUE)))
  expect_identical(conditionMessage(e), refusal)
})

test_that("points the move to the working CRS leaves without a location are refused, with or without a boundary", {
  # Latitudes 91 to 95 in EPSG:4326 have a location as given; the move to
  # EPSG:3857 places none of them.
  dl <- sf::st_as_sf(data.frame(x = seq(-90, -80, length.out = 20),
                                y = seq(91, 95, length.out = 20)),
                     coords = c("x", "y"), crs = 4326)
  bl <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(xmin = -91, ymin = 80, xmax = -79,
                                                         ymax = 89), crs = sf::st_crs(4326))))
  refusal <- function(f) paste0(f, ": the move of `points_sf` to the CRS the cells ",
                                "are built in placed none of its points, so no cell can be built.")
  calls <- list(
    function() create_voronoi_polygons(dl, quiet = TRUE),
    function() create_voronoi_polygons(dl, boundary = bl, quiet = TRUE),
    function() build_tessellation(dl, method = "voronoi", quiet = TRUE),
    function() build_tessellation(dl, boundary = bl, method = "voronoi", quiet = TRUE))
  callers <- rep(c("create_voronoi_polygons()", "build_tessellation(voronoi)"), each = 2L)
  for (k in seq_along(calls)) {
    texts <- character(0)
    e <- tryCatch(withCallingHandlers(calls[[k]](), warning = function(w) {
      texts <<- c(texts, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    expect_s3_class(e, "simpleError")
    expect_identical(conditionMessage(e), refusal(callers[k]), info = callers[k])
    expect_false(any(grepl("contains no geometries", texts, fixed = TRUE)), info = callers[k])
    expect_false(any(grepl("BufferOp", c(texts, conditionMessage(e)), fixed = TRUE)),
                 info = callers[k])
  }
})

test_that("one point without a location among others still gives its cell count and warning", {
  d <- .gs_pts(20)
  b <- .gs_square(0, 0, 1e4)
  d1 <- d
  sf::st_geometry(d1)[1] <- sf::st_sfc(sf::st_point(), crs = 32632)
  r <- .gs_catch(build_tessellation(d1, boundary = b, method = "voronoi", quiet = TRUE))
  un <- .gs_of_class(r, "spatialkit_unassigned")
  expect_length(un, 1L)
  expect_match(conditionMessage(un[[1L]]), "1 of 20 point(s) (5.0%) fall in no cell", fixed = TRUE)
  expect_identical(nrow(r$value$cells), 19L)
  # The hexagons of a layer with no location keep their warning and text.
  de <- d
  sf::st_geometry(de) <- sf::st_sfc(lapply(1:20, function(i) sf::st_point()), crs = 32632)
  r <- .gs_catch(build_tessellation(de, boundary = b, method = "hex", approx_n_cells = 10,
                                    quiet = TRUE))
  un <- .gs_of_class(r, "spatialkit_unassigned")
  expect_length(un, 1L)
  expect_match(conditionMessage(un[[1L]]), paste0(
    "20 of 20 point(s) (100.0%) fall in no cell and are NA in `index`; `params$n_outside` ",
    "holds the count. They have no location"), fixed = TRUE)
})

test_that("the index of points in no cell at all gives the spatialkit_unassigned count", {
  d <- .gs_pts(20)
  cells0 <- sf::st_sf(cell_id = integer(0), geometry = sf::st_sfc(crs = 32632))
  r <- .gs_catch(.build_point_cell_index(d, cells0))
  un <- .gs_of_class(r, "spatialkit_unassigned")
  expect_length(un, 1L)
  expect_identical(un[[1L]]$n, 20L)
  expect_identical(un[[1L]]$total, 20L)
  expect_match(conditionMessage(un[[1L]]), "20 of 20 point(s) (100.0%) fall in no cell",
               fixed = TRUE)
  expect_identical(as.integer(r$value), rep(NA_integer_, 20))
  expect_identical(attr(r$value, "snapped"),
                   list(n = 0L, which = integer(0), distance = numeric(0)))
  # Asked not to warn, it does not; and no point gives an empty index.
  r <- .gs_catch(.build_point_cell_index(d, cells0, warn = FALSE))
  expect_length(r$warnings, 0L)
  r <- .gs_catch(.build_point_cell_index(d[0, ], cells0))
  expect_length(r$warnings, 0L)
  expect_length(r$value, 0L)
})

# ---- BUG-175: get_voronoi_seeds() arguments a method does not read ---------

test_that("get_voronoi_seeds() reports `seeds` under kmeans and random, and `sample_points` under random", {
  bnd <- .gs_square(0, 0, 1e5, crs = 32617)
  sp <- withr::with_seed(4, sf::st_as_sf(
    data.frame(x = stats::runif(200, 0, 1e5), y = stats::runif(200, 0, 1e5)),
    coords = c("x", "y"), crs = 32617))
  seeds_msg <- function(m) sprintf(paste0("get_voronoi_seeds(): `seeds` is read by method = ",
                                          "\"provided\" only; method = \"%s\" ignores it."), m)
  sp_msg <- paste0("get_voronoi_seeds(): `sample_points` is not read by method = \"random\", ",
                   "which places its seeds in `boundary`; it is ignored.")
  for (m in c("kmeans", "random")) {
    plain <- suppressWarnings(get_voronoi_seeds(bnd, method = m, n = 6, sample_points = sp,
                                                set_seed = 1))
    r <- .gs_catch(get_voronoi_seeds(bnd, method = m, n = 6, sample_points = sp,
                                     seeds = sp[1:6, ], set_seed = 1))
    expect_identical(sum(.gs_texts(r) == seeds_msg(m)), 1L, info = m)
    expect_identical(r$value, plain, info = m)
    # `seed = 1` is matched to `seeds`, and is reported as such.
    r <- .gs_catch(get_voronoi_seeds(bnd, method = m, n = 6, sample_points = sp, seed = 1))
    expect_identical(sum(.gs_texts(r) == seeds_msg(m)), 1L, info = m)
  }
  r <- .gs_catch(get_voronoi_seeds(bnd, method = "random", n = 6, sample_points = sp,
                                   set_seed = 1))
  expect_identical(.gs_texts(r), sp_msg)
  expect_identical(r$value, get_voronoi_seeds(bnd, method = "random", n = 6, set_seed = 1))
  # Nothing new where the argument is read, or not passed.
  expect_length(.gs_catch(get_voronoi_seeds(bnd, method = "provided", seeds = sp[1:6, ]))$warnings, 0L)
  expect_length(.gs_catch(get_voronoi_seeds(bnd, method = "kmeans", n = 6, sample_points = sp,
                                            set_seed = 1))$warnings, 0L)
  expect_length(.gs_catch(get_voronoi_seeds(bnd, method = "random", n = 6,
                                            set_seed = 1))$warnings, 0L)
  # A call the method refuses is told only that: a random seeding without a
  # boundary, and a k-means seeding with neither a boundary nor a cloud.
  refused <- function(expr) {
    texts <- character(0)
    e <- tryCatch(withCallingHandlers(expr, warning = function(w) {
      texts <<- c(texts, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    list(error = e, warnings = texts)
  }
  r <- refused(get_voronoi_seeds(method = "random", n = 6, sample_points = sp))
  expect_identical(conditionMessage(r$error), "get_voronoi_seeds(): random seeds require 'boundary'.")
  expect_identical(r$warnings, character(0))
  r <- refused(get_voronoi_seeds(method = "random", n = 6, seeds = sp[1:6, ]))
  expect_identical(conditionMessage(r$error), "get_voronoi_seeds(): random seeds require 'boundary'.")
  expect_identical(r$warnings, character(0))
  r <- refused(get_voronoi_seeds(method = "kmeans", n = 6, seeds = sp[1:6, ]))
  expect_identical(conditionMessage(r$error),
                   "get_voronoi_seeds(): k-means seeds require 'boundary' (or provide 'sample_points').")
  expect_identical(r$warnings, character(0))
  # A k-means seeding from a cloud alone runs, and is told of `seeds`.
  r <- .gs_catch(get_voronoi_seeds(method = "kmeans", n = 6, sample_points = sp,
                                   seeds = sp[1:6, ], set_seed = 1))
  expect_identical(.gs_texts(r), seeds_msg("kmeans"))
  expect_identical(nrow(r$value), 6L)
})

# ---- BUG-213: the clip of a grid that leaves no cell -----------------------

.gs_wide <- function() sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(-179, -80), c(179, -80), c(179, 80), c(-179, 80), c(-179, -80)))), crs = 4326))

test_that("create_grid_polygons() refuses a clip that leaves zero cells", {
  suppressMessages(sf::sf_use_s2(TRUE))
  bw <- .gs_wide()
  refusal <- "create_grid_polygons(): the clip to `boundary` left zero cells."
  cause <- "In degrees s2 reads each edge of `boundary` as the shorter great-circle arc"
  e <- .gs_error(.gs_quiet(create_grid_polygons(bw, cellsize = c(358 / 3, 10), crs = 4326,
                                                quiet = TRUE)))
  expect_s3_class(e, "simpleError")
  expect_match(conditionMessage(e), refusal, fixed = TRUE)
  expect_match(conditionMessage(e), cause, fixed = TRUE)
  pw <- withr::with_seed(3, sf::st_as_sf(data.frame(x = stats::runif(50, -179, 179),
                                                     y = stats::runif(50, -80, 80)),
                                          coords = c("x", "y"), crs = 4326))
  e <- .gs_error(.gs_quiet(build_tessellation(pw, boundary = bw, method = "square",
                                              cellsize = c(358 / 3, 10), crs = 4326,
                                              quiet = TRUE)))
  expect_match(conditionMessage(e), refusal, fixed = TRUE)
  # Without the clip, and with s2 off, the grid is as it was.
  expect_identical(nrow(.gs_quiet(create_grid_polygons(bw, cellsize = c(358 / 3, 10),
                                                       crs = 4326, clip = FALSE,
                                                       quiet = TRUE))), 48L)
  local({
    suppressMessages(sf::sf_use_s2(FALSE))
    on.exit(suppressMessages(sf::sf_use_s2(TRUE)), add = TRUE)
    expect_identical(nrow(.gs_quiet(create_grid_polygons(bw, cellsize = c(358 / 3, 10),
                                                         crs = 4326, quiet = TRUE))), 48L)
  })
  # A one degree box at 0.1 degrees keeps its 100 cells.
  one <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(xmin = -93, ymin = 41, xmax = -92,
                                                          ymax = 42), crs = sf::st_crs(4326))))
  expect_identical(nrow(.gs_quiet(create_grid_polygons(one, cellsize = 0.1, crs = 4326,
                                                       quiet = TRUE))), 100L)
})

test_that("a projected grid whose clip leaves zero cells is refused without the s2 clause", {
  thin <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1e-9), c(0, 1e-9), c(0, 0)))), crs = 32632))
  e <- .gs_error(.gs_quiet(create_grid_polygons(thin, cellsize = 100, quiet = TRUE)))
  expect_s3_class(e, "simpleError")
  expect_identical(conditionMessage(e),
                   "create_grid_polygons(): the clip to `boundary` left zero cells.")
})

# ---- BUG-214 / BUG-215: the snap of cells in degrees -----------------------

# A 20 x 10 lattice at 0.05 degrees from (-93, 41), two points just south of
# its box, the box as the boundary.
.gs_degrees <- function() {
  g <- expand.grid(i = 0:19, j = 0:9)
  p <- sf::st_as_sf(data.frame(lon = -93 + g$i * 0.05, lat = 41 + g$j * 0.05, z = 1),
                    coords = c("lon", "lat"), crs = 4326)
  out <- sf::st_as_sf(data.frame(lon = c(-93 + 0.2, -93 + 0.4), lat = c(41 - 1e-6, 41 - 2e-6),
                                 z = 1), coords = c("lon", "lat"), crs = 4326)
  list(pts = rbind(p, out), box = sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(p))))
}
# The tolerance of the snap for 0.1 degree cells, as an arc in metres.
.gs_tol_m <- 0.1 / 1000 * pi / 180 * 6371010

test_that("with s2 on, cells in degrees snap a point within a thousandth of the cell width as an arc", {
  suppressMessages(sf::sf_use_s2(TRUE))
  L <- .gs_degrees()
  r <- .gs_catch(build_tessellation(L$pts, boundary = L$box, method = "square",
                                    cellsize = 0.1, crs = 4326, quiet = TRUE))
  tz <- r$value
  expect_identical(nrow(tz$cells), 50L)
  # The lattice point (-92.05, 41), 1.2 m from its cell on the sphere.
  k <- 20L
  expect_equal(unname(sf::st_coordinates(L$pts)[k, ]), c(-92.05, 41))
  d_all <- as.numeric(sf::st_distance(L$pts[k, ], tz$cells))
  expect_identical(tz$index[k], tz$cells$cell_id[which.min(d_all)])
  expect_identical(tz$params$snapped$n, 1L)
  expect_identical(tz$params$snapped$which, k)
  expect_equal(tz$params$snapped$distance, min(d_all))
  expect_gt(tz$params$snapped$distance, 1); expect_lt(tz$params$snapped$distance, 1.5)
  # Every other point in no cell lies 21 m or more outside, beyond the
  # tolerance of 11.1 m, and stays NA.
  na <- which(is.na(tz$index))
  expect_identical(na, c(2:19, 201L, 202L))
  d_na <- apply(as.matrix(sf::st_distance(L$pts[na, ], tz$cells)), 1L, min)
  expect_true(all(as.numeric(d_na) > 21))
  expect_true(all(as.numeric(d_na) > .gs_tol_m))
  un <- .gs_of_class(r, "spatialkit_unassigned")
  expect_identical(un[[1L]]$n, 20L)
})

test_that("with s2 off, cells in degrees are snapped by distances on the sphere, without lwgeom", {
  L <- .gs_degrees()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(TRUE)), add = TRUE)
  tz <- .gs_quiet(build_tessellation(L$pts, boundary = L$box, method = "square",
                                     cellsize = 0.1, crs = 4326, quiet = TRUE))
  expect_false(isTRUE(sf::sf_use_s2()))
  expect_identical(nrow(tz$cells), 50L)
  # The same cells indexed again with s2 off: the snap's distances are those
  # s2 measures.
  ix <- .gs_quiet(.build_point_cell_index(L$pts, tz$cells))
  expect_identical(as.integer(ix), tz$index)
  s <- attr(ix, "snapped")
  expect_identical(s, tz$params$snapped)
  suppressMessages(sf::sf_use_s2(TRUE))
  for (i in seq_along(s$which)) {
    d <- as.numeric(sf::st_distance(L$pts[s$which[i], ], tz$cells))
    expect_equal(s$distance[i], min(d))
    expect_lte(s$distance[i], .gs_tol_m)
    expect_identical(tz$index[s$which[i]], tz$cells$cell_id[which.min(d)])
  }
  suppressMessages(sf::sf_use_s2(FALSE))
  # The two points south of the box, 0.11 and 0.22 m below the planar edge.
  expect_identical(s$which, c(201L, 202L))
  expect_false(anyNA(tz$index))
})

test_that("projected cells keep their snap rule in the CRS's units", {
  L <- .gs_degrees()
  pts <- sf::st_transform(L$pts, 26915)
  box <- sf::st_transform(L$box, 26915)
  for (s2 in c(TRUE, FALSE)) {
    suppressMessages(sf::sf_use_s2(s2))
    tz <- .gs_quiet(build_tessellation(pts, boundary = box, method = "square",
                                       cellsize = 8000, quiet = TRUE))
    cw <- stats::median(vapply(sf::st_geometry(tz$cells), function(g) {
      bb <- sf::st_bbox(g); as.numeric(bb[["xmax"]] - bb[["xmin"]]) }, numeric(1)))
    d <- apply(as.matrix(sf::st_distance(pts, tz$cells)), 1L, min)
    s <- tz$params$snapped
    expect_true(all(d[s$which] > 0 & d[s$which] <= cw / 1000), info = paste("s2", s2))
    expect_equal(s$distance, unname(d[s$which]), info = paste("s2", s2))
    expect_true(all(d[is.na(tz$index)] > cw / 1000), info = paste("s2", s2))
  }
  suppressMessages(sf::sf_use_s2(TRUE))
})

test_that("with s2 off, the nearest cell of the snap in degrees is the nearest on the sphere", {
  # A 5 x 5 lattice of 0.05 degree cells at 60N without its middle cell, and
  # a point in that hole 8e-5 degrees east of the west cell and 6e-5 north
  # of the south cell: 4.4 m from the west cell on the sphere, 6.7 m from the
  # south cell, within the tolerance of 5.6 m of the one only.  In degrees
  # taken as planar the south cell is the nearer.
  sq <- function(x0, y0, w) sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w), c(x0, y0 + w), c(x0, y0))))
  g <- expand.grid(i = 0:4, j = 0:4)
  g <- g[!(g$i == 2 & g$j == 2), ]
  cells <- sf::st_sf(cell_id = seq_len(nrow(g)), geometry = sf::st_sfc(
    lapply(seq_len(nrow(g)), function(k) sq(10 + g$i[k] * 0.05, 60 + g$j[k] * 0.05, 0.05)),
    crs = 4326))
  west <- cells$cell_id[g$i == 1 & g$j == 2]
  south <- cells$cell_id[g$i == 2 & g$j == 1]
  p <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(10.1 + 8e-5, 60.1 + 6e-5)), crs = 4326))
  on.exit(suppressMessages(sf::sf_use_s2(TRUE)), add = TRUE)
  for (s2 in c(TRUE, FALSE)) {
    suppressMessages(sf::sf_use_s2(s2))
    ix <- .gs_quiet(.build_point_cell_index(p, cells))
    expect_identical(as.integer(ix), west, info = paste("s2", s2))
    s <- attr(ix, "snapped")
    expect_identical(s$n, 1L, info = paste("s2", s2))
    expect_gt(s$distance, 4); expect_lt(s$distance, 5)
  }
  expect_identical(.gs_quiet(sf::st_nearest_feature(p, cells)), south)
})
