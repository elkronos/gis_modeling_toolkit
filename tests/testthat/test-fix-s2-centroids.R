# tests/testthat/test-fix-s2-centroids.R
# ---------------------------------------------------------------------------
# Centroids in longitude and latitude where s2 refuses a polygon.
#
# s2 refuses a ring with a vertex that repeats the one before it (which GEOS
# accepts), and rings whose edges cross (a bowtie, a figure eight).  Taken
# in one call, one refused feature took the centroids of the others with it:
#   * .crs_distance_error() kept GEOS's interior point, on the far side of
#     the globe, for every feature across 180 degrees (BUG-133);
#   * .lonlat_centre() fell back to the mean of every vertex, in silence,
#     counting each ring's closing vertex twice (BUG-134);
#   * a fold probe saved without `kind` could not be retaken, and the
#     provenance check passed in silence (BUG-135);
#   * coerce_to_points(mode = "centroid") and prep_model_data(pointize =
#     "centroid") stopped with s2's bare message (BUG-130).
# The features s2 refuses are now repaired one by one (.s2_mend()): without
# their repeated vertices, which leaves the ring's points as they were, and
# for the CRS chooser also with crossing edges split.
# ---------------------------------------------------------------------------

# s2 on (or off) for the length of the calling test.
.sc_s2 <- function(on = TRUE, env = parent.frame()) {
  old <- suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}
# Warnings, messages and the log lines of the CRS chooser kept off the
# console.
.sc_quiet <- function(expr) {
  r <- NULL
  capture_spatialkit_log(r <- suppressWarnings(suppressMessages(expr)))
  r
}
.sc_sq <- function(x0, y0, d = 0.1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d),
                            c(x0, y0))))
.sc_rep <- function(x0, y0, d = 0.1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0), c(x0 + d, y0 + d),
                            c(x0, y0 + d), c(x0, y0))))
.sc_bow <- function(x0, y0, d = 0.1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0 + d), c(x0 + d, y0), c(x0, y0 + d),
                            c(x0, y0))))
.sc_fig8 <- function(x0, y0, d = 0.1) {
  m <- d / 2
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + m, y0 + m), c(x0 + d, y0), c(x0 + d, y0 + d),
                            c(x0 + m, y0 + m), c(x0, y0 + d), c(x0, y0))))
}
# nx by ny squares of side d from (x0, y0), row by row; rows `bad` drawn by
# `maker`.
.sc_grid <- function(nx, ny, bad = integer(0), maker = .sc_rep, x0 = -94, y0 = 39,
                     d = 0.1, crs = 4326) {
  g <- list(); k <- 0L
  for (j in seq_len(ny)) for (i in seq_len(nx)) {
    k <- k + 1L
    g[[k]] <- if (k %in% bad) maker(x0 + (i - 1) * d, y0 + (j - 1) * d, d)
              else .sc_sq(x0 + (i - 1) * d, y0 + (j - 1) * d, d)
  }
  sf::st_sf(poly_id = seq_along(g), geometry = sf::st_sfc(g, crs = crs))
}
# The log lines of `expr`, without their time stamps.
.sc_log <- function(expr) {
  l <- capture_spatialkit_log(suppressWarnings(suppressMessages(expr)))
  sub("^[A-Z]+ \\[[^]]*\\] ", "", l)
}

# ---- BUG-133: the wide features of the distance score -----------------------

.sc_fiji <- function(rep = FALSE, bow = FALSE) {
  m <- if (bow) rbind(c(177, -19), c(-178, -15), c(-178, -19), c(177, -15), c(177, -19))
       else rbind(c(177, -19), c(-178, -19), c(-178, -15), c(177, -15), c(177, -19))
  if (rep) m <- m[c(1, 2, 2, 3, 4, 5), ]
  sf::st_sf(id = 1:3, geometry = sf::st_sfc(sf::st_polygon(list(m)), .sc_sq(178, -18, 0.5),
                                            .sc_sq(-179, -17, 0.5), crs = 4326))
}

test_that("a feature across 180 with a repeated vertex scores as the same feature without it", {
  .sc_s2(TRUE)
  plain <- .sc_fiji(); rep <- .sc_fiji(rep = TRUE)
  expect_false(all(sf::st_is_valid(rep)))
  a <- .sc_quiet(ensure_projected(plain)); b <- .sc_quiet(ensure_projected(rep))
  expect_true(sf::st_crs(a) == sf::st_crs(b))
  ca <- attr(a, "crs_choice"); cb <- attr(b, "crs_choice")
  expect_identical(ca$name, cb$name)
  expect_identical(ca$crs, cb$crs)
  expect_equal(cb$distance_error, ca$distance_error, tolerance = 1e-9)
  expect_lt(cb$distance_error, 0.001)          # 0.2537 before
  ea <- .crs_distance_error(plain, sf::st_crs(a)); eb <- .crs_distance_error(rep, sf::st_crs(a))
  expect_equal(eb, ea, tolerance = 1e-9)
  # The same with s2 off in the session: the centroid is s2's whatever it says.
  .sc_s2(FALSE)
  expect_equal(.crs_distance_error(rep, sf::st_crs(a)), ea, tolerance = 1e-9)
})

test_that("a feature across 180 whose edges cross is scored on its rebuild", {
  .sc_s2(TRUE)
  bow <- .sc_fiji(bow = TRUE)
  g   <- sf::st_geometry(bow)[1]
  expect_false(sf::st_is_valid(g))
  # Its representative point is the centroid of its rebuild, which s2
  # accepts; GEOS's interior point lies at (-89.25, -17), on the far side.
  r <- .s2_centroids(g, split = TRUE)
  expect_identical(r$refused, integer(0))
  reb <- sf::st_make_valid(g, split_crossing_edges = TRUE)
  expect_equal(sf::st_coordinates(r$points), sf::st_coordinates(sf::st_centroid(reb)),
               tolerance = 1e-12)
  # Measured on every point (a sample large enough to keep them all), the
  # far-side point scored 9.5%.
  crs <- sf::st_crs(.sc_quiet(ensure_projected(.sc_fiji())))
  expect_lt(.crs_distance_error(bow, crs, max_n = 400L), 0.001)
})

test_that("only a wide feature s2 refuses after both retries keeps GEOS's point", {
  .sc_s2(TRUE)
  local_mocked_bindings(.split_crossing_edges = function(g1) NULL)
  g <- sf::st_geometry(rbind(.sc_fiji(rep = TRUE)[1, ], .sc_fiji(bow = TRUE)[1, ],
                             .sc_fiji()[1, ]))
  r <- .s2_centroids(g, split = TRUE)
  expect_identical(r$refused, 2L)
  expect_true(sf::st_is_empty(r$points[2]))
  ref <- sf::st_centroid(g[3])
  expect_equal(sf::st_coordinates(r$points[1]), sf::st_coordinates(ref), tolerance = 1e-12)
  expect_equal(sf::st_coordinates(r$points[3]), sf::st_coordinates(ref), tolerance = 1e-12)
  # In the score the box with the repeat gets its centroid; with every wide
  # feature sent to GEOS's point it scored 25%.
  two <- sf::st_sf(id = 1:4, geometry = c(sf::st_geometry(.sc_fiji(rep = TRUE)), g[2]))
  one <- sf::st_sf(id = 1:4, geometry = c(sf::st_geometry(.sc_fiji()), g[2]))
  crs <- sf::st_crs(.sc_quiet(ensure_projected(.sc_fiji())))
  expect_equal(.crs_distance_error(two, crs), .crs_distance_error(one, crs), tolerance = 1e-9)
  # The refused bowtie keeps GEOS's point in the score, not an empty point.
  # Measured on every point, that far-side point scores 9.5%; left out, the
  # score read 0.04%.  The score built by hand: the box's centroid and the
  # bowtie's GEOS point handed in as the wide features' points.
  got <- .crs_distance_error(one, crs, max_n = 400L)
  expect_gt(got, 0.05)
  geos <- suppressWarnings(sf::st_point_on_surface(g[2]))
  box  <- sf::st_centroid(sf::st_geometry(.sc_fiji())[1])
  local_mocked_bindings(.s2_centroids = function(g, split = TRUE)
    list(points = c(box, geos), refused = integer(0), error = NULL))
  expect_equal(got, .crs_distance_error(one, crs, max_n = 400L), tolerance = 1e-12)
})

test_that("layers whose wide features s2 accepts score as before", {
  .sc_s2(TRUE)
  x <- .sc_fiji()
  crs <- sf::st_crs(.sc_quiet(ensure_projected(x)))
  # The score built by hand as the old code took it: the spherical centroid
  # of the wide feature in one call.
  g <- sf::st_geometry(x)
  r <- .s2_centroids(g[1], split = TRUE)
  expect_identical(r$points, suppressWarnings(sf::st_centroid(g[1])))
  expect_identical(r$refused, integer(0))
  expect_null(r$error)
})

# ---- BUG-134: the centre of the layer --------------------------------------

# 100 one-degree squares over -104..-94 by 30..40, the north-east square drawn
# with 4000 vertices; row `bad` drawn by `maker`.
.sc_hundred <- function(bad = integer(0), maker = .sc_bow) {
  x <- .sc_grid(10, 10, bad = bad, maker = maker, x0 = -104, y0 = 30, d = 1)
  t <- seq(0, 1, length.out = 1001)[-1001]
  m <- rbind(cbind(-95 + t, 39), cbind(-94, 39 + t), cbind(-94 - t, 40), cbind(-95, 40 - t))
  sf::st_geometry(x)[100] <- sf::st_polygon(list(rbind(m, m[1, ])))
  x
}
# The layer with row `i` replaced by the rebuild s2 makes of it with crossing
# edges split: what a caller who repaired it would pass.
.sc_rebuilt <- function(x, i) {
  g <- sf::st_geometry(x)
  sf::st_geometry(x)[i] <- sf::st_make_valid(g[i], split_crossing_edges = TRUE)
  x
}

test_that("one bowtie does not move the centre to the mean of the vertices", {
  .sc_s2(TRUE)
  ctrl <- .sc_hundred()
  c0 <- .lonlat_centre(ctrl)
  d0 <- .sc_quiet(ensure_projected(ctrl))
  for (bad in c(1L, 55L)) {
    x  <- .sc_hundred(bad)
    cx <- .lonlat_centre(x)
    # The bowtie is half the area of the square it replaces, so the centre is
    # that of the layer with the bowtie rebuilt, close to the control's.
    expect_equal(cx, .lonlat_centre(.sc_rebuilt(x, bad)), tolerance = 1e-9)
    expect_lt(max(abs(cx - c0)), 0.05)     # 4 degrees before
    # s2 accepts the rebuild, so no feature is left out and nothing is logged.
    expect_length(.sc_log(.lonlat_centre(x)), 0L)
    dx <- .sc_quiet(ensure_projected(x))
    expect_true(sf::st_crs(dx) == sf::st_crs(d0))           # UTM zone 14
    ax <- .sc_quiet(ensure_projected(x, purpose = "area"))
    expect_match(sf::st_crs(ax)$proj4string, "+proj=laea", fixed = TRUE)
  }
  # On 100 plain squares a bowtie moved the centre by 0.1 degree.
  p0 <- .lonlat_centre(.sc_grid(10, 10, x0 = -104, y0 = 30, d = 1))
  p1 <- .lonlat_centre(.sc_grid(10, 10, bad = 55L, maker = .sc_bow, x0 = -104, y0 = 30, d = 1))
  expect_lt(max(abs(p1 - p0)), 0.005)
})

test_that("a bowtie gives one centre from whichever corner its ring starts", {
  .sc_s2(TRUE)
  m0 <- rbind(c(-94.2, 39.5), c(-93.8, 39.9), c(-93.8, 39.5), c(-94.2, 39.9))
  from <- function(k) { m <- m0[c((k:3) + 1, if (k) 1:k), ]; sf::st_polygon(list(rbind(m, m[1, ]))) }
  ctr <- lapply(0:3, function(k)
    .lonlat_centre(sf::st_sfc(from(k), .sc_sq(-93, 39.5, 0.4), crs = 4326)))
  for (k in 2:4) expect_equal(ctr[[k]], ctr[[1]], tolerance = 1e-9)
  crs <- lapply(0:3, function(k) sf::st_crs(.sc_quiet(ensure_projected(
    sf::st_sf(id = 1:2, geometry = sf::st_sfc(from(k), .sc_sq(-93, 39.5, 0.4), crs = 4326)),
    purpose = "area"))))
  for (k in 2:4) expect_true(crs[[k]] == crs[[1]])
  # Where s2 refuses every feature, the vertex mean counts each ring's
  # closing vertex once, so it does not depend on the start either.
  local_mocked_bindings(.split_crossing_edges = function(g1) NULL)
  vm <- lapply(0:3, function(k) .lonlat_centre(sf::st_sfc(from(k), crs = 4326)))
  for (k in 2:4) expect_equal(vm[[k]], vm[[1]], tolerance = 1e-12)
  expect_equal(unname(vm[[1]][1, ]), c(-94, 39.7), tolerance = 1e-3)
})

test_that("features s2 refuses even repaired are left out of the centre, with a log line", {
  .sc_s2(TRUE)
  local_mocked_bindings(.split_crossing_edges = function(g1) NULL)
  x <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(.sc_bow(-94, 39, 0.4), .sc_sq(-93, 39.5, 0.4),
                                                 .sc_sq(-92, 39.5, 0.4), crs = 4326))
  lines <- .sc_log(ctr <- .lonlat_centre(x))
  ref <- .with_s2(sf::st_coordinates(sf::st_centroid(sf::st_union(sf::st_geometry(x)[2:3]))))
  expect_equal(unname(ctr[1, ]), unname(ref[1, 1:2]), tolerance = 1e-12)
  expect_true(any(grepl("s2 refuses 1 of the 3 features in longitude and latitude even repaired (row 1)",
                        lines, fixed = TRUE)))
  expect_true(any(grepl("taken over the other 2", lines, fixed = TRUE)))
  # A layer s2 refuses every feature of falls back to the vertex mean, and
  # says so.
  lines <- .sc_log(ctr <- .lonlat_centre(x[1, ]))
  expect_true(all(is.finite(ctr)))
  expect_true(any(grepl("s2 refuses every feature in longitude and latitude, even repaired",
                        lines, fixed = TRUE)))
  expect_true(any(grepl("mean of the vertices of the layer", lines, fixed = TRUE)))
  # An empty row holds no location: the counts are of the features that are
  # not empty, and add up.
  y <- sf::st_sf(id = 1:4, geometry = sf::st_sfc(.sc_bow(-94, 39, 0.4), sf::st_polygon(),
                                                 .sc_sq(-93, 39.5, 0.4), .sc_sq(-92, 39.5, 0.4),
                                                 crs = 4326))
  lines <- .sc_log(cy <- .lonlat_centre(y))
  expect_equal(unname(cy[1, ]), unname(ref[1, 1:2]), tolerance = 1e-12)
  expect_true(any(grepl(paste0("s2 refuses 1 of the 3 features that are not empty in longitude ",
                               "and latitude even repaired (row 1)"), lines, fixed = TRUE)))
  expect_true(any(grepl("taken over the other 2", lines, fixed = TRUE)))
})

test_that("with no union of the accepted features, the vertex mean leaves the refused ones out", {
  .sc_s2(TRUE)
  local_mocked_bindings(.split_crossing_edges = function(g1) NULL)
  x <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(.sc_bow(-94, 39, 0.4), .sc_sq(-93, 39.5, 0.4),
                                                 .sc_sq(-92, 39.5, 0.4), crs = 4326))
  local_mocked_bindings(st_union = function(x, y, ...) stop("no union here"), .package = "sf")
  lines <- .sc_log(ctr <- .lonlat_centre(x))
  # The mean of the unit vectors of the eight corners of the two squares.
  xy  <- expand.grid(x = c(-93, -92.6, -92, -91.6), y = c(39.5, 39.9))
  lam <- xy$x * pi / 180; phi <- xy$y * pi / 180
  v   <- c(sum(cos(phi) * cos(lam)), sum(cos(phi) * sin(lam)), sum(sin(phi)))
  ref <- c(atan2(v[2], v[1]), atan2(v[3], sqrt(v[1]^2 + v[2]^2))) * 180 / pi
  expect_equal(unname(ctr[1, ]), ref, tolerance = 1e-12)
  expect_true(any(grepl(paste0("s2 takes no union of the 2 features it accepts, so the centre ",
                               "that places the projection is the mean of the vertices of ",
                               "those features."), lines, fixed = TRUE)))
})

test_that("the vertex mean reads a layer of mixed types and counts closing vertices once", {
  .sc_s2(TRUE)
  local_mocked_bindings(.split_crossing_edges = function(g1) NULL)
  bow <- sf::st_polygon(list(rbind(c(-110, 30), c(-90, 45), c(-90, 30), c(-110, 45), c(-110, 30))))
  bow2 <- sf::st_polygon(list(rbind(c(-100, 32), c(-95, 40), c(-95, 32), c(-100, 40), c(-100, 32))))
  m  <- sf::st_sfc(bow, sf::st_multipolygon(list(bow2)), crs = 4326)
  mc <- sf::st_cast(m, "MULTIPOLYGON")
  expect_s3_class(m, "sfc_GEOMETRY")
  ca <- .lonlat_centre(m); cb <- .lonlat_centre(mc)
  expect_true(all(is.finite(ca)))
  expect_equal(ca, cb, tolerance = 1e-12)
  # The mean of the eight distinct corners, each once.
  xy  <- rbind(c(-110, 30), c(-90, 45), c(-90, 30), c(-110, 45),
               c(-100, 32), c(-95, 40), c(-95, 32), c(-100, 40)) * pi / 180
  v   <- colSums(cbind(cos(xy[, 2]) * cos(xy[, 1]), cos(xy[, 2]) * sin(xy[, 1]), sin(xy[, 2])))
  ref <- c(atan2(v[2], v[1]), atan2(v[3], sqrt(v[1]^2 + v[2]^2))) * 180 / pi
  expect_equal(unname(ca[1, ]), ref, tolerance = 1e-12)
})

test_that("layers s2 accepts get their centre and CRS as before", {
  .sc_s2(TRUE)
  for (x in list(.sc_hundred(), .sc_grid(5, 4), .sc_grid(5, 4, bad = 3L, maker = .sc_rep),
                 .sc_grid(5, 4, bad = 3L, maker = .sc_fig8))) {
    g <- sf::st_geometry(x)
    old <- tryCatch(sf::st_coordinates(sf::st_centroid(sf::st_union(g)))[1L, 1:2, drop = FALSE],
                    error = function(e)
                      sf::st_coordinates(sf::st_centroid(sf::st_union(.safe_make_valid(g))))[1L, 1:2, drop = FALSE])
    lines <- .sc_log(ctr <- .lonlat_centre(x))
    expect_identical(ctr, old)
    expect_false(any(grepl("centre that places the projection", lines, fixed = TRUE)))
  }
})

# ---- BUG-135: a fold probe saved without `kind` ----------------------------

.sc_folds <- function() {
  A <- .sc_grid(6, 5, x0 = -94, y0 = 39)
  A$v <- seq_len(nrow(A)); A$..row_id <- seq_len(nrow(A))
  f <- .sc_quiet(make_folds(A, k = 3, method = "random_kfold", seed = 1))
  legacy <- f
  legacy$params$row_probe$kind <- NULL
  list(A = A, folds = f, legacy = legacy)
}
.sc_other <- function(maker) {
  B <- .sc_grid(6, 5, bad = if (is.null(maker)) integer(0) else 7L,
                maker = if (is.null(maker)) .sc_rep else maker, x0 = -92, y0 = 39)
  B$..row_id <- seq_len(nrow(B))
  B
}

test_that("a legacy probe on data with a polygon s2 refuses still refuses other data's folds", {
  F <- .sc_folds()
  for (maker in list(.sc_rep, .sc_bow)) {
    B <- .sc_other(maker)
    .sc_s2(TRUE)
    expect_error(.check_fold_probe(F$legacy, B, "cv_spatial"),
                 "built from different data -- 30 of 30", fixed = TRUE)
    p <- .fold_row_probe(B, nrow(B), legacy = TRUE)
    expect_false(is.null(p))
    expect_identical(p$kind, "planar_centroid")
    # The same as the probe taken with s2 off, as a 2.0.0 probe was.
    q <- .fold_row_probe(B, nrow(B), legacy = FALSE)
    expect_identical(p[c("x", "y")], q[c("x", "y")])
    .sc_s2(FALSE)
    expect_error(.check_fold_probe(F$legacy, B, "cv_spatial"),
                 "built from different data -- 30 of 30", fixed = TRUE)
  }
})

test_that("legacy and current probes behave as before where s2 accepts the layer", {
  F <- .sc_folds()
  B <- .sc_other(NULL)
  .sc_s2(TRUE)
  p <- .fold_row_probe(B, nrow(B), legacy = TRUE)
  expect_identical(p$kind, "session_centroid")
  expect_error(.check_fold_probe(F$legacy, B, "cv_spatial"), "built from different data")
  expect_error(.check_fold_probe(F$folds, B, "cv_spatial"), "built from different data")
  expect_silent(.check_fold_probe(F$folds, F$A, "cv_spatial"))
  .sc_s2(FALSE)
  expect_silent(.check_fold_probe(F$legacy, F$A, "cv_spatial"))
  expect_identical(.fold_row_probe(B, nrow(B), legacy = TRUE)$kind, "session_centroid")
})

test_that("a legacy probe is retaken planar only for s2's refusal of a ring", {
  F <- .sc_folds()
  B <- .sc_other(NULL)
  .sc_s2(TRUE)
  # Another failure of the centroid is not s2's refusal of a polygon: the
  # check is skipped and says why, rather than taken on planar centroids.
  planar <- 0L
  local_mocked_bindings(.planar_centroid = function(g) {
    planar <<- planar + 1L
    suppressMessages(sf::sf_use_s2(FALSE)); on.exit(suppressMessages(sf::sf_use_s2(TRUE)))
    suppressWarnings(sf::st_centroid(g))
  })
  local_mocked_bindings(st_centroid = function(x, ...) stop("centroid failed for another reason"),
                        .package = "sf")
  expect_null(.fold_row_probe(B, nrow(B), legacy = TRUE))
  lines <- .sc_log(r <- .check_fold_probe(F$legacy, B, "cv_spatial"))
  expect_null(r)
  expect_identical(planar, 0L)
  expect_true(any(grepl("(\"centroid failed for another reason\"), so the supplied `folds`",
                        lines, fixed = TRUE)))
})

test_that("a provenance check that cannot be made is logged as skipped", {
  F <- .sc_folds()
  B <- .sc_other(NULL)
  .sc_s2(TRUE)
  local_mocked_bindings(.planar_centroid = function(g) stop("no centroid here"))
  expect_null(.fold_row_probe(B, nrow(B)))
  lines <- .sc_log(r <- .check_fold_probe(F$folds, B, "cv_spatial"))
  expect_null(r)
  expect_true(any(grepl(paste0("cv_spatial(): the locations of the checked rows of this data ",
                               "could not be taken (\"no centroid here\"), so the supplied ",
                               "`folds` cannot be compared with them; skipping the provenance ",
                               "check."), lines, fixed = TRUE)))
})

# ---- BUG-130: coerce_to_points(mode = "centroid") ---------------------------

.sc_layer <- function(maker = NULL) {
  x <- .sc_grid(5, 4, bad = if (is.null(maker)) integer(0) else 6L,
                maker = if (is.null(maker)) .sc_rep else maker)
  x$v <- seq_len(nrow(x)) / 7; x$x1 <- cos(seq_len(nrow(x)))
  x
}

test_that("a polygon with a repeated vertex gets the centroid of the same ring", {
  .sc_s2(TRUE)
  plain <- .sc_layer(); rep <- .sc_layer(.sc_rep)
  for (mode in c("centroid", "line_midpoint")) {
    a <- sf::st_coordinates(coerce_to_points(plain, mode))
    b <- sf::st_coordinates(coerce_to_points(rep, mode))
    expect_equal(b, a, tolerance = 1e-12)
    expect_identical(b[-6, ], a[-6, ])
  }
  pa <- prep_model_data(plain, "v", "x1", pointize = "centroid")
  pb <- prep_model_data(rep, "v", "x1", pointize = "centroid")
  expect_equal(sf::st_coordinates(pb), sf::st_coordinates(pa), tolerance = 1e-12)
})

test_that("a polygon s2 refuses even then stops the call with the features' class", {
  .sc_s2(TRUE)
  for (maker in list(.sc_bow, .sc_fig8)) {
    x <- .sc_layer(maker)
    for (mode in c("centroid", "line_midpoint")) {
      e <- tryCatch(coerce_to_points(x, mode), error = function(e) e)
      expect_s3_class(e, "spatialkit_features_degenerate")
      expect_identical(e$rows, 6L)
      expect_identical(e$n, 20L)
      expect_match(conditionMessage(e),
                   "^coerce_to_points\\(\\): s2, .* refuses 1 of the 20 features of `x` \\(row 6; ")
      expect_match(conditionMessage(e), "Transform `x` to a projected CRS", fixed = TRUE)
      expect_match(conditionMessage(e), "Or pass mode = \"auto\"", fixed = TRUE)
      expect_false(grepl("corners", conditionMessage(e)))
    }
    e <- tryCatch(prep_model_data(x, "v", "x1", pointize = "centroid"), error = function(e) e)
    expect_s3_class(e, "spatialkit_features_degenerate")
    expect_identical(e$rows, 6L)
    expect_match(conditionMessage(e),
                 "^prep_model_data\\(\\): s2, .* refuses 1 of the 20 features of `data_sf` \\(row 6; ")
    expect_match(conditionMessage(e), "Transform `data_sf` to a projected CRS", fixed = TRUE)
    expect_match(conditionMessage(e), "Or pass pointize = \"auto\"", fixed = TRUE)
    # The remedies run.  Moved to a projected CRS, the bowtie's lobes differ
    # in area and GEOS puts its centroid far outside it: that point is
    # refused, naming the row, where it was returned without a word (review
    # finding R-fb180-1 of batch found-bug-180); repaired, it runs.  The
    # figure eight's centroid lies in it.
    moved <- .sc_quiet(ensure_projected(x))
    if (identical(maker, .sc_bow)) {
      e <- tryCatch(prep_model_data(moved, "v", "x1", pointize = "centroid"),
                    error = function(e) e)
      expect_s3_class(e, "spatialkit_features_degenerate")
      expect_identical(e$site, "geos_centroid")
      expect_identical(e$rows, 6L)
      moved <- sf::st_make_valid(moved)
    }
    expect_no_error(prep_model_data(moved, "v", "x1", pointize = "centroid"))
    expect_no_error(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "auto")))
    expect_no_error(.sc_quiet(coerce_to_points(x, "auto")))
  }
})

test_that("a collection s2 refuses under mode auto names its row and point_on_surface", {
  .sc_s2(TRUE)
  x <- .sc_layer()
  sf::st_geometry(x)[3] <- sf::st_geometrycollection(list(.sc_bow(-93.8, 39)))
  e <- tryCatch(.sc_quiet(coerce_to_points(x, "auto")), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 3L)
  expect_match(conditionMessage(e), "Or pass mode = \"point_on_surface\"", fixed = TRUE)
  expect_no_error(.sc_quiet(coerce_to_points(x, "point_on_surface")))
})

test_that("a collection holding a polygon with a repeated vertex gets the centroid of the same ring", {
  .sc_s2(TRUE)
  x <- .sc_layer(); y <- .sc_layer()
  sf::st_geometry(x)[3] <- sf::st_geometrycollection(list(.sc_rep(-93.8, 39)))
  sf::st_geometry(y)[3] <- sf::st_geometrycollection(list(.sc_sq(-93.8, 39)))
  for (mode in c("auto", "centroid")) {
    a <- sf::st_coordinates(.sc_quiet(coerce_to_points(y, mode)))
    b <- sf::st_coordinates(.sc_quiet(coerce_to_points(x, mode)))
    expect_equal(b, a, tolerance = 1e-12)
    expect_identical(b[-3, ], a[-3, ])
  }
  pa <- .sc_quiet(prep_model_data(y, "v", "x1"))
  pb <- .sc_quiet(prep_model_data(x, "v", "x1"))
  expect_equal(sf::st_coordinates(pb), sf::st_coordinates(pa), tolerance = 1e-12)
  # The other members of a collection are kept as they were.
  z <- .sc_layer()
  sf::st_geometry(z)[3] <- sf::st_geometrycollection(list(
    .sc_rep(-93.8, 39), sf::st_point(c(-93.75, 39.05)),
    sf::st_linestring(rbind(c(-93.8, 39), c(-93.8, 39), c(-93.7, 39.1)))))
  d <- .drop_repeated_vertices(sf::st_geometry(z)[3])[[1]]
  expect_identical(unclass(d[[1]]), unclass(.sc_sq(-93.8, 39)))
  expect_identical(d[[2]], sf::st_geometry(z)[[3]][[2]])
  expect_identical(d[[3]], sf::st_geometry(z)[[3]][[3]])
})

test_that("a refusal quotes s2's error for the rows it names", {
  .sc_s2(TRUE)
  x <- .sc_grid(5, 4, bad = c(2L, 6L), maker = .sc_rep)
  sf::st_geometry(x)[6] <- .sc_bow(-93.5, 39)
  e <- tryCatch(coerce_to_points(x, "centroid"), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 6L)
  # Row 2, mended, is not the cause; its "duplicate vertex" is not quoted.
  expect_false(grepl("duplicate vertex", conditionMessage(e), fixed = TRUE))
  expect_false(grepl("duplicate vertex", conditionMessage(e$parent), fixed = TRUE))
})

test_that("layers s2 accepts, and any layer with s2 off, get their centroids as before", {
  .sc_s2(TRUE)
  x <- .sc_layer()
  g <- sf::st_geometry(x)
  expect_identical(sf::st_geometry(coerce_to_points(x, "centroid")),
                   suppressWarnings(sf::st_centroid(g)))
  .sc_s2(FALSE)
  for (maker in list(.sc_rep, .sc_bow, .sc_fig8)) {
    y <- .sc_layer(maker)
    expect_identical(sf::st_geometry(coerce_to_points(y, "centroid")),
                     suppressWarnings(sf::st_centroid(sf::st_geometry(y))))
  }
})

test_that("compare_models_cv() and the cv_*() functions name `data_sf` and its rows", {
  skip_if_not_installed("ranger")
  .sc_s2(TRUE)
  x <- .sc_layer(.sc_bow)
  e <- tryCatch(.sc_quiet(compare_models_cv(x, "v", "x1", models = "RF", k = 3,
                                            pointize = "centroid", quiet = TRUE)),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 6L)
  expect_match(conditionMessage(e),
               "^compare_models_cv\\(\\): s2, .* refuses 1 of the 20 features of `data_sf` \\(row 6; ")
  expect_match(conditionMessage(e), "Or pass pointize = \"auto\"", fixed = TRUE)
  e <- tryCatch(.sc_quiet(cv_rf(x, "v", "x1", k = 3, pointize = "centroid")),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 6L)
  expect_match(conditionMessage(e), "of `data_sf` (row 6; ", fixed = TRUE)
  # With a repeated vertex the call runs, on the centroid of the same ring.
  expect_no_error(.sc_quiet(cv_rf(.sc_layer(.sc_rep), "v", "x1", k = 3, pointize = "centroid")))
})

# ---- BUG-190: the remedy offered names a mode that runs on the rows refused --
# Under mode "centroid" (pointize = "centroid") the refusal offered mode
# "auto" for every row.  "auto" takes a point on the surface of a POLYGON or
# MULTIPOLYGON, but the centroid of a collection, so a collection row it was
# offered for was refused again with the same reason.

.sc_coll <- function(maker, rows = 3L) {
  x <- .sc_layer()
  for (r in rows) {
    bb <- sf::st_bbox(sf::st_geometry(x)[r])
    sf::st_geometry(x)[r] <- sf::st_geometrycollection(list(maker(bb[["xmin"]], bb[["ymin"]])))
  }
  x
}

test_that("a collection row refused under centroid offers point_on_surface, which runs", {
  .sc_s2(TRUE)
  for (maker in list(.sc_bow, .sc_fig8)) {
    x <- .sc_coll(maker)
    e <- tryCatch(coerce_to_points(x, "centroid"), error = function(e) e)
    expect_s3_class(e, "spatialkit_features_degenerate")
    expect_identical(e$rows, 3L)
    expect_identical(e$n, 20L)
    expect_false(e$auto_helps)
    expect_match(conditionMessage(e), paste0(
      "Or pass mode = \"point_on_surface\", which takes a point on each geometry ",
      "with GEOS.$"))
    expect_false(grepl("mode = \"auto\"", conditionMessage(e), fixed = TRUE))
    expect_match(conditionMessage(e), "Transform `x` to a projected CRS", fixed = TRUE)
    expect_no_error(.sc_quiet(coerce_to_points(x, "point_on_surface")))
    e <- tryCatch(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "centroid")),
                  error = function(e) e)
    expect_s3_class(e, "spatialkit_features_degenerate")
    expect_identical(e$rows, 3L)
    expect_false(e$auto_helps)
    expect_match(conditionMessage(e), paste0(
      "^prep_model_data\\(\\): s2, .* refuses 1 of the 20 features of `data_sf` ",
      "\\(row 3; "))
    expect_match(conditionMessage(e), paste0(
      "Or pass pointize = \"point_on_surface\", which takes a point on each geometry ",
      "with GEOS.$"))
    expect_false(grepl("pointize = \"auto\"", conditionMessage(e), fixed = TRUE))
    expect_no_error(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "point_on_surface")))
  }
})

test_that("compare_models_cv() offers pointize = \"point_on_surface\" for a collection row", {
  skip_if_not_installed("ranger")
  .sc_s2(TRUE)
  x <- .sc_coll(.sc_bow)
  e <- tryCatch(.sc_quiet(compare_models_cv(x, "v", "x1", models = "RF", k = 3,
                                            pointize = "centroid", quiet = TRUE)),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 3L)
  expect_match(conditionMessage(e),
               "^compare_models_cv\\(\\): s2, .* refuses 1 of the 20 features of `data_sf` \\(row 3; ")
  expect_match(conditionMessage(e), "Or pass pointize = \"point_on_surface\"", fixed = TRUE)
  expect_false(grepl("pointize = \"auto\"", conditionMessage(e), fixed = TRUE))
  expect_no_error(.sc_quiet(compare_models_cv(x, "v", "x1", models = "RF", k = 3,
                                              pointize = "point_on_surface", quiet = TRUE)))
})

test_that("POLYGON and MULTIPOLYGON rows refused under centroid still offer auto, which runs", {
  .sc_s2(TRUE)
  mbow <- function(x0, y0, d = 0.1) sf::st_multipolygon(list(unclass(.sc_bow(x0, y0, d))))
  for (maker in list(.sc_bow, .sc_fig8, mbow)) {
    x <- .sc_layer(maker)
    for (mode in c("centroid", "line_midpoint")) {
      e <- tryCatch(coerce_to_points(x, mode), error = function(e) e)
      expect_s3_class(e, "spatialkit_features_degenerate")
      expect_identical(e$rows, 6L)
      expect_true(e$auto_helps)
      expect_match(conditionMessage(e), paste0(
        "Or pass mode = \"auto\", which takes a point on the surface of each polygon ",
        "with GEOS.$"))
    }
    e <- tryCatch(prep_model_data(x, "v", "x1", pointize = "centroid"), error = function(e) e)
    expect_true(e$auto_helps)
    expect_match(conditionMessage(e), paste0(
      "Or pass pointize = \"auto\", which takes a point on the surface of each polygon ",
      "with GEOS.$"))
    expect_no_error(.sc_quiet(coerce_to_points(x, "auto")))
    expect_no_error(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "auto")))
  }
})

test_that("a refused POLYGON row beside a refused collection row offers point_on_surface", {
  .sc_s2(TRUE)
  x <- .sc_coll(.sc_bow)
  sf::st_geometry(x)[6] <- .sc_bow(-93.5, 39)
  e <- tryCatch(coerce_to_points(x, "centroid"), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, c(3L, 6L))
  expect_false(e$auto_helps)
  expect_match(conditionMessage(e), "Or pass mode = \"point_on_surface\"", fixed = TRUE)
  expect_false(grepl("mode = \"auto\"", conditionMessage(e), fixed = TRUE))
  e2 <- tryCatch(.sc_quiet(coerce_to_points(x, "auto")), error = function(e) e)
  expect_identical(e2$rows, 3L)
  e <- tryCatch(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "centroid")),
                error = function(e) e)
  expect_identical(e$rows, c(3L, 6L))
  expect_match(conditionMessage(e), "Or pass pointize = \"point_on_surface\"", fixed = TRUE)
  expect_false(grepl("pointize = \"auto\"", conditionMessage(e), fixed = TRUE))
  expect_no_error(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "point_on_surface")))
})

test_that("a collection row s2 accepts does not change the offer for a refused POLYGON row", {
  .sc_s2(TRUE)
  # Row 6 a bowtie POLYGON (refused), row 3 a collection holding a plain
  # square (accepted): only the refused rows decide the offer.
  x <- .sc_layer(.sc_bow)
  bb <- sf::st_bbox(sf::st_geometry(x)[3])
  sf::st_geometry(x)[3] <- sf::st_geometrycollection(list(.sc_sq(bb[["xmin"]], bb[["ymin"]])))
  e <- tryCatch(coerce_to_points(x, "centroid"), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 6L)
  expect_true(e$auto_helps)
  expect_match(conditionMessage(e), paste0(
    "Or pass mode = \"auto\", which takes a point on the surface of each polygon ",
    "with GEOS.$"))
  e <- tryCatch(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "centroid")),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 6L)
  expect_true(e$auto_helps)
  expect_match(conditionMessage(e), paste0(
    "Or pass pointize = \"auto\", which takes a point on the surface of each polygon ",
    "with GEOS.$"))
  expect_no_error(.sc_quiet(coerce_to_points(x, "auto")))
  expect_no_error(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "auto")))
})

test_that("the refusal under mode auto, and a condition without auto_helps, read as before", {
  .sc_s2(TRUE)
  x <- .sc_coll(.sc_bow)
  e <- tryCatch(.sc_quiet(coerce_to_points(x, "auto")), error = function(e) e)
  expect_false(e$auto_helps)
  expect_match(conditionMessage(e), paste0(
    "Or pass mode = \"point_on_surface\", which takes a point on each geometry ",
    "with GEOS.$"))
  e <- tryCatch(.sc_quiet(prep_model_data(x, "v", "x1", pointize = "auto")),
                error = function(e) e)
  expect_match(conditionMessage(e), paste0(
    "Or pass pointize = \"point_on_surface\", which takes a point on each geometry ",
    "with GEOS.$"))
  e <- tryCatch(.sc_quiet(prep_model_data(x, "v", "x1")), error = function(e) e)
  expect_match(conditionMessage(e), "Or pass pointize = \"point_on_surface\"", fixed = TRUE)
  # The callers that take no `pointize` keep their words.
  e <- tryCatch(.sc_quiet(make_folds(x, k = 3, method = "random")), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_match(conditionMessage(e), paste0(
    "Or reduce `points_sf` to points first, for example with ",
    "coerce_to_points(points_sf, \"point_on_surface\")"), fixed = TRUE)
  # A condition raised without the field (another caller) keeps the offer
  # of "auto" under another pointize.
  parent <- simpleError("Loop 0 is not valid: Edge 0 crosses edge 2")
  old <- errorCondition("x", class = "spatialkit_features_degenerate", rows = 2L, n = 5L,
                        reasons = "Loop 0 is not valid: Edge 0 crosses edge 2",
                        parent = parent)
  e <- tryCatch(.pointize_refused(old, "prep_model_data", "centroid"), error = function(e) e)
  expect_null(e$auto_helps)
  expect_match(conditionMessage(e), paste0(
    "Or pass pointize = \"auto\", which takes a point on the surface of each polygon ",
    "with GEOS.$"))
})
