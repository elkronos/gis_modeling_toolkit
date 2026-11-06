# BUG-194: a projected row far beyond the Earth is not moved to EPSG:4326.
#
# PROJ 9.5.1 inverts Web Mercator in time proportional to x: a row of an
# EPSG:3857 layer at x = 1e20 never came back from the move to EPSG:4326
# that the fold probe, the row fingerprint of compare_models(), the scored
# seeds, the area distortion check, harmonize_crs() and the shared moves of
# one layer into another's CRS (.align_crs(), .transform_or_stamp()) take.
# Such a row (a coordinate beyond 1e10 in the CRS's own unit) is now
# emptied before a move into a CRS in degrees and read as a row without a
# location.  Each test first installs a move that
# refuses such a row, so a regression fails here and does not hang.

# sf::st_transform(), refusing to move a projected set with a coordinate
# beyond 1e10 to a CRS in degrees.  Each refusal is counted in `seen`, since
# some callers catch the error.
.fpf_guard <- function(seen = new.env()) {
  real <- sf::st_transform
  function(x, crs, ...) {
    to <- if (inherits(crs, "crs")) crs else sf::st_crs(crs)
    if (isFALSE(sf::st_is_longlat(x)) && isTRUE(sf::st_is_longlat(to))) {
      # sf reads no coordinates of a set of several geometry types; every
      # number of its geometries is read then.
      xy <- tryCatch({
        xy <- suppressWarnings(sf::st_coordinates(x))
        xy[, intersect(colnames(xy), c("X", "Y")), drop = FALSE]
      }, error = function(e) unlist(sf::st_geometry(x)))
      if (any(!is.na(xy) & abs(xy) > 1e10)) {
        seen$refused <- (seen$refused %||% 0L) + 1L
        stop("moved a projected row beyond 1e10 to degrees")
      }
    }
    real(x, crs, ...)
  }
}

# 150 points of EPSG:32617 moved to EPSG:3857, row 1 put at (x, 5e6), or
# emptied for x = NULL, the others shifted by `dx`.
.fpf_layer <- function(x = NULL, dx = 0, n = 150L) {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  d <- sf::st_as_sf(data.frame(x = 5e5 + stats::runif(n, 0, 1e5) + dx,
                               y = 5e6 + stats::runif(n, 0, 1e5),
                               w = stats::rnorm(n)),
                    coords = c("x", "y"), crs = 32617)
  d$z <- 5 + 2 * d$w + stats::rnorm(n, 0, 0.5)
  d <- sf::st_transform(d, 3857)
  g <- sf::st_geometry(d)
  g[1] <- sf::st_sfc(if (is.null(x)) sf::st_point() else sf::st_point(c(x, 5e6)),
                     crs = 3857)
  sf::st_set_geometry(d, g)
}

.fpf_ids <- function(d) { d$..row_id <- seq_len(nrow(d)); d }

test_that("the fold probe gives NA for a projected row beyond 1e10 and leaves the others as they are", {
  local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
  ctl <- .fold_row_probe(.fpf_ids(.fpf_layer()), max_probe = 150L, catch = FALSE)
  expect_true(is.na(ctl$x[1L]) && is.na(ctl$y[1L]))
  for (x in c(1e20, -1e20, 1e12)) {
    pr <- .fold_row_probe(.fpf_ids(.fpf_layer(x)), max_probe = 150L, catch = FALSE)
    expect_identical(pr, ctl, label = format(x))
  }
})

test_that("make_folds() returns on a projected row beyond 1e10, with the folds the probe does not shape", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  far <- .fpf_layer(1e20)
  for (m in c("random_kfold", "block_kfold")) {
    ref <- local({
      local_mocked_bindings(.fold_row_probe = function(...) NULL)
      suppressWarnings(make_folds(far, k = 5, method = m, seed = 1))
    })
    f <- local({
      local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
      suppressWarnings(make_folds(far, k = 5, method = m, seed = 1))
    })
    expect_identical(f$folds, ref$folds, label = m)
    tested <- sort(unlist(lapply(f$folds, `[[`, "test")))
    expect_identical(as.integer(tested), seq_len(nrow(far)), label = m)
    pr <- f$params$row_probe
    expect_true(is.na(pr$x[pr$row_id == 1L]) && is.na(pr$y[pr$row_id == 1L]), label = m)
  }
})

test_that("supplied folds are checked on a layer with a projected row beyond 1e10", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
  far <- .fpf_layer(1e20)
  f <- suppressWarnings(make_folds(far, k = 5, method = "random_kfold", seed = 1))
  cv <- suppressWarnings(cv_spatial(far, "z", "w", fit_fn = lm_spatial_fit, folds = f))
  expect_true(is.list(cv))
  # Folds built on the same layer 1 km to the east are refused, as before.
  other <- suppressWarnings(make_folds(.fpf_layer(1e20, dx = 1000), k = 5,
                                       method = "random_kfold", seed = 1))
  expect_error(suppressWarnings(cv_spatial(far, "z", "w", fit_fn = lm_spatial_fit,
                                           folds = other)),
               "built from different data")
})

test_that("the row fingerprint gives NA for a projected row beyond 1e10", {
  local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
  ctl <- .fit_rowset(lm_spatial_fit(.fpf_layer(), "z", "w"))
  expect_false(is.null(ctl))
  for (x in c(1e20, 1e12)) {
    expect_identical(.fit_rowset(lm_spatial_fit(.fpf_layer(x), "z", "w")), ctl,
                     label = format(x))
  }
})

test_that("compare_models() compares fits whose rows differ only beyond 1e10 as before", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
  stub <- function(data, looic) {
    fit <- lm_spatial_fit(data, "z", "w")
    fit$info$looic <- looic
    class(fit) <- c("lmsurf_fit", "bayesian_fit", "spatial_fit")
    fit
  }
  quiet <- function(expr) {
    w <- character(0)
    v <- withCallingHandlers(expr, warning = function(c) {
      w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
    })
    list(value = v, warnings = w)
  }
  a <- stub(.fpf_layer(1e20), 54.8)
  # Row 1 at another far point: the same rows to the fingerprint, kept.
  b <- stub(.fpf_layer(3e20), 55.1)
  res <- quiet(compare_models(list(A = a, B = b)))
  expect_false(any(grepl("LOOIC is a sum", res$warnings)))
  expect_equal(res$value$LOOIC[match(c("A", "B"), res$value$model)], c(54.8, 55.1))
  # The other rows 1 km to the east: different rows, blanked with the warning.
  c2 <- stub(.fpf_layer(1e20, dx = 1000), 55.1)
  res2 <- quiet(compare_models(list(A = a, C = c2)))
  expect_true(any(grepl("LOOIC is a sum over the rows", res2$warnings)))
  expect_true(all(is.na(res2$value$LOOIC)))
})

test_that("the fold probe gives NA for a row at x = 1e16 without the guard", {
  # PROJ moves this row in about 0.16 s, so the real move is taken here.
  # The fingerprint the same.
  pr <- .fold_row_probe(.fpf_ids(.fpf_layer(1e16)), max_probe = 150L, catch = FALSE)
  ctl <- .fold_row_probe(.fpf_ids(.fpf_layer()), max_probe = 150L, catch = FALSE)
  expect_true(is.na(pr$x[1L]) && is.na(pr$y[1L]))
  expect_identical(pr, ctl)
  expect_identical(.fit_rowset(lm_spatial_fit(.fpf_layer(1e16), "z", "w")),
                   .fit_rowset(lm_spatial_fit(.fpf_layer(), "z", "w")))
})

test_that("the probe and the fingerprint of layers within the bound, in degrees or without a CRS are the plain move", {
  plain <- function(g) sf::st_transform(g, 4326)
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  nc$z <- nc$BIR74; nc$w <- nc$SID74
  nc32119 <- sf::st_transform(nc, 32119)
  near <- .fpf_layer(-8.7e6)
  ll <- sf::st_transform(near, 4326)
  crsless <- sf::st_set_crs(near, NA)
  ful <- local({
    set.seed(77L, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    x <- stats::runif(120, 0, 1e5); y <- stats::runif(120, 0, 1e5)
    sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, z = stats::rnorm(120),
                            w = stats::rnorm(120)), coords = c("x", "y"), crs = 32617)
  })
  kt <- local({
    xy <- withr::with_seed(4, cbind(stats::runif(120, 0, 5000), stats::runif(120, 0, 5000)))
    sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], z = stats::rnorm(120),
                            w = stats::rnorm(120)), coords = c("X", "Y"), crs = 32615)
  })
  layers <- list(nc32119 = nc32119, nc = nc, near = near, ll = ll, crsless = crsless,
                 ful3857 = sf::st_transform(ful, 3857), ful = ful, kt = kt)
  now_pr <- lapply(layers, function(d) .fold_row_probe(.fpf_ids(d), max_probe = nrow(d)))
  now_rs <- lapply(layers, function(d) .fit_rowset(lm_spatial_fit(d, "z", "w")))
  local_mocked_bindings(.lonlat_of_placed = plain)
  for (nm in names(layers)) {
    d <- layers[[nm]]
    expect_false(is.null(now_pr[[nm]]), label = nm)
    expect_identical(now_pr[[nm]], .fold_row_probe(.fpf_ids(d), max_probe = nrow(d)),
                     label = nm)
    expect_false(is.null(now_rs[[nm]]), label = nm)
    expect_identical(now_rs[[nm]], .fit_rowset(lm_spatial_fit(d, "z", "w")), label = nm)
  }
})

test_that("scored seeds moved to CRS-less degrees leave a projected seed beyond 1e10 empty", {
  local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
  like <- withr::with_seed(1, sf::st_as_sf(
    data.frame(x = stats::runif(50, -81, -79), y = stats::runif(50, 40, 41)),
    coords = c("x", "y")))
  seeds <- function(x) sf::st_as_sf(
    data.frame(seed_id = 1:3, x = c(x, -8.9e6, -8.8e6), y = c(5e6, 5.0e6, 5.1e6)),
    coords = c("x", "y"), crs = 3857)
  ref <- sf::st_set_crs(sf::st_transform(seeds(-8.7e6), 4326), NA)
  for (x in c(1e20, 1e12)) {
    s <- .seeds_in_crs_of(seeds(x), like)
    expect_true(is.na(sf::st_crs(s)), label = format(x))
    expect_true(sf::st_is_empty(s)[1L], label = format(x))
    expect_identical(sf::st_geometry(s)[2:3], sf::st_geometry(ref)[2:3], label = format(x))
  }
})

test_that("the area distortion is not measured on a projected polygon with a vertex beyond 1e10", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  # .crs_area_error() turns an error into NA, so the refusals are counted.
  seen <- new.env()
  local_mocked_bindings(st_transform = .fpf_guard(seen), .package = "sf")
  sq <- function(x0, y0, s, far = NULL) {
    m <- rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))
    if (!is.null(far)) m[2, ] <- c(far, y0)
    sf::st_polygon(list(m))
  }
  cells <- function(far = NULL) sf::st_sf(cell_id = 1:3, geometry = sf::st_sfc(
    sq(-8.9e6, 5.0e6, 1e4), sq(-8.8e6, 5.0e6, 1e4), sq(-8.7e6, 5.1e6, 1e4, far = far),
    crs = 3857))
  expect_true(is.finite(.crs_area_error(cells())))
  for (x in c(1e20, 1e12)) {
    expect_identical(.crs_area_error(cells(x)), NA_real_, label = format(x))
    expect_identical(.crs_area_error(cells(x), grid = TRUE), NA_real_, label = format(x))
  }
  lines <- capture_spatialkit_log(
    e <- ensure_projected(cells(1e20), purpose = "area"))
  expect_identical(sf::st_crs(e), sf::st_crs(3857))
  expect_true(log_has(lines, "could not be measured"))
  pts <- withr::with_seed(1, sf::st_as_sf(
    data.frame(x = -8.9e6 + stats::runif(30, 0, 1e4), y = 5e6 + stats::runif(30, 0, 1e4),
               v = 1, cell_id = c(rep(1L, 29), 3L)),
    coords = c("x", "y"), crs = 3857))
  expect_error(summarize_by_cell(pts, "v", cells_sf = cells(1e20), id_col = "cell_id",
                                 area = TRUE),
               "could not be measured")
  expect_null(seen$refused)
})

test_that("get_voronoi_seeds() returns a scored centre beyond 1e10 empty for CRS-less points in degrees", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  far <- .fpf_layer(1e20)
  rp <- suppressWarnings(resolution_profile(far, "z", n_levels = 3, seed = 1))
  ll <- sf::st_coordinates(sf::st_transform(far[-1, ], 4326))
  sp <- sf::st_as_sf(data.frame(x = ll[, 1], y = ll[, 2]), coords = c("x", "y"))
  local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
  s <- suppressWarnings(get_voronoi_seeds(n = rp, sample_points = sp))
  expect_true(is.na(sf::st_crs(s)))
  expect_identical(sum(sf::st_is_empty(s)), 1L)
})

test_that("two projected rows beyond 1e10 both give NA in the probe and the fingerprint", {
  local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
  two <- function(far) {
    d <- .fpf_layer(far[1])
    g <- sf::st_geometry(d)
    g[2] <- sf::st_sfc(if (is.null(far)) sf::st_point() else sf::st_point(c(far[2], 5e6)),
                       crs = 3857)
    sf::st_set_geometry(d, g)
  }
  ctl <- .fold_row_probe(.fpf_ids(two(NULL)), max_probe = 150L, catch = FALSE)
  pr  <- .fold_row_probe(.fpf_ids(two(c(1e20, -1e20))), max_probe = 150L, catch = FALSE)
  expect_true(all(is.na(pr$x[1:2])) && all(is.na(pr$y[1:2])))
  expect_identical(pr, ctl)
  expect_identical(.fit_rowset(lm_spatial_fit(two(c(1e20, -1e20)), "z", "w")),
                   .fit_rowset(lm_spatial_fit(two(NULL), "z", "w")))
})

test_that("a layer in degrees with a longitude beyond 1e10 takes the plain move", {
  d <- sf::st_transform(.fpf_layer(-8.7e6), 4326)
  g <- sf::st_geometry(d)
  g[1] <- sf::st_sfc(sf::st_point(c(1e12, 40)), crs = 4326)
  d <- .fpf_ids(sf::st_set_geometry(d, g))
  now <- .fold_row_probe(d, max_probe = 150L, catch = FALSE)
  expect_true(is.finite(now$x[1L]))
  local_mocked_bindings(.lonlat_of_placed = function(g, crs = 4326) sf::st_transform(g, crs))
  expect_identical(now, .fold_row_probe(d, max_probe = 150L, catch = FALSE))
})

test_that("the area distortion is not measured on a multipolygon with one part beyond 1e10", {
  seen <- new.env()
  local_mocked_bindings(st_transform = .fpf_guard(seen), .package = "sf")
  sq <- function(x0, y0, s) list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                                       c(x0, y0 + s), c(x0, y0)))
  mp <- sf::st_multipolygon(list(sq(-8.7e6, 5.1e6, 1e4), sq(1e20, 5.1e6, 1e4)))
  g <- sf::st_sfc(sf::st_polygon(sq(-8.9e6, 5e6, 1e4)), sf::st_polygon(sq(-8.8e6, 5e6, 1e4)),
                  mp, crs = 3857)
  expect_identical(.crs_area_error(g), NA_real_)
  expect_null(seen$refused)
  # One part beyond the bound is enough, in a polygon of a set or of a
  # collection, and the set moved to degrees has that row empty.
  expect_identical(.beyond_lonlat_move(g), c(FALSE, FALSE, TRUE))
  gc <- sf::st_sfc(sf::st_geometrycollection(list(sf::st_point(c(-8.7e6, 5e6)),
                                                  sf::st_linestring(rbind(c(-8.7e6, 5e6), c(1e20, 5e6))))),
                   sf::st_geometrycollection(list(sf::st_point(c(-8.7e6, 5e6)))), crs = 3857)
  expect_identical(.beyond_lonlat_move(gc), c(TRUE, FALSE))
  ll <- .align_crs(g, sf::st_sfc(sf::st_point(c(-80, 40)), crs = 4326))
  expect_identical(sf::st_is_empty(ll), c(FALSE, FALSE, TRUE))
  expect_identical(ll[1:2], sf::st_transform(g[1:2], 4326))
  expect_null(seen$refused)
})

test_that("scored seeds moved to a layer in EPSG:4326 leave a projected seed beyond 1e10 empty", {
  local_mocked_bindings(st_transform = .fpf_guard(), .package = "sf")
  like <- withr::with_seed(1, sf::st_as_sf(
    data.frame(x = stats::runif(50, -81, -79), y = stats::runif(50, 40, 41)),
    coords = c("x", "y"), crs = 4326))
  seeds <- function(x) sf::st_as_sf(
    data.frame(seed_id = 1:3, x = c(x, -8.9e6, -8.8e6), y = c(5e6, 5.0e6, 5.1e6)),
    coords = c("x", "y"), crs = 3857)
  ref <- sf::st_transform(seeds(-8.7e6), 4326)
  for (x in c(1e20, -1e20, 1e12)) {
    s <- .seeds_in_crs_of(seeds(x), like)
    expect_identical(sf::st_crs(s), sf::st_crs(4326), label = format(x))
    expect_true(sf::st_is_empty(s)[1L], label = format(x))
    expect_identical(sf::st_geometry(s)[2:3], sf::st_geometry(ref)[2:3], label = format(x))
  }
  # Seeds within the bound take the plain move.
  expect_identical(.seeds_in_crs_of(seeds(-8.7e6), like), ref)
})

test_that("get_voronoi_seeds() and voronoi_seeds_kmeans() return a scored centre beyond 1e10 empty for points in EPSG:4326", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  far <- .fpf_layer(1e20)
  rp <- suppressWarnings(resolution_profile(far, "z", n_levels = 3, seed = 1))
  sp <- sf::st_transform(far[-1, ], 4326)
  seen <- new.env()
  local_mocked_bindings(st_transform = .fpf_guard(seen), .package = "sf")
  s <- suppressWarnings(get_voronoi_seeds(n = rp, sample_points = sp))
  expect_identical(sf::st_crs(s), sf::st_crs(4326))
  expect_identical(sum(sf::st_is_empty(s)), 1L)
  k <- suppressWarnings(voronoi_seeds_kmeans(sp, k = rp))
  expect_identical(sf::st_crs(k), sf::st_crs(4326))
  expect_identical(sum(sf::st_is_empty(k)), 1L)
  expect_null(seen$refused)
})

test_that("a projected layer moved into a CRS in degrees leaves a row beyond 1e10 empty", {
  seen <- new.env()
  local_mocked_bindings(st_transform = .fpf_guard(seen), .package = "sf")
  ll  <- sf::st_transform(.fpf_layer(-8.7e6), 4326)
  ref <- sf::st_transform(.fpf_layer(NULL), 4326)
  for (x in c(1e20, 1e12)) {
    a <- .align_crs(.fpf_layer(x), ll)
    expect_identical(a, ref, label = format(x))
    expect_identical(suppressWarnings(.transform_or_stamp(.fpf_layer(x), sf::st_crs(4326))),
                     ref, label = format(x))
    h <- harmonize_crs(.fpf_layer(x), ll, prefer = "b")
    expect_identical(h$a, ref, label = format(x))
    h <- harmonize_crs(.fpf_layer(x), ll, target_crs = 4326)
    expect_identical(h$a, ref, label = format(x))
  }
  # Into a projected CRS the row is moved as before.
  utm <- sf::st_transform(ll, 32617)
  expect_identical(.align_crs(.fpf_layer(1e12), utm),
                   sf::st_transform(.fpf_layer(1e12), 32617))
  expect_null(seen$refused)
})

test_that("get_voronoi_seeds(method = \"kmeans\") drops a sample point beyond 1e10 moved into a boundary in degrees", {
  old_q <- spatialkit_quiet(TRUE); on.exit(spatialkit_quiet(old_q), add = TRUE)
  seen <- new.env()
  local_mocked_bindings(st_transform = .fpf_guard(seen), .package = "sf")
  b <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-81, 40), c(-79, 40), c(-79, 41), c(-81, 41), c(-81, 40)))), crs = 4326))
  sp <- withr::with_seed(1, sf::st_transform(sf::st_as_sf(
    data.frame(x = stats::runif(60, -81, -79), y = stats::runif(60, 40, 41)),
    coords = c("x", "y"), crs = 4326), 3857))
  g <- sf::st_geometry(sp)
  g[1] <- sf::st_sfc(sf::st_point(c(1e20, 5e6)), crs = 3857)
  sp <- sf::st_set_geometry(sp, g)
  expect_warning(s <- get_voronoi_seeds(n = 5, method = "kmeans", boundary = b,
                                        sample_points = sp, set_seed = 1),
                 "dropping 1 point")
  expect_identical(nrow(s), 5L)
  expect_identical(sf::st_crs(s), sf::st_crs(4326))
  expect_identical(attr(s, "kmeans")$rows, 2:60)
  expect_null(seen$refused)
})
