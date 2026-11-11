# tests/testthat/test-fix-mixed-zm-profile-folds.R
# ---------------------------------------------------------------------------
# sf reads a layer in the dimension of its first geometry: st_coordinates()
# of a point layer that mixes XY and XYZ points lays the values out of step
# (XYZ first) or stops on its "dimnames" (XY first), st_transform() of one
# with an XYZ geometry first stops with "OGR error", and GEOS stops on mixed
# polygons ("ParseException", "Points of LinearRing do not form a closed
# linestring").  rbind() of an XY and an XYZ layer builds such a layer
# without a word.  resolution_profile(), determine_optimal_levels(), the
# fold provenance probe, make_folds()'s `boundary`, `blocks` and
# `prediction_points`, fold_separation()'s `data_sf` and
# `prediction_points`, and prep_model_data()'s `boundary` are now read in x
# and y.  None of them returns the layer, so no warning is raised.  Each
# mixed layer is held to the same rows in XY, to the bit, warnings
# included.
# ---------------------------------------------------------------------------

.mpf_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value,
       classes = vapply(got, function(w) class(w)[1L], character(1)),
       texts = vapply(got, conditionMessage, character(1)))
}

# `got` equals `ref` in value, warning classes and warning texts, and raises
# no warning about Z and M.
.mpf_same <- function(got, ref, info = NULL) {
  expect_identical(got$value, ref$value, info = info)
  expect_identical(got$classes, ref$classes, info = info)
  expect_identical(got$texts, ref$texts, info = info)
  expect_false("spatialkit_zm_dropped" %in% got$classes, info = info)
}

# Feature `i` of `l` given a third coordinate at its own x and y: Z, or M
# (or both, for "XYZM").
.mpf_row <- function(l, i, dim = "XYZ") {
  sfc <- inherits(l, "sfc")
  g  <- if (sfc) l else sf::st_geometry(l)
  gi <- sf::st_zm(g[[i]], drop = TRUE, what = "ZM")
  extra <- nchar(dim) - 2L
  add <- function(m) cbind(m, matrix(3, nrow(m), extra))
  co <- unclass(gi)
  gi <- switch(class(gi)[2L],
    POINT = sf::st_point(c(co, rep(3, extra)), dim = dim),
    POLYGON = sf::st_polygon(lapply(co, add), dim = dim),
    MULTIPOLYGON = sf::st_multipolygon(lapply(co, function(p) lapply(p, add)), dim = dim))
  lst <- as.list(g); lst[[i]] <- gi
  out <- sf::st_sfc(lst, crs = sf::st_crs(g), recompute_bbox = TRUE)
  if (sfc) out else sf::st_set_geometry(l, out)
}
.mpf_allz <- function(l, dim = "XYZ") {
  for (i in seq_len(nrow(l))) l <- .mpf_row(l, i, dim)
  l
}

# Rows 1 and 5 (or 1 and 3 for the four quadrants), each given Z or M.
.mpf_cases <- function(rows = c(1L, 5L))
  expand.grid(row = rows, dim = c("XYZ", "XYM"), stringsAsFactors = FALSE)

.mpf_store <- new.env(parent = emptyenv())
.mpf_fixtures <- function() {
  if (is.null(.mpf_store$f)) {
    d <- local({
      withr::local_seed(6802, .rng_kind = "Mersenne-Twister",
                        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
      n <- 200
      x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
      z <- sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 0.3)
      data.frame(x = 5e5 + x, y = 5e6 + y, z = z, a1 = stats::rnorm(n),
                 a2 = 0.5 * z + stats::rnorm(n, 0, 0.3))
    })
    pts   <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
    polys <- sf::st_set_geometry(pts, sf::st_buffer(sf::st_geometry(pts), 10,
                                                    nQuadSegs = 2))
    bb <- sf::st_bbox(pts)
    mx <- (bb[["xmin"]] + bb[["xmax"]]) / 2; my <- (bb[["ymin"]] + bb[["ymax"]]) / 2
    q <- function(x0, y0, x1, y1)
      sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
    bnd4 <- sf::st_sf(id = 1:4, geometry = sf::st_sfc(
      q(bb[["xmin"]], bb[["ymin"]], mx, my), q(mx, bb[["ymin"]], bb[["xmax"]], my),
      q(bb[["xmin"]], my, mx, bb[["ymax"]]), q(mx, my, bb[["xmax"]], bb[["ymax"]]),
      crs = 32617))
    pp  <- pts[seq(1, nrow(pts), by = 4), "a1"]
    ppp <- sf::st_set_geometry(pp, sf::st_buffer(sf::st_geometry(pp), 10, nQuadSegs = 2))
    .mpf_store$f <- list(pts = pts, pts_ll = sf::st_transform(pts, 4326), polys = polys,
                         bnd4 = bnd4, pp = pp, ppp = ppp)
  }
  .mpf_store$f
}

test_that("the premise: sf reads a mixed layer in the dimension of its first geometry", {
  F <- .mpf_fixtures()
  utils::capture.output(expect_error(sf::st_transform(.mpf_row(F$pts, 1), 32616)))
  expect_error(suppressWarnings(sf::st_coordinates(.mpf_row(F$pts, 5))))
  # sf raises GEOS's error from an on.exit() handler, which expect_error()
  # does not see; tryCatch() does.
  expect_s3_class(tryCatch(sf::st_union(.mpf_row(F$bnd4, 3)), error = function(e) e),
                  "error")
})

# ---- BUG-68 ---------------------------------------------------------------

test_that("resolution_profile() reads mixed points as the same points in XY", {
  F <- .mpf_fixtures()
  rp <- function(l) .mpf_catch(resolution_profile(l, "z", n_levels = 4L, seed = 1))
  for (lay in c("pts", "pts_ll")) {
    ref <- rp(F[[lay]])
    cs <- if (lay == "pts") .mpf_cases() else .mpf_cases()[1:2, ]
    for (j in seq_len(nrow(cs)))
      .mpf_same(rp(.mpf_row(F[[lay]], cs$row[j], cs$dim[j])), ref,
                info = paste(lay, cs$row[j], cs$dim[j]))
  }
})

test_that("resolution_profile() reads mixed polygons with no warning about Z", {
  F <- .mpf_fixtures()
  rp <- function(l) .mpf_catch(resolution_profile(l, "z", n_levels = 4L, seed = 1))
  ref <- rp(F$polys)
  for (i in c(1L, 5L)) .mpf_same(rp(.mpf_row(F$polys, i)), ref, info = paste(i))
})

test_that("determine_optimal_levels() reads mixed points and polygons as the same layer in XY", {
  F <- .mpf_fixtures()
  dl <- function(l) .mpf_catch(determine_optimal_levels(l, max_levels = 6L))
  for (lay in c("pts", "pts_ll")) {
    ref <- dl(F[[lay]])
    cs <- if (lay == "pts") .mpf_cases() else .mpf_cases()[1:2, ]
    for (j in seq_len(nrow(cs)))
      .mpf_same(dl(.mpf_row(F[[lay]], cs$row[j], cs$dim[j])), ref,
                info = paste(lay, cs$row[j], cs$dim[j]))
  }
  ref <- dl(F$polys)
  for (i in c(1L, 5L)) .mpf_same(dl(.mpf_row(F$polys, i)), ref, info = paste("polys", i))
  # With a response, the model-aware pass reads the same rows.
  dm <- function(l) .mpf_catch(determine_optimal_levels(l, max_levels = 6L,
                                                        response_var = "z",
                                                        criterion = "combined"))
  ref <- dm(F$pts)
  .mpf_same(dm(.mpf_row(F$pts, 5)), ref, info = "combined")
})

test_that("make_folds() reads mixed points in degrees as the same points in XY", {
  F <- .mpf_fixtures()
  mk <- function(l) .mpf_catch(make_folds(l, k = 5, method = "random_kfold", seed = 1))
  ref <- mk(F$pts_ll)
  for (i in c(1L, 5L)) .mpf_same(mk(.mpf_row(F$pts_ll, i)), ref, info = paste(i))
})

# ---- BUG-70 ---------------------------------------------------------------

test_that("the provenance probe reads mixed rows as the same rows in XY", {
  F <- .mpf_fixtures()
  for (crs in c(TRUE, FALSE)) {
    p <- F$pts
    if (!crs) sf::st_crs(p) <- NA
    p$..row_id <- seq_len(nrow(p))
    ref <- .fold_row_probe(p)
    expect_false(is.null(ref))
    cs <- .mpf_cases()
    for (j in seq_len(nrow(cs))) {
      got <- .fold_row_probe(.mpf_row(p, cs$row[j], cs$dim[j]))
      expect_identical(got, ref, info = paste(crs, cs$row[j], cs$dim[j]))
      expect_identical(.fold_row_probe(.mpf_row(p, cs$row[j], cs$dim[j]),
                                       max_probe = nrow(p)),
                       .fold_row_probe(p, max_probe = nrow(p)),
                       info = paste("all rows", crs, cs$row[j], cs$dim[j]))
    }
  }
  # Polygons: the probe takes their centroids in XY.
  p <- F$polys; p$..row_id <- seq_len(nrow(p))
  for (i in c(1L, 5L))
    expect_identical(.fold_row_probe(.mpf_row(p, i)), .fold_row_probe(p), info = paste(i))
  # make_folds() records the XY layer's probe on a mixed layer.
  ref <- make_folds(F$pts, k = 5, method = "random_kfold", seed = 1)
  got <- make_folds(.mpf_row(F$pts, 1), k = 5, method = "random_kfold", seed = 1)
  expect_false(is.null(got$params$row_probe))
  expect_identical(got$params$row_probe, ref$params$row_probe)
})

test_that("folds built on other data are refused on mixed rows, whichever row comes first", {
  F <- .mpf_fixtures()
  folds <- make_folds(F$pts, k = 5, method = "random_kfold", seed = 1)
  other <- F$pts
  sf::st_geometry(other) <- sf::st_geometry(F$pts) + c(5e4, 0)
  sf::st_crs(other) <- 32617
  other$poly_id <- 1L + (seq_len(nrow(other)) %% 4L)
  cells <- F$bnd4; cells$poly_id <- 1:4
  msg <- "built from different data"
  expect_error(fold_separation(folds, other), msg)
  cs <- .mpf_cases()
  for (j in seq_len(nrow(cs))) {
    l <- .mpf_row(other, cs$row[j], cs$dim[j])
    expect_error(suppressWarnings(fold_separation(folds, l)), msg,
                 info = paste(cs$row[j], cs$dim[j]))
    # kriging_adequacy() needs gstat (Suggests) before it reads the folds.
    if (requireNamespace("gstat", quietly = TRUE))
      expect_error(suppressWarnings(kriging_adequacy(l, "z", cells, folds = folds)),
                   msg, info = paste(cs$row[j], cs$dim[j]))
  }
  skip_if_not_installed("ranger")
  for (j in seq_len(nrow(cs)))
    expect_error(suppressWarnings(cv_rf(.mpf_row(other, cs$row[j], cs$dim[j]), "z",
                                        c("a1", "a2"), folds = folds, num_trees = 20,
                                        seed = 1)),
                 msg, info = paste(cs$row[j], cs$dim[j]))
})

test_that("folds built on the same rows without a CRS are accepted on mixed rows", {
  skip_if_not_installed("ranger")
  F <- .mpf_fixtures()
  p0 <- F$pts; sf::st_crs(p0) <- NA
  folds <- make_folds(p0, k = 5, method = "random_kfold", seed = 1)
  cv <- function(l) .mpf_catch(cv_rf(l, "z", c("a1", "a2"), folds = folds,
                                     num_trees = 20, seed = 1))
  ref <- cv(p0)
  for (i in c(1L, 5L)) {
    got <- cv(.mpf_row(p0, i))
    expect_identical(got$value$predictions, ref$value$predictions, info = paste(i))
    expect_identical(got$value$overall, ref$value$overall, info = paste(i))
    expect_identical(got$classes, ref$classes, info = paste(i))
    expect_false("spatialkit_zm_dropped" %in% got$classes, info = paste(i))
  }
})

# ---- BUG-71 ---------------------------------------------------------------

test_that("make_folds() reads a mixed boundary or mixed blocks as the same layer in XY", {
  F <- .mpf_fixtures()
  cs <- .mpf_cases(c(1L, 3L))
  bk <- function(p, b) .mpf_catch(make_folds(p, k = 5, method = "block_kfold", seed = 1,
                                             block_size = 20000, boundary = b))
  for (lay in c("pts", "pts_ll")) {
    ref <- bk(F[[lay]], F$bnd4)
    for (j in seq_len(nrow(cs)))
      .mpf_same(bk(F[[lay]], .mpf_row(F$bnd4, cs$row[j], cs$dim[j])), ref,
                info = paste("boundary", lay, cs$row[j], cs$dim[j]))
  }
  # An sfc boundary, and a boundary without a CRS around points without one.
  ref <- bk(F$pts, sf::st_geometry(F$bnd4))
  .mpf_same(bk(F$pts, .mpf_row(sf::st_geometry(F$bnd4), 1)), ref, info = "sfc")
  p0 <- F$pts; sf::st_crs(p0) <- NA
  b0 <- F$bnd4; sf::st_crs(b0) <- NA
  ref <- bk(p0, b0)
  for (i in c(1L, 3L)) .mpf_same(bk(p0, .mpf_row(b0, i)), ref, info = paste("CRS-less", i))
  bl <- function(p, b) .mpf_catch(make_folds(p, k = 4, method = "block_kfold", seed = 1,
                                             blocks = b))
  for (lay in c("pts", "pts_ll")) {
    ref <- bl(F[[lay]], F$bnd4)
    for (j in seq_len(nrow(cs)))
      .mpf_same(bl(F[[lay]], .mpf_row(F$bnd4, cs$row[j], cs$dim[j])), ref,
                info = paste("blocks", lay, cs$row[j], cs$dim[j]))
  }
})

test_that("the cv_*() wrappers pass a mixed boundary on in x and y", {
  F <- .mpf_fixtures()
  expect_identical(.cv_boundary_crs(.mpf_row(F$bnd4, 1), "cv_rf"), F$bnd4)
  b0 <- F$bnd4; sf::st_crs(b0) <- NA
  p0 <- F$pts; sf::st_crs(p0) <- NA
  expect_identical(suppressWarnings(.cv_boundary_crs(.mpf_row(b0, 3), "cv_rf", to = F$pts)),
                   suppressWarnings(.cv_boundary_crs(b0, "cv_rf", to = F$pts)))
  skip_if_not_installed("ranger")
  cv <- function(p, b) .mpf_catch(cv_rf(p, "z", c("a1", "a2"), k = 5, seed = 1, boundary = b,
                                        block_size = 20000, num_trees = 20))
  for (lay in list(list(F$pts, F$bnd4, "projected"), list(p0, b0, "CRS-less"))) {
    ref <- cv(lay[[1]], lay[[2]])
    for (i in c(1L, 3L)) {
      got <- cv(lay[[1]], .mpf_row(lay[[2]], i))
      expect_identical(got$value$predictions, ref$value$predictions, info = paste(lay[[3]], i))
      expect_identical(got$value$overall, ref$value$overall, info = paste(lay[[3]], i))
      expect_identical(got$classes, ref$classes, info = paste(lay[[3]], i))
      expect_identical(got$texts, ref$texts, info = paste(lay[[3]], i))
      expect_false("spatialkit_zm_dropped" %in% got$classes, info = paste(lay[[3]], i))
    }
  }
})

# ---- BUG-72 ---------------------------------------------------------------

test_that("make_folds(method = 'nndm' | 'knndm') reads mixed prediction points as in XY", {
  F <- .mpf_fixtures()
  nn <- function(l) .mpf_catch(make_folds(F$pts, k = 5, method = "nndm", seed = 1,
                                          prediction_points = l, max_test = 30))
  ref <- nn(F$pp)
  cs <- .mpf_cases()
  for (j in seq_len(nrow(cs)))
    .mpf_same(nn(.mpf_row(F$pp, cs$row[j], cs$dim[j])), ref, info = paste(cs$row[j], cs$dim[j]))
  pll <- sf::st_transform(F$pp, 4326)
  ref <- nn(pll)
  for (i in c(1L, 5L)) .mpf_same(nn(.mpf_row(pll, i)), ref, info = paste("degrees", i))
  ref <- nn(sf::st_geometry(F$pp))
  .mpf_same(nn(.mpf_row(sf::st_geometry(F$pp), 1)), ref, info = "sfc")
  ref <- nn(F$ppp)
  for (i in c(1L, 5L)) .mpf_same(nn(.mpf_row(F$ppp, i)), ref, info = paste("polygons", i))
  kn <- function(l) .mpf_catch(make_folds(F$pts, k = 5, method = "knndm", seed = 1,
                                          prediction_points = l))
  ref <- kn(F$pp)
  for (i in c(1L, 5L)) .mpf_same(kn(.mpf_row(F$pp, i)), ref, info = paste("knndm", i))
})

# ---- BUG-73 ---------------------------------------------------------------

test_that("fold_separation() reads mixed data_sf and prediction points as in XY", {
  F <- .mpf_fixtures()
  folds <- make_folds(F$pts, k = 5, method = "random_kfold", seed = 1)
  fs <- function(...) .mpf_catch(fold_separation(folds, ...))
  ref <- fs(F$pts, prediction_points = F$pp)
  cs <- .mpf_cases()
  for (j in seq_len(nrow(cs)))
    .mpf_same(fs(F$pts, prediction_points = .mpf_row(F$pp, cs$row[j], cs$dim[j])), ref,
              info = paste("pp", cs$row[j], cs$dim[j]))
  ref <- fs(F$pts, prediction_points = F$ppp)
  for (i in c(1L, 5L))
    .mpf_same(fs(F$pts, prediction_points = .mpf_row(F$ppp, i)), ref,
              info = paste("pp polygons", i))
  sac <- structure(15000, crs = sf::st_crs(32617))
  p16 <- sf::st_transform(F$pts, 32616)
  ref <- fs(p16, sac)
  for (j in seq_len(nrow(cs)))
    .mpf_same(fs(.mpf_row(p16, cs$row[j], cs$dim[j]), sac), ref,
              info = paste("32616", cs$row[j], cs$dim[j]))
  ref <- fs(F$pts_ll)
  for (i in c(1L, 5L)) .mpf_same(fs(.mpf_row(F$pts_ll, i)), ref, info = paste("degrees", i))
  ref <- fs(F$polys)
  for (i in c(1L, 5L)) .mpf_same(fs(.mpf_row(F$polys, i)), ref, info = paste("polygons", i))
})

# ---- found in the build: prep_model_data()'s boundary ----------------------

test_that("prep_model_data() reads a mixed boundary in degrees with no warning about Z", {
  F <- .mpf_fixtures()
  bll <- sf::st_transform(F$bnd4, 4326)
  pr <- function(b) .mpf_catch(prep_model_data(F$pts, "z", c("a1", "a2"), boundary = b))
  ref <- pr(bll)
  for (i in c(1L, 3L)) .mpf_same(pr(.mpf_row(bll, i)), ref, info = paste(i))
})

# ---- a layer XYZ throughout whose points come out mixed --------------------

# The first half of `p`'s rows as points and the rest as `q`'s polygons (or
# the reverse), XYZ throughout or XY throughout.  The layer XYZ throughout
# is in one dimension, but the point on the surface of a POLYGON Z is XY
# beside the POINT Z rows kept, so the points it is reduced to mix them.
.mpf_pts_polys <- function(p, q, polys_first = FALSE, z = TRUE) {
  h <- seq_len(nrow(p) %/% 2L); t <- setdiff(seq_len(nrow(p)), h)
  if (z) { p <- .mpf_allz(p); q <- .mpf_allz(q) }
  if (polys_first) rbind(q[h, ], p[t, ]) else rbind(p[h, ], q[t, ])
}

test_that("a layer XYZ throughout that mixes points and polygons is read as in XY", {
  F <- .mpf_fixtures()
  folds <- make_folds(F$pts, k = 5, method = "random_kfold", seed = 1)
  fs <- function(...) .mpf_catch(fold_separation(folds, ...))
  rp <- function(l) .mpf_catch(resolution_profile(l, "z", n_levels = 4L, seed = 1))
  dl <- function(l) .mpf_catch(determine_optimal_levels(l, max_levels = 6L))
  nn <- function(l) .mpf_catch(make_folds(F$pts, k = 5, method = "nndm", seed = 1,
                                          prediction_points = l, max_test = 30))
  kn <- function(l) .mpf_catch(make_folds(F$pts, k = 5, method = "knndm", seed = 1,
                                          prediction_points = l))
  sac16 <- structure(15000, crs = sf::st_crs(32616))
  for (pf in c(FALSE, TRUE)) {
    info <- if (pf) "polygons first" else "points first"
    dz <- .mpf_pts_polys(F$pts, F$polys, pf)
    dx <- .mpf_pts_polys(F$pts, F$polys, pf, z = FALSE)
    # The premise: the layer is in one dimension and its points are not.
    expect_false(.mixes_zm(dz), info = info)
    expect_true(.mixes_zm(suppressWarnings(coerce_to_points(dz, "auto"))), info = info)
    .mpf_same(rp(dz), rp(dx), info = paste("resolution_profile()", info))
    .mpf_same(dl(dz), dl(dx), info = paste("determine_optimal_levels()", info))
    .mpf_same(fs(dz, sac16), fs(dx, sac16), info = paste("fold_separation() data, sac", info))
    .mpf_same(fs(sf::st_transform(dz, 4326)), fs(sf::st_transform(dx, 4326)),
              info = paste("fold_separation() data in degrees", info))
    pz <- .mpf_pts_polys(F$pp, F$ppp, pf)
    px <- .mpf_pts_polys(F$pp, F$ppp, pf, z = FALSE)
    .mpf_same(fs(F$pts, prediction_points = pz), fs(F$pts, prediction_points = px),
              info = paste("fold_separation() prediction points", info))
    .mpf_same(nn(pz), nn(px), info = paste("nndm", info))
    .mpf_same(kn(pz), kn(px), info = paste("knndm", info))
  }
})

# ---- layers whose geometries carry M ---------------------------------------

test_that("layers whose geometries carry M are read as the same layers in XY", {
  F <- .mpf_fixtures()
  pm <- .mpf_allz(F$polys, "XYM")
  bm <- .mpf_allz(F$bnd4, "XYM")
  ppm <- .mpf_allz(F$ppp, "XYM")
  rp <- function(l) .mpf_catch(resolution_profile(l, "z", n_levels = 4L, seed = 1))
  .mpf_same(rp(.mpf_allz(F$pts, "XYM")), rp(F$pts), info = "profile, points XYM")
  .mpf_same(rp(.mpf_allz(F$pts, "XYZM")), rp(F$pts), info = "profile, points XYZM")
  .mpf_same(rp(pm), rp(F$polys), info = "profile, polygons XYM")
  dl <- function(l) .mpf_catch(determine_optimal_levels(l, max_levels = 6L))
  .mpf_same(dl(.mpf_allz(F$pts, "XYM")), dl(F$pts), info = "levels, points XYM")
  .mpf_same(dl(pm), dl(F$polys), info = "levels, polygons XYM")
  # make_folds(): the layer, the boundary, the blocks and the prediction points.
  mk <- function(l) .mpf_catch(make_folds(l, k = 5, method = "random_kfold", seed = 1))
  .mpf_same(mk(pm), mk(F$polys), info = "make_folds(), polygons XYM")
  .mpf_same(mk(.mpf_allz(F$pts, "XYM")), mk(F$pts), info = "make_folds(), points XYM")
  bk <- function(b) .mpf_catch(make_folds(F$pts, k = 5, method = "block_kfold", seed = 1,
                                          block_size = 20000, boundary = b))
  .mpf_same(bk(bm), bk(F$bnd4), info = "boundary XYM")
  .mpf_same(bk(sf::st_geometry(bm)), bk(sf::st_geometry(F$bnd4)), info = "sfc boundary XYM")
  bl <- function(b) .mpf_catch(make_folds(F$pts, k = 4, method = "block_kfold", seed = 1,
                                          blocks = b))
  .mpf_same(bl(bm), bl(F$bnd4), info = "blocks XYM")
  nn <- function(l) .mpf_catch(make_folds(F$pts, k = 5, method = "nndm", seed = 1,
                                          prediction_points = l, max_test = 30))
  .mpf_same(nn(ppm), nn(F$ppp), info = "nndm, polygons XYM")
  kn <- function(l) .mpf_catch(make_folds(F$pts, k = 5, method = "knndm", seed = 1,
                                          prediction_points = l))
  .mpf_same(kn(ppm), kn(F$ppp), info = "knndm, polygons XYM")
  # fold_separation(): the layer and the prediction points.
  fp <- make_folds(F$polys, k = 5, method = "random_kfold", seed = 1)
  .mpf_same(.mpf_catch(fold_separation(fp, pm)), .mpf_catch(fold_separation(fp, F$polys)),
            info = "fold_separation(), polygons XYM")
  folds <- make_folds(F$pts, k = 5, method = "random_kfold", seed = 1)
  .mpf_same(.mpf_catch(fold_separation(folds, F$pts, prediction_points = ppm)),
            .mpf_catch(fold_separation(folds, F$pts, prediction_points = F$ppp)),
            info = "fold_separation(), prediction polygons XYM")
  # The provenance probe takes its sample of polygons XYM, and the cv_*()
  # wrappers pass a boundary XYM on in XY.
  p <- F$polys; p$..row_id <- seq_len(nrow(p))
  q <- pm; q$..row_id <- seq_len(nrow(q))
  expect_false(is.null(.fold_row_probe(p)))
  expect_identical(.fold_row_probe(q), .fold_row_probe(p))
  expect_identical(.cv_boundary_crs(bm, "cv_rf"), F$bnd4)
  # prep_model_data() reads polygons XYM, as the cv_*() wrappers pass them.
  .mpf_same(.mpf_catch(prep_model_data(pm, "z", c("a1", "a2"))),
            .mpf_catch(prep_model_data(F$polys, "z", c("a1", "a2"))),
            info = "prep_model_data(), polygons XYM")
  skip_if_not_installed("gstat")
  .mpf_same(.mpf_catch(estimate_sac_range(pm, "z")),
            .mpf_catch(estimate_sac_range(F$polys, "z")),
            info = "estimate_sac_range(), polygons XYM")
  skip_if_not_installed("ranger")
  cv <- function(l, ...) .mpf_catch(cv_rf(l, "z", c("a1", "a2"), num_trees = 20, seed = 1, ...))
  got <- cv(pm, folds = fp); ref <- cv(F$polys, folds = fp)
  expect_identical(got$value$predictions, ref$value$predictions)
  expect_identical(got$classes, ref$classes)
  expect_false("spatialkit_zm_dropped" %in% got$classes)
  got <- cv(F$pts, k = 5, boundary = bm, block_size = 20000)
  ref <- cv(F$pts, k = 5, boundary = F$bnd4, block_size = 20000)
  expect_identical(got$value$predictions, ref$value$predictions)
  expect_identical(got$classes, ref$classes)
  expect_false("spatialkit_zm_dropped" %in% got$classes)
})

# ---- the old path ---------------------------------------------------------

test_that("layers in one dimension take the path they took", {
  F <- .mpf_fixtures()
  # A layer XYZ throughout is not flattened: the profile's hull keeps its Z.
  rz <- suppressWarnings(resolution_profile(.mpf_allz(F$pts), "z", n_levels = 4L, seed = 1))
  rx <- suppressWarnings(resolution_profile(F$pts, "z", n_levels = 4L, seed = 1))
  expect_identical(class(attr(rz, "hull")[[1]])[1], "XYZ")
  attr(rz, "hull") <- attr(rx, "hull")
  expect_identical(rz, rx)
  expect_identical(.mixes_zm(.mpf_allz(F$bnd4)), FALSE)
  # make_folds() reads boundary, blocks and prediction points XYZ throughout
  # as in XY, as it did.
  expect_identical(
    make_folds(F$pts, k = 5, method = "block_kfold", seed = 1, block_size = 20000,
               boundary = .mpf_allz(F$bnd4)),
    make_folds(F$pts, k = 5, method = "block_kfold", seed = 1, block_size = 20000,
               boundary = F$bnd4))
  expect_identical(
    suppressWarnings(make_folds(F$pts, k = 5, method = "nndm", seed = 1,
                                prediction_points = .mpf_allz(F$pp), max_test = 30)),
    suppressWarnings(make_folds(F$pts, k = 5, method = "nndm", seed = 1,
                                prediction_points = F$pp, max_test = 30)))
  # .xy_if_mixed() hands a layer in one dimension back as it came.
  for (l in list(F$pts, .mpf_allz(F$pp), F$bnd4, .mpf_allz(F$bnd4)))
    expect_identical(.xy_if_mixed(l), l)
})

test_that("the help says these layers are read in x and y", {
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  pts_text <- paste0("A layer that mixes XY and XYZ (or XYM) geometries is read in x and y, ",
                     "before and after it is reduced to points, since sf reads a layer in ",
                     "the dimension of its first geometry.")
  expect_match(page("resolution_profile.Rd"), pts_text, fixed = TRUE)
  expect_match(page("determine_optimal_levels.Rd"), pts_text, fixed = TRUE)
  mf <- page("make_folds.Rd")
  expect_match(mf, pts_text, fixed = TRUE)
  expect_match(mf, paste0("A boundary that mixes XY and XYZ (or XYM) geometries is read in x ",
                          "and y"), fixed = TRUE)
  expect_match(mf, paste0("Blocks that mix XY and XYZ (or XYM) polygons are read in x and y, ",
                          "since sf reads a layer in the dimension of its first geometry, and ",
                          "params$blocks holds them in XY."), fixed = TRUE)
  fs <- page("fold_separation.Rd")
  expect_match(fs, paste0("A layer that mixes XY and XYZ (or XYM) geometries is read in x and ",
                          "y, for the distances and for that check"), fixed = TRUE)
  expect_match(fs, paste0("one that mixes XY and XYZ (or XYM) geometries is read in x and y ",
                          "first."), fixed = TRUE)
  expect_match(page("prep_model_data.Rd"),
               "One that mixes XY and XYZ (or XYM) geometries is read in x and y",
               fixed = TRUE)
  # Layers that carry M, and the points a layer is reduced to.
  m_text <- "So is a layer whose geometries carry M (XYM or XYZM), which GEOS refuses."
  red_text <- paste0("The points a layer is reduced to are read in x and y when they mix ",
                     "dimensions, as those of a layer XYZ throughout that holds points and ",
                     "polygons do.")
  for (pg in list(page("resolution_profile.Rd"), page("determine_optimal_levels.Rd"), mf, fs)) {
    expect_match(pg, m_text, fixed = TRUE)
    expect_match(pg, red_text, fixed = TRUE)
  }
  expect_match(mf, paste0("and so is one whose polygons carry M (XYM or XYZM), which GEOS ",
                          "refuses."), fixed = TRUE)
  expect_match(mf, "So are blocks whose polygons carry M (XYM or XYZM), which GEOS refuses.",
               fixed = TRUE)
  expect_match(fs, paste0("So is one whose geometries carry M (XYM or XYZM), and so are the ",
                          "points it is reduced to when they mix dimensions."), fixed = TRUE)
})
