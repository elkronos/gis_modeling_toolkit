# tests/testthat/test-fix-mixed-zm-models-plots.R
# ---------------------------------------------------------------------------
# sf reads a layer in the dimension of its first geometry: st_coordinates()
# of a point layer that mixes XY and XYZ points lays the values out of step
# (XYZ first) or stops on its "dimnames" (XY first), st_transform() of one
# with an XYZ geometry first stops with "OGR error", and GEOS stops on mixed
# polygons and on any geometry that carries M.  area_of_applicability(),
# predict_surface(), the predict() methods, plot_tessellation_map(),
# plot_folds() and residual_morans_i() now read such layers in x and y, and
# coerce_to_points() reduces a geometry that carries M in x and y with a
# warning.  Each mixed layer is held to the same rows in XY, to the bit,
# warnings included.  A warning of class spatialkit_zm_dropped is raised only
# where a copy of the caller's rows comes back in XY.
# ---------------------------------------------------------------------------

.mmp_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value,
       classes = vapply(got, function(w) class(w)[1L], character(1)),
       texts = vapply(got, conditionMessage, character(1)),
       rows = lapply(got, function(w) w$rows))
}

# `got` equals `ref` in value, warning classes and warning texts, and raises
# no warning about Z and M.
.mmp_same <- function(got, ref, info = NULL) {
  expect_identical(got$value, ref$value, info = info)
  expect_identical(got$classes, ref$classes, info = info)
  expect_identical(got$texts, ref$texts, info = info)
  expect_false("spatialkit_zm_dropped" %in% got$classes, info = info)
}

# `got` equals `ref`, save one warning of class spatialkit_zm_dropped that
# names `caller` and `arg` and carries `rows`.
.mmp_same_but_zm <- function(got, ref, caller, arg, rows, info = NULL) {
  expect_identical(got$value, ref$value, info = info)
  zm <- got$classes == "spatialkit_zm_dropped"
  expect_identical(sum(zm), 1L, info = info)
  expect_identical(got$classes[!zm], ref$classes, info = info)
  expect_identical(got$texts[!zm], ref$texts, info = info)
  expect_match(got$texts[zm], sprintf("%s(): `%s` mixes geometries", caller, arg),
               fixed = TRUE, info = info)
  expect_identical(got$rows[zm][[1L]], rows, info = info)
}

# Feature `i` of `l` given a third (or fourth) coordinate at its own x and y.
.mmp_row <- function(l, i, dim = "XYZ") {
  sfc <- inherits(l, "sfc")
  g  <- if (sfc) l else sf::st_geometry(l)
  gi <- sf::st_zm(g[[i]], drop = TRUE, what = "ZM")
  extra <- nchar(dim) - 2L
  add <- function(m) cbind(m, matrix(3, nrow(m), extra))
  co <- unclass(gi)
  gi <- switch(class(gi)[2L],
    POINT = sf::st_point(c(co, rep(3, extra)), dim = dim),
    MULTIPOINT = sf::st_multipoint(add(co), dim = dim),
    LINESTRING = sf::st_linestring(add(co), dim = dim),
    POLYGON = sf::st_polygon(lapply(co, add), dim = dim),
    MULTIPOLYGON = sf::st_multipolygon(lapply(co, function(p) lapply(p, add)), dim = dim))
  lst <- as.list(g); lst[[i]] <- gi
  out <- sf::st_sfc(lst, crs = sf::st_crs(g), recompute_bbox = TRUE)
  if (sfc) out else sf::st_set_geometry(l, out)
}
.mmp_allz <- function(l, dim = "XYZ") {
  for (i in seq_len(if (inherits(l, "sfc")) length(l) else nrow(l))) l <- .mmp_row(l, i, dim)
  l
}
.mmp_dims <- function(l) {
  g <- if (inherits(l, "sfc")) l else sf::st_geometry(l)
  unique(vapply(g, function(s) class(s)[1L], character(1)))
}

# Rows 1 and 5 (or 1 and 3 for the four quadrants), each given Z or M.
.mmp_cases <- function(rows = c(1L, 5L))
  expand.grid(row = rows, dim = c("XYZ", "XYM"), stringsAsFactors = FALSE)

.mmp_store <- new.env(parent = emptyenv())
.mmp_fixtures <- function() {
  if (is.null(.mmp_store$f)) {
    d <- local({
      withr::local_seed(7411, .rng_kind = "Mersenne-Twister",
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
    grid <- sf::st_sf(id = 1:16, geometry = sf::st_make_grid(bnd4, n = c(4, 4)))
    gpts <- sf::st_sf(id = 1:16, geometry = sf::st_centroid(sf::st_geometry(grid)))
    nd <- seq(1, nrow(pts), by = 4)
    .mmp_store$f <- list(pts = pts, polys = polys, bnd4 = bnd4, grid = grid,
                         gpts = gpts, nd = nd)
  }
  .mmp_store$f
}

# The values of an aoa object, and its rows without their geometry.
.mmp_aoa_vals <- function(r) {
  r <- unclass(r)
  r$aoa <- sf::st_drop_geometry(r$aoa)
  r
}

test_that("the premise: sf reads a mixed layer in the dimension of its first geometry", {
  F <- .mmp_fixtures()
  utils::capture.output(expect_error(sf::st_transform(.mmp_row(F$pts, 1), 4326)))
  expect_error(suppressWarnings(sf::st_coordinates(.mmp_row(F$pts, 5))))
  expect_s3_class(tryCatch(sf::st_point_on_surface(sf::st_geometry(.mmp_allz(F$polys[1:5, ], "XYM"))),
                           error = function(e) e), "error")
})

# ---- BUG-74 ---------------------------------------------------------------

test_that("area_of_applicability() reads mixed newdata as the same rows in XY", {
  skip_if_not_installed("ranger")
  F <- .mmp_fixtures()
  fitc <- fit_rf_model(F$pts, "z", c("a1", "a2"), include_coords = TRUE,
                       num_trees = 50, seed = 1)
  aoa <- function(l) .mmp_catch(area_of_applicability(l, fitc,
                                                      predictor_vars = c("a1", "a2"),
                                                      seed = 1))
  vals <- function(r) { r$value <- .mmp_aoa_vals(r$value); r }
  nd <- F$pts[F$nd, ]
  ref <- vals(aoa(nd))
  cs <- .mmp_cases()
  for (j in seq_len(nrow(cs))) {
    l <- .mmp_row(nd, cs$row[j], cs$dim[j])
    got <- aoa(l)
    .mmp_same(vals(got), ref, info = paste("training CRS", cs$row[j], cs$dim[j]))
    # Not moved: the rows come back with the geometry they had.
    expect_identical(sf::st_geometry(got$value$aoa), sf::st_geometry(l))
  }
  .mmp_same(vals(aoa(.mmp_allz(nd))), ref, info = "XYZ throughout")
  # Moved into the training CRS: the rows come back in XY, with the warning.
  nd_ll <- sf::st_transform(nd, 4326)
  ref <- aoa(nd_ll)
  for (j in seq_len(nrow(cs))) {
    got <- aoa(.mmp_row(nd_ll, cs$row[j], cs$dim[j]))
    .mmp_same_but_zm(got, ref, "area_of_applicability", "newdata", cs$row[j],
                     info = paste("degrees", cs$row[j], cs$dim[j]))
  }
  .mmp_same(aoa(.mmp_allz(nd_ll))[c("classes", "texts")],
            ref[c("classes", "texts")], info = "degrees XYZ throughout")
  expect_identical(.mmp_aoa_vals(aoa(.mmp_allz(nd_ll))$value), .mmp_aoa_vals(ref$value))
  # Polygons are reduced to points in XY, with no warning about Z or M.
  ndp <- F$polys[F$nd, ]
  ref <- aoa(ndp)
  for (j in seq_len(nrow(cs)))
    .mmp_same(aoa(.mmp_row(ndp, cs$row[j], cs$dim[j])), ref,
              info = paste("polygons", cs$row[j], cs$dim[j]))
})

test_that("area_of_applicability() reads mixed train_sf and newdata for the location index in XY", {
  F <- .mmp_fixtures()
  nd <- F$pts[F$nd, ]
  tr <- F$pts[-F$nd, ]
  folds <- rep(1:5, length.out = nrow(tr))
  aoa <- function(n, t) .mmp_catch(area_of_applicability(
    n, train_sf = t, predictor_vars = c("a1", "a2"), coords = "yes",
    folds = folds, seed = 1))
  vals <- function(r) { r$value <- .mmp_aoa_vals(r$value); r }
  cs <- .mmp_cases()
  for (kind in c("pts", "polys")) {
    trk <- F[[kind]][-F$nd, ]
    ndk <- F[[kind]][F$nd, ]
    ref_t <- aoa(nd, trk)
    ref_n <- vals(aoa(ndk, tr))
    for (j in seq_len(nrow(cs))) {
      .mmp_same(aoa(nd, .mmp_row(trk, cs$row[j], cs$dim[j])), ref_t,
                info = paste(kind, "train_sf", cs$row[j], cs$dim[j]))
      l <- .mmp_row(ndk, cs$row[j], cs$dim[j])
      got <- aoa(l, tr)
      .mmp_same(vals(got), ref_n, info = paste(kind, "newdata", cs$row[j], cs$dim[j]))
      # The location index returns newdata as passed.
      expect_identical(sf::st_geometry(got$value$aoa), sf::st_geometry(l))
    }
    .mmp_same(aoa(nd, .mmp_allz(trk)), ref_t, info = paste(kind, "train_sf XYZ throughout"))
  }
  # train_sf in the training data's own CRS-less space, newdata with a CRS:
  # the training layer is moved, in x and y.
  tr0 <- sf::st_set_crs(tr, NA)
  ref <- aoa(nd, tr0)
  for (i in c(1L, 5L))
    .mmp_same(aoa(nd, .mmp_row(tr0, i)), ref, info = paste("CRS-less train_sf", i))
})

test_that("area_of_applicability() with include_coords reads CRS-less mixed layers in XY", {
  skip_if_not_installed("ranger")
  F <- .mmp_fixtures()
  fitc <- fit_rf_model(F$pts, "z", c("a1", "a2"), include_coords = TRUE,
                       num_trees = 50, seed = 1)
  aoa <- function(l, ...) .mmp_catch(area_of_applicability(
    l, fitc, predictor_vars = c("a1", "a2"), seed = 1, ...))
  vals <- function(r) { r$value <- .mmp_aoa_vals(r$value); r }
  nd <- F$pts[F$nd, ]
  cs <- .mmp_cases()
  # newdata without a CRS that looks like lon/lat is reprojected into the
  # training CRS: its rows come back in XY, with the warning.
  nd_ll0 <- sf::st_set_crs(sf::st_transform(nd, 4326), NA)
  ref <- aoa(nd_ll0)
  for (j in seq_len(nrow(cs)))
    .mmp_same_but_zm(aoa(.mmp_row(nd_ll0, cs$row[j], cs$dim[j])), ref,
                     "area_of_applicability", "newdata", cs$row[j],
                     info = paste("CRS-less degrees", cs$row[j], cs$dim[j]))
  # newdata without a CRS that is stamped with the training CRS is not
  # moved: its rows come back as given, Z and M included.
  nd0 <- sf::st_set_crs(nd, NA)
  ref <- vals(aoa(nd0))
  for (j in seq_len(nrow(cs))) {
    l <- .mmp_row(nd0, cs$row[j], cs$dim[j])
    got <- aoa(l)
    .mmp_same(vals(got), ref, info = paste("CRS-less stamped", cs$row[j], cs$dim[j]))
    expect_identical(lapply(sf::st_geometry(got$value$aoa), unclass),
                     lapply(sf::st_geometry(l), unclass))
    expect_identical(.mmp_dims(got$value$aoa),
                     if (cs$row[j] == 1L) c(cs$dim[j], "XY") else c("XY", cs$dim[j]))
    expect_identical(sf::st_crs(got$value$aoa), sf::st_crs(F$pts))
  }
  # A training layer without a CRS that looks like lon/lat, given beside the
  # forest, is moved into newdata's CRS in x and y, with no word about Z or M.
  tr_ll0 <- sf::st_set_crs(sf::st_transform(F$pts, 4326), NA)
  ref <- vals(aoa(nd, train_sf = tr_ll0))
  for (j in seq_len(nrow(cs)))
    .mmp_same(vals(aoa(nd, train_sf = .mmp_row(tr_ll0, cs$row[j], cs$dim[j]))), ref,
              info = paste("CRS-less degrees train_sf", cs$row[j], cs$dim[j]))
})

# ---- BUG-81 ---------------------------------------------------------------

test_that("predict_surface() reads a mixed boundary and covariates as the same layers in XY", {
  F <- .mmp_fixtures()
  fit <- lm_spatial_fit(F$pts, "z", c("a1", "a2"))
  ps <- function(...) .mmp_catch(predict_surface(fit, n_cells = 100, ...))
  ref <- ps(boundary = F$bnd4, covariates = F$pts)
  cs <- .mmp_cases(c(1L, 3L))
  for (j in seq_len(nrow(cs)))
    .mmp_same(ps(boundary = .mmp_row(F$bnd4, cs$row[j], cs$dim[j]), covariates = F$pts),
              ref, info = paste("boundary", cs$row[j], cs$dim[j]))
  .mmp_same(ps(boundary = .mmp_allz(F$bnd4), covariates = F$pts), ref,
            info = "boundary XYZ throughout")
  .mmp_same(ps(boundary = sf::st_geometry(.mmp_row(F$bnd4, 1)), covariates = F$pts), ref,
            info = "boundary sfc")
  bnd_ll <- sf::st_transform(F$bnd4, 4326)
  ref_ll <- ps(boundary = bnd_ll, covariates = F$pts)
  for (j in seq_len(nrow(cs)))
    .mmp_same(ps(boundary = .mmp_row(bnd_ll, cs$row[j], cs$dim[j]), covariates = F$pts),
              ref_ll, info = paste("boundary degrees", cs$row[j], cs$dim[j]))

  ref <- ps(covariates = F$pts)
  cs <- .mmp_cases()
  for (j in seq_len(nrow(cs)))
    .mmp_same(ps(covariates = .mmp_row(F$pts, cs$row[j], cs$dim[j])), ref,
              info = paste("covariates", cs$row[j], cs$dim[j]))
  .mmp_same(ps(covariates = .mmp_allz(F$pts)), ref, info = "covariates XYZ throughout")
  cov_ll <- sf::st_transform(F$pts, 4326)
  ref_ll <- ps(covariates = cov_ll)
  for (j in seq_len(nrow(cs)))
    .mmp_same(ps(covariates = .mmp_row(cov_ll, cs$row[j], cs$dim[j])), ref_ll,
              info = paste("covariates degrees", cs$row[j], cs$dim[j]))
})

test_that("predict_surface() reads a mixed grid in XY and returns its rows as given or in XY", {
  F <- .mmp_fixtures()
  fit <- lm_spatial_fit(F$pts, "z", c("a1", "a2"))
  ps <- function(g, ...) .mmp_catch(predict_surface(fit, grid = g, covariates = F$pts, ...))
  cs <- .mmp_cases()
  # Polygons are reduced to points in XY, with no warning about Z or M.
  ref <- ps(F$grid)
  for (j in seq_len(nrow(cs)))
    .mmp_same(ps(.mmp_row(F$grid, cs$row[j], cs$dim[j])), ref,
              info = paste("polygon grid", cs$row[j], cs$dim[j]))
  .mmp_same(ps(.mmp_allz(F$grid)), ref, info = "polygon grid XYZ throughout")
  grid_ll <- sf::st_transform(F$grid, 4326)
  ref_ll <- ps(grid_ll)
  for (i in c(1L, 5L))
    .mmp_same(ps(.mmp_row(grid_ll, i)), ref_ll, info = paste("polygon grid degrees", i))
  # Points in the fit's CRS: the rows come back with the geometry they had.
  nogeom <- function(r) { r$value <- sf::st_drop_geometry(r$value); r }
  for (bnd in list(NULL, F$bnd4)) {
    ref <- nogeom(ps(F$gpts, boundary = bnd))
    for (j in seq_len(nrow(cs))) {
      l <- .mmp_row(F$gpts, cs$row[j], cs$dim[j])
      got <- ps(l, boundary = bnd)
      .mmp_same(nogeom(got), ref, info = paste("point grid", cs$row[j], cs$dim[j]))
      expect_identical(sf::st_geometry(got$value),
                       sf::st_geometry(l)[if (is.null(bnd)) TRUE else
                         lengths(sf::st_intersects(F$gpts, bnd)) > 0L])
    }
    l <- .mmp_allz(F$gpts)
    got <- ps(l, boundary = bnd)
    .mmp_same(nogeom(got), ref, info = "point grid XYZ throughout")
    expect_identical(.mmp_dims(got$value), "XYZ")
  }
  # Points moved into the fit's CRS: in XY, with the warning.
  gpts_ll <- sf::st_transform(F$gpts, 4326)
  ref <- ps(gpts_ll)
  for (j in seq_len(nrow(cs)))
    .mmp_same_but_zm(ps(.mmp_row(gpts_ll, cs$row[j], cs$dim[j])), ref,
                     "predict_surface", "grid", cs$row[j],
                     info = paste("point grid degrees", cs$row[j], cs$dim[j]))
  got <- ps(.mmp_allz(gpts_ll))
  expect_identical(sf::st_drop_geometry(got$value), sf::st_drop_geometry(ref$value))
  expect_false("spatialkit_zm_dropped" %in% got$classes)
})

test_that("predict_surface() reads a POINT grid that carries M, or has no CRS, in XY", {
  F <- .mmp_fixtures()
  fit <- lm_spatial_fit(F$pts, "z", c("a1", "a2"))
  ps <- function(g, ...) .mmp_catch(predict_surface(fit, grid = g, covariates = F$pts, ...))
  nogeom <- function(r) { r$value <- sf::st_drop_geometry(r$value); r }
  # A POINT grid XYM or XYZM throughout in the fit's CRS keeps its rows.
  for (bnd in list(NULL, F$bnd4)) {
    ref <- nogeom(ps(F$gpts, boundary = bnd))
    keep <- if (is.null(bnd)) TRUE else lengths(sf::st_intersects(F$gpts, bnd)) > 0L
    for (d in c("XYM", "XYZM")) {
      l <- .mmp_allz(F$gpts, d)
      got <- ps(l, boundary = bnd)
      .mmp_same(nogeom(got), ref, info = paste("point grid", d, is.null(bnd)))
      expect_identical(lapply(sf::st_geometry(got$value), unclass),
                       lapply(sf::st_geometry(l)[keep], unclass))
      expect_identical(.mmp_dims(got$value), d)
    }
  }
  cs <- .mmp_cases()
  # Without a CRS and looking like lon/lat, the grid is reprojected into the
  # fit's CRS: in XY, with the warning.
  g_ll0 <- sf::st_set_crs(sf::st_transform(F$gpts, 4326), NA)
  ref <- ps(g_ll0)
  for (j in seq_len(nrow(cs)))
    .mmp_same_but_zm(ps(.mmp_row(g_ll0, cs$row[j], cs$dim[j])), ref,
                     "predict_surface", "grid", cs$row[j],
                     info = paste("CRS-less degrees", cs$row[j], cs$dim[j]))
  # Without a CRS and stamped with the fit's, the grid is not moved: its
  # rows come back as given, Z and M included, as area_of_applicability()
  # returns stamped newdata.
  g0 <- sf::st_set_crs(F$gpts, NA)
  ref <- nogeom(ps(g0))
  for (j in seq_len(nrow(cs))) {
    l <- .mmp_row(g0, cs$row[j], cs$dim[j])
    got <- ps(l)
    .mmp_same(nogeom(got), ref, info = paste("CRS-less stamped", cs$row[j], cs$dim[j]))
    expect_identical(lapply(sf::st_geometry(got$value), unclass),
                     lapply(sf::st_geometry(l), unclass))
    expect_identical(sf::st_crs(got$value), sf::st_crs(F$pts))
  }
})

# ---- BUG-83 ---------------------------------------------------------------

test_that("predict() reads mixed newdata in another CRS as the same rows in XY", {
  skip_if_not_installed("ranger")
  F <- .mmp_fixtures()
  fitc <- fit_rf_model(F$pts, "z", c("a1", "a2"), include_coords = TRUE,
                       num_trees = 50, seed = 1)
  nd <- F$pts[F$nd, ]
  nd_ll <- sf::st_transform(nd, 4326)
  cs <- .mmp_cases()
  for (l0 in list(nd, nd_ll)) {
    ref <- .mmp_catch(predict(fitc, l0))
    for (j in seq_len(nrow(cs)))
      .mmp_same(.mmp_catch(predict(fitc, .mmp_row(l0, cs$row[j], cs$dim[j]))), ref,
                info = paste("rf", sf::st_crs(l0)$epsg, cs$row[j], cs$dim[j]))
    .mmp_same(.mmp_catch(predict(fitc, .mmp_allz(l0))), ref, info = "rf XYZ throughout")
  }
  ref <- .mmp_catch(evaluate_insample(list(rf = fitc), newdata = nd_ll))
  for (i in c(1L, 5L))
    .mmp_same(.mmp_catch(evaluate_insample(list(rf = fitc), newdata = .mmp_row(nd_ll, i))),
              ref, info = paste("evaluate_insample", i))
  ref <- .mmp_catch(compare_models(list(rf = fitc), newdata = nd_ll))
  .mmp_same(.mmp_catch(compare_models(list(rf = fitc), newdata = .mmp_row(nd_ll, 1))),
            ref, info = "compare_models")
})

test_that("predict() on a GWR reads mixed newdata in another CRS as the same rows in XY", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  F <- .mmp_fixtures()
  fitg <- suppressWarnings(fit_gwr_model(F$pts, "z", c("a1", "a2")))
  nd_ll <- sf::st_transform(F$pts[F$nd, ], 4326)
  ref <- .mmp_catch(predict(fitg, nd_ll))
  cs <- .mmp_cases()
  for (j in seq_len(nrow(cs)))
    .mmp_same(.mmp_catch(predict(fitg, .mmp_row(nd_ll, cs$row[j], cs$dim[j]))), ref,
              info = paste("gwr", cs$row[j], cs$dim[j]))
})

test_that("predict() on a bayesian_fit reads mixed newdata in another CRS in XY", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan smoke tests")
  skip_if_not_installed("brms")
  F <- .mmp_fixtures()
  fit <- NULL
  utils::capture.output(suppressWarnings(suppressMessages(
    fit <- fit_bayesian_spatial_model(F$pts[1:60, ], response_var = "z",
                                      predictor_vars = "a1", chains = 1,
                                      iter = 200, warmup = 100, cores = 1,
                                      compute_loo = FALSE,
                                      check_convergence = FALSE, seed = 1234))),
    type = "output")
  nd_ll <- sf::st_transform(F$pts[F$nd, ], 4326)
  ref <- .mmp_catch(predict(fit, nd_ll))
  for (i in c(1L, 5L))
    .mmp_same(.mmp_catch(predict(fit, .mmp_row(nd_ll, i))), ref, info = paste("bayes", i))
})

# ---- BUG-82 ---------------------------------------------------------------

test_that("plot_tessellation_map() draws mixed layers as the same layers in XY", {
  skip_if_not_installed("ggplot2")
  F <- .mmp_fixtures()
  cells <- sf::st_sf(cell_id = 1:16, geometry = sf::st_geometry(F$grid))
  seeds <- F$gpts
  built <- function(p) ggplot2::ggplot_build(p)$data
  pb <- function(...) .mmp_catch(built(plot_tessellation_map(...)))
  nogeom <- function(r) {
    r$value <- lapply(r$value, function(d) d[setdiff(names(d), "geometry")])
    r
  }
  cs <- .mmp_cases()
  for (tcrs in list(NULL, 4326)) {
    ref <- pb(cells, target_crs = tcrs, labels = TRUE)
    for (j in seq_len(nrow(cs)))
      .mmp_same(pb(.mmp_row(cells, cs$row[j], cs$dim[j]), target_crs = tcrs, labels = TRUE),
                ref, info = paste("cells", tcrs, cs$row[j], cs$dim[j]))
    .mmp_same(nogeom(pb(.mmp_allz(cells), target_crs = tcrs, labels = TRUE)), nogeom(ref),
              info = "cells XYZ throughout")
    # Cells that carry M are labelled at the points of the same cells in XY.
    for (d in c("XYM", "XYZM"))
      .mmp_same(nogeom(pb(.mmp_allz(cells, d), target_crs = tcrs, labels = TRUE)),
                nogeom(ref), info = paste("cells", d, "throughout"))
    ref <- pb(cells, seeds_sf = seeds, features_sf = F$pts, boundary = F$bnd4,
              target_crs = tcrs)
    for (j in seq_len(nrow(cs))) {
      r <- cs$row[j]; d <- cs$dim[j]
      .mmp_same(pb(cells, seeds_sf = .mmp_row(seeds, r, d), features_sf = F$pts,
                   boundary = F$bnd4, target_crs = tcrs), ref,
                info = paste("seeds", tcrs, r, d))
      .mmp_same(pb(cells, seeds_sf = seeds, features_sf = .mmp_row(F$pts, r, d),
                   boundary = F$bnd4, target_crs = tcrs), ref,
                info = paste("features", tcrs, r, d))
      .mmp_same(pb(cells, seeds_sf = seeds, features_sf = F$pts,
                   boundary = .mmp_row(F$bnd4, min(r, 3L), d), target_crs = tcrs), ref,
                info = paste("boundary", tcrs, r, d))
    }
  }
})

test_that("plot_folds() draws mixed points and boundary as the same layers in XY", {
  skip_if_not_installed("ggplot2")
  F <- .mmp_fixtures()
  folds <- make_folds(F$pts, k = 4, method = "block_kfold", seed = 1)
  pb <- function(...) .mmp_catch(ggplot2::ggplot_build(plot_folds(folds, ...))$data)
  nogeom <- function(r) {
    r$value <- lapply(r$value, function(d) d[setdiff(names(d), "geometry")])
    r
  }
  bnd_ll <- sf::st_transform(F$bnd4, 4326)
  ref <- pb(F$pts, boundary = bnd_ll)
  cs <- .mmp_cases()
  for (j in seq_len(nrow(cs))) {
    .mmp_same(pb(.mmp_row(F$pts, cs$row[j], cs$dim[j]), boundary = bnd_ll), ref,
              info = paste("points", cs$row[j], cs$dim[j]))
    .mmp_same(pb(F$pts, boundary = .mmp_row(bnd_ll, min(cs$row[j], 3L), cs$dim[j])), ref,
              info = paste("boundary", cs$row[j], cs$dim[j]))
  }
  .mmp_same(pb(F$pts, boundary = sf::st_geometry(.mmp_row(bnd_ll, 1))), ref,
            info = "boundary sfc")
  .mmp_same(nogeom(pb(.mmp_allz(F$pts), boundary = .mmp_allz(bnd_ll))), nogeom(ref),
            info = "XYZ throughout")
})

test_that("plot() of a custom fit draws its mixed layer as the same layer in XY", {
  skip_if_not_installed("ggplot2")
  F <- .mmp_fixtures()
  pb <- function(l, ...) .mmp_catch(ggplot2::ggplot_build(
    plot(lm_spatial_fit(l, "z", c("a1", "a2")), ...))$data)
  cs <- .mmp_cases()
  for (kind in c("pts", "polys")) {
    for (l0 in list(F[[kind]][1:100, ], sf::st_transform(F[[kind]][1:100, ], 4326))) {
      ref <- pb(l0)
      for (j in seq_len(nrow(cs)))
        .mmp_same(pb(.mmp_row(l0, cs$row[j], cs$dim[j])), ref,
                  info = paste(kind, sf::st_crs(l0)$epsg, cs$row[j], cs$dim[j]))
    }
  }
  if (requireNamespace("gstat", quietly = TRUE)) {
    ref <- pb(F$pts, type = "variogram")
    .mmp_same(pb(.mmp_row(F$pts, 1), type = "variogram"), ref, info = "variogram")
  }
})

# ---- BUG-173 --------------------------------------------------------------

test_that("residual_morans_i() reads a custom fit's mixed layer as the same layer in XY", {
  F <- .mmp_fixtures()
  rmi <- function(l) .mmp_catch(residual_morans_i(lm_spatial_fit(l, "z", c("a1", "a2"))))
  cs <- .mmp_cases()
  for (kind in c("pts", "polys")) {
    for (l0 in list(F[[kind]], sf::st_transform(F[[kind]], 4326))) {
      ref <- rmi(l0)
      expect_true(is.list(ref$value))
      for (j in seq_len(nrow(cs)))
        .mmp_same(rmi(.mmp_row(l0, cs$row[j], cs$dim[j])), ref,
                  info = paste(kind, sf::st_crs(l0)$epsg, cs$row[j], cs$dim[j]))
      .mmp_same(rmi(.mmp_allz(l0)), ref, info = paste(kind, "XYZ throughout"))
    }
  }
})

# ---- BUG-178 --------------------------------------------------------------

# `got` equals `ref` (the XY copy's result, which raises no warning about Z
# or M), save one warning of class spatialkit_zm_dropped about M.
.mmp_same_m <- function(got, ref, rows, n_m = length(rows), info = NULL) {
  expect_identical(got$value, ref$value, info = info)
  zm <- got$classes == "spatialkit_zm_dropped"
  expect_identical(sum(zm), 1L, info = info)
  expect_identical(got$classes[!zm], ref$classes, info = info)
  expect_identical(got$texts[!zm], ref$texts, info = info)
  expect_match(got$texts[zm],
               sprintf("coerce_to_points(): %d of the %d geometries of `x` carry M",
                       n_m, nrow(got$value)), fixed = TRUE, info = info)
  expect_identical(got$rows[zm][[1L]], rows, info = info)
}

test_that("coerce_to_points() reduces geometries that carry M in x and y, with a warning", {
  F <- .mmp_fixtures()
  polys <- F$polys[1:12, ]
  mpts <- sf::st_set_geometry(polys, sf::st_cast(sf::st_geometry(polys), "MULTIPOINT"))
  lines <- sf::st_set_geometry(polys, sf::st_cast(sf::st_geometry(polys), "LINESTRING"))
  pts <- F$pts[1:12, ]
  modes <- c("auto", "centroid", "point_on_surface", "line_midpoint", "bbox_center")
  layers <- list(polys = polys, mpts = mpts, lines = lines, pts = pts)
  for (crs in c(32617, 4326)) for (kind in names(layers)) {
    l0 <- sf::st_transform(layers[[kind]], crs)
    for (m in modes) {
      if (m == "line_midpoint" && kind == "mpts") next
      ref <- .mmp_catch(coerce_to_points(l0, m))
      # XYZ throughout takes its old path: no warning about Z or M.
      expect_false("spatialkit_zm_dropped" %in%
                     .mmp_catch(coerce_to_points(.mmp_allz(l0), m))$classes,
                   info = paste(crs, kind, m, "XYZ"))
      for (d in c("XYM", "XYZM")) {
        info <- paste(crs, kind, m, d)
        l <- .mmp_allz(l0, d)
        got <- .mmp_catch(coerce_to_points(l, m))
        if (m == "bbox_center") {
          .mmp_same(got, ref, info = info)
        } else if (m == "auto" && kind == "pts") {
          # POINT rows under "auto" are returned as they are, with their M
          # (the layer's z and m ranges are recomputed, as they were).
          expect_identical(sf::st_drop_geometry(got$value), sf::st_drop_geometry(l),
                           info = info)
          expect_identical(lapply(sf::st_geometry(got$value), unclass),
                           lapply(sf::st_geometry(l), unclass), info = info)
          expect_identical(.mmp_dims(got$value), d, info = info)
          expect_identical(got$classes, character(0), info = info)
        } else {
          .mmp_same_m(got, .mmp_catch(coerce_to_points(sf::st_zm(l), m)),
                      seq_len(nrow(l)), info = info)
          expect_identical(.mmp_dims(got$value), "XY", info = info)
        }
      }
    }
  }
  # An empty geometry in M among XY rows is not counted; in an XYM layer it
  # is left out of the rows.
  e <- polys
  g <- sf::st_geometry(e); g[[3]] <- sf::st_polygon(dim = "XYM")
  e <- sf::st_set_geometry(e, sf::st_sfc(as.list(g), crs = sf::st_crs(e)))
  got <- .mmp_catch(coerce_to_points(e, "auto"))
  expect_false("spatialkit_zm_dropped" %in% got$classes)
  expect_identical(got$value, coerce_to_points(sf::st_zm(e), "auto"))
  em <- .mmp_allz(polys, "XYM")
  g <- sf::st_geometry(em); g[[3]] <- sf::st_polygon(dim = "XYM")
  em <- sf::st_set_geometry(em, sf::st_sfc(as.list(g), crs = sf::st_crs(em)))
  got <- .mmp_catch(coerce_to_points(em, "auto"))
  .mmp_same_m(got, .mmp_catch(coerce_to_points(sf::st_zm(em), "auto")),
              setdiff(seq_len(nrow(em)), 3L), n_m = nrow(em) - 1L)
  # A layer that mixes dimensions keeps BUG-69's warning.
  got <- .mmp_catch(coerce_to_points(.mmp_row(polys, 2, "XYM"), "auto"))
  expect_identical(got$classes, "spatialkit_zm_dropped")
  expect_match(got$texts, "`x` mixes geometries", fixed = TRUE)
})

test_that("coerce_to_points() reads empty geometries that carry M in x and y, without a word", {
  F <- .mmp_fixtures()
  modes <- c("auto", "centroid", "point_on_surface")
  # Polygons that are all empty and carry M.
  for (d in c("XYM", "XYZM")) {
    ae <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(
      list(sf::st_polygon(dim = d), sf::st_polygon(dim = d)), crs = 32617))
    for (m in modes) {
      got <- .mmp_catch(coerce_to_points(ae, m))
      .mmp_same(got, .mmp_catch(coerce_to_points(sf::st_zm(ae), m)),
                info = paste("all empty", d, m))
    }
  }
  # POINT rows that carry M beside an empty polygon that carries M, under
  # "auto": the points keep their M; the empty row is read in x and y.
  for (d in c("XYM", "XYZM")) {
    g <- as.list(sf::st_geometry(.mmp_allz(F$pts[1:5, ], d)))
    g[[6]] <- sf::st_polygon(dim = d)
    x <- sf::st_sf(id = 1:6, geometry = sf::st_sfc(g, crs = 32617))
    got <- .mmp_catch(coerce_to_points(x, "auto"))
    expect_identical(got$classes, character(0), info = d)
    expect_identical(lapply(sf::st_geometry(got$value)[1:5], unclass),
                     lapply(sf::st_geometry(x)[1:5], unclass), info = d)
    expect_identical(sf::st_geometry(got$value)[[6]],
                     sf::st_geometry(coerce_to_points(sf::st_zm(x), "auto"))[[6]], info = d)
    expect_identical(sf::st_drop_geometry(got$value), sf::st_drop_geometry(x), info = d)
  }
})

test_that("coerce_to_points() names every row that carries M, points included, under auto", {
  F <- .mmp_fixtures()
  for (d in c("XYM", "XYZM")) {
    mix <- .mmp_allz(rbind(F$pts[1:3, ], F$polys[4:6, ]), d)
    got <- .mmp_catch(coerce_to_points(mix, "auto"))
    .mmp_same_m(got, .mmp_catch(coerce_to_points(sf::st_zm(mix), "auto")), 1:6, info = d)
    expect_identical(.mmp_dims(got$value), "XY", info = d)
  }
})

test_that("kriging_adequacy() reads multipoints that carry M in x and y, without a word", {
  skip_if_not_installed("gstat")
  F <- .mmp_fixtures()
  cells <- sf::st_sf(cell_id = 1:16, geometry = sf::st_geometry(F$grid))
  asg <- assign_features_to_polygons(F$pts[1:150, ], cells, polygon_id_col = "cell_id",
                                     quiet = TRUE)
  g <- sf::st_geometry(asg)
  mp <- sf::st_sfc(lapply(seq_along(g), function(i)
    sf::st_multipoint(rbind(unclass(g[[i]]), unclass(g[[i]]) + 1))), crs = 32617)
  asgm <- sf::st_set_geometry(asg, mp)
  sac <- suppressWarnings(estimate_sac_range(F$pts, "z", predictor_vars = "a2",
                                             cutoff = 0.3))
  expect_true(is.finite(as.numeric(sac)))
  ka <- function(l) .mmp_catch(kriging_adequacy(l, "z", cells, sac = sac,
                                                id_col = "cell_id", seed = 1))
  ref <- ka(asgm)
  for (d in c("XYM", "XYZM"))
    .mmp_same(ka(.mmp_allz(asgm, d)), ref, info = d)
})

test_that("the callers of coerce_to_points() read polygons that carry M as in XY", {
  F <- .mmp_fixtures()
  m_cases <- c("XYM", "XYZM")
  same_xy <- function(f, l, info) {
    ref <- .mmp_catch(f(l))
    for (d in m_cases) .mmp_same(.mmp_catch(f(.mmp_allz(l, d))), ref, info = paste(info, d))
  }
  # select_features_forward(select_on = "split")
  same_xy(function(l) select_features_forward(l, "z", c("a1", "a2"),
                                              fit_fn = function(tr, v) lm_spatial_fit(tr, "z", v),
                                              k = 3, quiet = TRUE, select_on = "split"),
          F$polys, "select_features_forward")
  # residual_morans_i() of a custom fit
  same_xy(function(l) residual_morans_i(lm_spatial_fit(l, "z", c("a1", "a2"))),
          F$polys, "residual_morans_i")
  # predict_surface(grid = polygons)
  fit <- lm_spatial_fit(F$pts, "z", c("a1", "a2"))
  same_xy(function(g) predict_surface(fit, grid = g, covariates = F$pts), F$grid,
          "predict_surface")
  # area_of_applicability(): the location index returns newdata as passed.
  nd <- F$polys[F$nd, ]
  tr <- F$polys[-F$nd, ]
  folds <- rep(1:5, length.out = nrow(tr))
  aoa <- function(n, t) .mmp_aoa_vals(area_of_applicability(
    n, train_sf = t, predictor_vars = c("a1", "a2"), coords = "yes",
    folds = folds, seed = 1))
  same_xy(function(l) aoa(l, tr), nd, "area_of_applicability newdata")
  same_xy(function(l) aoa(F$pts[F$nd, ], l), tr, "area_of_applicability train_sf")
  # cell_values_at() and summarize_by_cell()
  cells <- sf::st_sf(cell_id = 1:16, geometry = sf::st_geometry(F$grid))
  asg <- assign_features_to_polygons(F$pts[1:150, ], cells, polygon_id_col = "cell_id",
                                     quiet = TRUE)
  pooled <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = cells,
                              estimate = "pooled")
  asgp <- sf::st_set_geometry(asg, sf::st_buffer(sf::st_geometry(asg), 10, nQuadSegs = 2))
  same_xy(function(l) cell_values_at(l, pooled, asg, "z", id_col = "cell_id"),
          F$polys[151:190, ], "cell_values_at newdata")
  same_xy(function(l) cell_values_at(F$pts[151:190, ], pooled, l, "z", id_col = "cell_id"),
          asgp, "cell_values_at points")
  same_xy(function(l) summarize_by_cell(l, "z", id_col = "cell_id", cells_sf = cells,
                                        estimate = "pooled"),
          asgp, "summarize_by_cell pooled")
  if (requireNamespace("gstat", quietly = TRUE)) {
    sac <- suppressWarnings(estimate_sac_range(F$pts, "z", predictor_vars = "a2",
                                               cutoff = 0.3))
    expect_true(is.finite(as.numeric(sac)))
    expect_false(any(grepl("fallback", suppressWarnings(summarize_by_cell(
      asgp, "z", id_col = "cell_id", estimate = "mean", deff = "variogram",
      sac = sac))$deff_method)))
    same_xy(function(l) summarize_by_cell(l, "z", id_col = "cell_id", estimate = "mean",
                                          deff = "variogram", sac = sac),
            asgp, "summarize_by_cell deff")
  }
  if (requireNamespace("ranger", quietly = TRUE))
    same_xy(function(l) compare_models_cv(l, "z", c("a1", "a2"), models = "RF", k = 3,
                                          seed = 11, rf_args = list(num_trees = 20),
                                          quiet = TRUE)$comparison,
            F$polys, "compare_models_cv")
})

test_that("the help says these layers are read in x and y", {
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  mixed <- "mix XY with XYZ or XYM (as rbind() of an XY and an XYZ layer gives)"
  aoa <- page("area_of_applicability.Rd")
  expect_match(aoa, "Both layers are read in x and y for the location", fixed = TRUE)
  expect_match(aoa, "When location is measured it is read in x and y, as newdata is.",
               fixed = TRUE)
  ps <- page("predict_surface.Rd")
  expect_match(ps, "The grid is read in x and y", fixed = TRUE)
  expect_match(ps, paste0("A POINT grid whose coordinates are not moved (one in the fit's ",
                          "CRS, or one without a CRS that is stamped with it) is returned ",
                          "with its rows as given, Z and M included"), fixed = TRUE)
  expect_match(page("kriging_adequacy.Rd"),
               paste0("Multipoints whose geometries carry M (XYM or XYZM) are reduced to ",
                      "points in x and y, since GEOS does not read M."), fixed = TRUE)
  expect_match(ps, paste0("Read in x and y: a boundary whose geometries mix XY with XYZ ",
                          "or XYM, or carry M, is read without its Z and M."), fixed = TRUE)
  for (nm in c("predict.rf_fit.Rd", "predict.gwr_fit.Rd", "predict.bayesian_fit.Rd"))
    expect_match(page(nm), paste0("A layer whose geometries ", mixed, " is read in x and y."),
                 fixed = TRUE, info = nm)
  expect_match(page("plot_tessellation_map.Rd"),
               paste0("Every layer is drawn in x and y: a layer whose geometries ", mixed,
                      " is drawn without its Z and M. The labels are placed in x and y, ",
                      "also on cells whose geometries carry M (XYM or XYZM)."), fixed = TRUE)
  expect_match(page("plot_folds.Rd"), "and so are boundary and the blocks.", fixed = TRUE)
  expect_match(page("plot.spatial_fit.Rd"),
               "The maps are drawn in x and y: a data_sf whose geometries mix XY with XYZ or XYM",
               fixed = TRUE)
  expect_match(page("residual_morans_i.Rd"),
               "or carry M (XYM or XYZM), is read in x and y, with no warning.", fixed = TRUE)
  expect_match(page("coerce_to_points.Rd"), "GEOS does not read M coordinates.", fixed = TRUE)
  expect_match(page("compare_models_cv.Rd"),
               "or carry M (XYM or XYZM), is read in x and y.", fixed = TRUE)
  expect_match(page("select_features_forward.Rd"),
               "or carry M (XYM or XYZM), is read in x and y.", fixed = TRUE)
  expect_match(page("cell_values_at.Rd"),
               "whose geometries carry M (XYM or XYZM), which GEOS does not read.", fixed = TRUE)
  sbc <- page("summarize_by_cell.Rd")
  expect_match(sbc, paste0("A layer whose geometries carry M (XYM or XYZM), which GEOS does ",
                           "not read, is read in x and y too."), fixed = TRUE)
  expect_match(sbc, paste0("is read without its Z and M. So is a layer whose geometries carry ",
                           "M (XYM or XYZM), which GEOS does not read."), fixed = TRUE)
})
