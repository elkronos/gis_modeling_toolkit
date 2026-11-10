# tests/testthat/test-fix-geometry-class-points.R
# ---------------------------------------------------------------------------
# A point layer can come in a geometry column of the general class
# sfc_GEOMETRY: sf::st_cast(x, "GEOMETRY") builds one, and so do c() of such a
# column with points and a point set with [[<- in place of the last other
# geometry of a mixed column.  sf::st_coordinates() refuses that
# class ("not implemented for objects of class sfc_GEOMETRY"), and
# sf::st_dimension() and sf::st_transform() stop with Rcpp's "Not compatible
# with STRSXP" on the column sf::st_cast() returns when it lacks the
# "classes" attribute sf::st_sfc() sets.  The fits, predict(), the metrics,
# the cross-validations, the folds, fold_separation(), the selection, the
# block-size sweep, area_of_applicability(), every move of such a layer to
# another CRS, and (last section) the cell summaries, the level choice and
# profile, the tessellations, the seeders and predict_surface() stopped with
# those messages.  Each now reads such a layer as the same points in
# sfc_POINT, and its answer is held to that layer's.
#
# The layer: 60 synthetic points in EPSG:32616 (seed 42).  It is given in
# two forms of the sfc_GEOMETRY column: as sf::st_cast() returns it, and with
# the "classes" attribute set, as sf::st_sfc() sets it; which of the two
# sf::st_cast() returns depends on the sf version.
# ---------------------------------------------------------------------------

.gcp_store <- new.env(parent = emptyenv())
.gcp_layers <- function() {
  if (is.null(.gcp_store$layers)) {
    d <- local({
      withr::local_seed(42, .rng_kind = "Mersenne-Twister",
                        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
      n  <- 60
      xy <- cbind(500000 + stats::runif(n, 0, 20000),
                  4000000 + stats::runif(n, 0, 20000))
      df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                       x3 = stats::rnorm(n))
      df$y <- 1 + 2 * df$x1 - df$x2 + stats::rnorm(n, sd = 0.3)
      cbind(df, X = xy[, 1], Y = xy[, 2])
    })
    pts <- sf::st_as_sf(d, coords = c("X", "Y"), crs = 32616)
    cast <- sf::st_cast(pts, "GEOMETRY")
    g <- sf::st_geometry(cast)
    attr(g, "classes") <- rep("POINT", length(g))
    classed <- cast
    suppressWarnings(sf::st_geometry(classed) <- g)
    .gcp_store$layers <- list(pts = pts, cast = cast, classed = classed)
  }
  .gcp_store$layers
}

# The two sfc_GEOMETRY forms of the layer.
.gcp_forms <- function() .gcp_layers()[c("cast", "classed")]

.gcp_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

.gcp_rf <- function(d, vars = c("x1", "x2"))
  fit_rf_model(d, "y", vars, num_trees = 20, num_threads = 1, seed = 1)

# What a GWR fit computes.  Its formula's environment is the frame of the
# call, which holds the layer as the caller gave it, and its info carries a
# cache environment; neither is part of the answer.
.gcp_gwr_core <- function(g) {
  info <- g$info
  info$.cache <- NULL
  list(data_sf = g$data_sf, n = g$n, formula = deparse(g$formula),
       response_var = g$response_var, predictor_vars = g$predictor_vars,
       info = info, sdf = as.data.frame(g$engine$SDF),
       fitted = g$engine$fitted, diagnostic = g$engine$GW.diagnostic,
       lm = stats::coef(g$engine$lm), centres = g$engine$centres)
}


test_that("the layers are the same points in an sfc_GEOMETRY column", {
  L <- .gcp_layers()
  for (f in names(.gcp_forms())) {
    d <- L[[f]]
    expect_s3_class(sf::st_geometry(d), "sfc_GEOMETRY")
    expect_true(all(sf::st_geometry_type(d, by_geometry = TRUE) == "POINT"))
  }
  expect_false(is.null(attr(sf::st_geometry(L$classed), "classes")))
})


test_that(".typed_points() casts points in sfc_GEOMETRY and leaves other layers", {
  L <- .gcp_layers()
  for (d in .gcp_forms()) {
    expect_identical(.typed_points(d), L$pts)
    expect_identical(.typed_points(sf::st_geometry(d)), sf::st_geometry(L$pts))
  }
  # Another name for the geometry column is kept.
  renamed <- sf::st_set_geometry(L$cast, "geom")
  expect_identical(.typed_points(renamed), sf::st_set_geometry(L$pts, "geom"))
  # An sfc_POINT layer is returned as it came.
  expect_identical(.typed_points(L$pts), L$pts)
  # Points with Z, and an empty point among the points, keep both.
  gz <- sf::st_sfc(sf::st_point(c(1, 2, 3)), sf::st_point(c(4, 5, 6)), crs = 32616)
  expect_identical(.typed_points(sf::st_cast(gz, "GEOMETRY")), gz)
  ge <- sf::st_sfc(sf::st_point(c(1, 2)), sf::st_point(), crs = 32616)
  ge_g <- ge
  class(ge_g) <- c("sfc_GEOMETRY", "sfc")
  expect_identical(.typed_points(ge_g), ge)
  # A layer that mixes points with a polygon is returned as it came.
  mixed <- L$pts
  sf::st_geometry(mixed) <- sf::st_sfc(
    c(list(sf::st_buffer(sf::st_geometry(L$pts)[[1]], 50)),
      unclass(sf::st_geometry(L$pts))[-1]), crs = 32616)
  expect_s3_class(sf::st_geometry(mixed), "sfc_GEOMETRY")
  expect_identical(.typed_points(mixed), mixed)
  # Anything that is not sf or sfc, and a layer of no rows, as they came.
  expect_identical(.typed_points(data.frame(a = 1)), data.frame(a = 1))
  z0 <- L$cast[0, ]
  expect_identical(.typed_points(z0), z0)
})


test_that("prep_model_data() reads the points as in sfc_POINT", {
  L <- .gcp_layers()
  ref <- prep_model_data(L$pts, "y", c("x1", "x2"))
  for (d in .gcp_forms())
    expect_identical(prep_model_data(d, "y", c("x1", "x2")), ref)
})


test_that("the fits give the sfc_POINT layer's answer", {
  L <- .gcp_layers()
  # A fit records its run time and its engine's environment, so two fits of
  # one layer are not identical(); what they compute is compared.
  rf_ref  <- .gcp_rf(L$pts)
  gwr_ref <- .gcp_quiet(fit_gwr_model(L$pts, "y", c("x1", "x2"),
                                      bandwidth = 20, adaptive = TRUE))
  sel_ref <- .gcp_quiet(gwr_model_selection(L$pts, "y", c("x1", "x2"),
                                            bandwidth = 20, adaptive = TRUE,
                                            quiet = TRUE))
  for (d in .gcp_forms()) {
    r <- .gcp_rf(d)
    expect_identical(fitted(r), fitted(rf_ref))
    expect_identical(predict(r, newdata = L$pts), predict(rf_ref, newdata = L$pts))
    expect_equal(r, rf_ref)
    g <- .gcp_quiet(fit_gwr_model(d, "y", c("x1", "x2"), bandwidth = 20,
                                  adaptive = TRUE))
    expect_identical(fitted(g), fitted(gwr_ref))
    expect_identical(coef(g), coef(gwr_ref))
    expect_identical(.gcp_gwr_core(g), .gcp_gwr_core(gwr_ref))
    expect_identical(.gcp_quiet(gwr_model_selection(d, "y", c("x1", "x2"),
                                                    bandwidth = 20, adaptive = TRUE,
                                                    quiet = TRUE)), sel_ref)
  }
})


test_that("fit_bayesian_spatial_model() hands brms the sfc_POINT layer's data", {
  skip_if_not_installed("brms")
  L <- .gcp_layers()
  # brms::brm() replaced by a recorder: the Stan data is what is compared.
  capture <- function(d) {
    cap <- new.env()
    local_mocked_bindings(
      brm = function(...) { cap$args <- list(...); structure(list(), class = "gcp_stub") },
      .package = "brms")
    .gcp_quiet(fit_bayesian_spatial_model(d, "y", "x1", gp_k = 5L,
                                          compute_loo = FALSE,
                                          check_convergence = FALSE))
    cap$args
  }
  ref <- capture(L$pts)
  for (d in .gcp_forms()) {
    got <- capture(d)
    expect_identical(got$data, ref$data)
    expect_identical(got$formula, ref$formula)
  }
})


test_that("the metrics and the residual Moran's I give the sfc_POINT layer's answer", {
  L <- .gcp_layers()
  rf_ref <- .gcp_rf(L$pts)
  mm_ref  <- .gcp_quiet(model_metrics(rf_ref, newdata = L$pts))
  ins_ref <- .gcp_quiet(evaluate_insample(list(a = rf_ref), newdata = L$pts))
  cmp_ref <- .gcp_quiet(compare_models(list(a = rf_ref), newdata = L$pts))
  mor_ref <- .gcp_quiet(residual_morans_i(rf_ref))
  mor_res <- .gcp_quiet(residual_morans_i(rf_ref, null = "residual"))
  for (d in .gcp_forms()) {
    expect_identical(.gcp_quiet(model_metrics(rf_ref, newdata = d)), mm_ref)
    r <- .gcp_rf(d)
    expect_identical(.gcp_quiet(evaluate_insample(list(a = r), newdata = d)), ins_ref)
    expect_identical(.gcp_quiet(compare_models(list(a = r), newdata = d)), cmp_ref)
    expect_identical(.gcp_quiet(residual_morans_i(r)), mor_ref)
    expect_identical(.gcp_quiet(residual_morans_i(r, null = "residual")), mor_res)
  }
})


test_that("make_folds() gives the sfc_POINT layer's folds", {
  L <- .gcp_layers()
  for (m in c("block_kfold", "random_kfold")) {
    ref <- make_folds(L$pts, k = 4, method = m, seed = 1)
    for (d in .gcp_forms())
      expect_identical(make_folds(d, k = 4, method = m, seed = 1), ref)
  }
})


test_that("the cross-validations give the sfc_POINT layer's answer", {
  L <- .gcp_layers()
  folds <- make_folds(L$pts, k = 4, method = "block_kfold", seed = 1)
  cvs <- function(d, ...)
    .gcp_quiet(cv_spatial(d, "y", c("x1", "x2"), fit_fn = function(tr) .gcp_rf(tr),
                          seed = 1, ...))
  cvr <- function(d) .gcp_quiet(cv_rf(d, "y", c("x1", "x2"), k = 4, seed = 1,
                                      num_trees = 20, num_threads = 1))
  cvg <- function(d) .gcp_quiet(cv_gwr(d, "y", c("x1", "x2"), k = 4, seed = 1,
                                       bandwidth = 20, adaptive = TRUE))
  cmp <- function(d) .gcp_quiet(compare_models_cv(
    d, "y", c("x1", "x2"), models = c("GWR", "RF"), k = 4, seed = 1,
    gwr_args = list(bandwidth = 20, adaptive = TRUE),
    rf_args = list(num_trees = 20, num_threads = 1), quiet = TRUE))
  ref <- list(cvs(L$pts, k = 4), cvs(L$pts, folds = folds), cvr(L$pts),
              cvg(L$pts), cmp(L$pts))
  for (d in .gcp_forms()) {
    expect_identical(cvs(d, k = 4), ref[[1]])
    expect_identical(cvs(d, folds = folds), ref[[2]])
    expect_identical(cvr(d), ref[[3]])
    expect_identical(cvg(d), ref[[4]])
    expect_identical(cmp(d), ref[[5]])
  }
})


test_that("select_features_forward() selects as on the sfc_POINT layer", {
  L <- .gcp_layers()
  sff <- function(d, ...) .gcp_quiet(select_features_forward(
    d, "y", c("x1", "x2", "x3"), fit_fn = function(tr, vars) .gcp_rf(tr, vars),
    k = 4, seed = 1, quiet = TRUE, ...))
  ref <- sff(L$pts)
  ref_se <- sff(L$pts, rule = "se")
  expect_true("x1" %in% ref$selected)
  for (d in .gcp_forms()) {
    expect_identical(sff(d), ref)
    expect_identical(sff(d, rule = "se"), ref_se)
  }
})


test_that("cv_block_size_sweep() gives the sfc_POINT layer's sweep", {
  L <- .gcp_layers()
  sw <- function(d) .gcp_quiet(cv_block_size_sweep(
    d, "y", c("x1", "x2"), fit_fn = function(tr) .gcp_rf(tr),
    block_sizes = c(4000, 8000), k = 4, seed = 1, quiet = TRUE))
  ref <- sw(L$pts)
  for (d in .gcp_forms()) expect_identical(sw(d), ref)
})


test_that("area_of_applicability() reads either layer as in sfc_POINT", {
  L <- .gcp_layers()
  aoa <- function(nd, tr) .gcp_quiet(area_of_applicability(
    nd, train_sf = tr, predictor_vars = c("x1", "x2"), coords = "yes"))
  ref <- aoa(L$pts, L$pts)
  for (d in .gcp_forms()) {
    expect_identical(aoa(L$pts, d), ref)
    # `aoa` returns the rows of `newdata` as they came, in their own class;
    # everything it computes is the sfc_POINT layer's.
    got <- aoa(d, L$pts)
    expect_identical(sf::st_geometry(got$aoa), sf::st_geometry(d))
    sf::st_geometry(got$aoa) <- sf::st_geometry(L$pts)
    expect_identical(got, ref)
  }
})


test_that("area_of_applicability() reads either layer for a forest fitted on its coordinates", {
  L <- .gcp_layers()
  # A forest fitted with include_coords = TRUE: area_of_applicability()
  # reads the coordinates as the forest's predictors ..x and ..y
  # (.aoa_add_coords()), in the training data and in `newdata`.
  rfc <- .gcp_quiet(fit_rf_model(L$pts, "y", c("x1", "x2"), num_trees = 20,
                                 num_threads = 1, seed = 1,
                                 include_coords = TRUE))
  for (coords in c("auto", "yes")) {
    ref    <- .gcp_quiet(area_of_applicability(L$pts, model = rfc, coords = coords))
    ref_tr <- .gcp_quiet(area_of_applicability(L$pts, model = rfc, train_sf = L$pts,
                                               coords = coords))
    for (d in .gcp_forms()) {
      expect_identical(.gcp_quiet(area_of_applicability(
        L$pts, model = rfc, train_sf = d, coords = coords)), ref_tr)
      # The rows of `newdata` come back in their own class.
      got <- .gcp_quiet(area_of_applicability(d, model = rfc, coords = coords))
      expect_identical(sf::st_geometry(got$aoa), sf::st_geometry(d))
      sf::st_geometry(got$aoa) <- sf::st_geometry(L$pts)
      expect_identical(got, ref)
    }
  }
})


test_that("`prediction_points` in an sfc_GEOMETRY column is read as in sfc_POINT", {
  L <- .gcp_layers()
  # The prediction points: a 7 x 7 grid over the layer's box, in the same
  # CRS.  sf::st_transform() stops on the column sf::st_cast() returns even
  # when the CRS does not change, so the reader casts it first.
  grid_pts <- sf::st_as_sf(
    expand.grid(X = seq(500000, 520000, length.out = 7),
                Y = seq(4000000, 4020000, length.out = 7)),
    coords = c("X", "Y"), crs = 32616)
  grid_cast <- sf::st_cast(grid_pts, "GEOMETRY")
  gg <- sf::st_geometry(grid_cast)
  attr(gg, "classes") <- rep("POINT", length(gg))
  grid_classed <- grid_cast
  suppressWarnings(sf::st_geometry(grid_classed) <- gg)
  forms <- list(grid_cast, grid_classed, sf::st_geometry(grid_cast))
  for (m in c("nndm", "knndm")) {
    ref <- .gcp_quiet(make_folds(L$pts, k = 4, method = m,
                                 prediction_points = grid_pts, seed = 1))
    for (pp in forms)
      expect_identical(.gcp_quiet(make_folds(L$pts, k = 4, method = m,
                                             prediction_points = pp, seed = 1)),
                       ref)
  }
  folds <- make_folds(L$pts, k = 4, method = "block_kfold", seed = 1)
  ref <- .gcp_quiet(fold_separation(folds, L$pts, prediction_points = grid_pts))
  for (pp in forms)
    expect_identical(.gcp_quiet(fold_separation(folds, L$pts,
                                                prediction_points = pp)), ref)
})


test_that("a layer that mixes points and a polygon is read as before", {
  L <- .gcp_layers()
  mixed <- L$pts
  sf::st_geometry(mixed) <- sf::st_sfc(
    c(list(sf::st_buffer(sf::st_geometry(L$pts)[[1]], 50)),
      unclass(sf::st_geometry(L$pts))[-1]), crs = 32616)
  expect_identical(.typed_points(mixed), mixed)
  # The polygon is reduced to a point inside it, and the points are kept.
  p <- prep_model_data(mixed, "y", c("x1", "x2"))
  ref <- prep_model_data(L$pts, "y", c("x1", "x2"))
  expect_s3_class(sf::st_geometry(p), "sfc_POINT")
  expect_identical(sf::st_coordinates(p)[-1, ], sf::st_coordinates(ref)[-1, ])
  expect_true(sf::st_intersects(sf::st_geometry(p)[1], sf::st_geometry(mixed)[1],
                                sparse = FALSE)[1, 1])
  expect_identical(make_folds(mixed, k = 4, method = "random_kfold", seed = 1)$folds,
                   make_folds(L$pts, k = 4, method = "random_kfold", seed = 1)$folds)
})


test_that("a layer in an sfc_GEOMETRY column is moved to another CRS as in sfc_POINT", {
  L <- .gcp_layers()
  # The same 60 points in degrees.  sf::st_transform() stops on the column
  # sf::st_cast() returns, so every move of a caller's layer
  # (.transform_keep_unlocated()) casts it first, and the layer comes back
  # in sfc_POINT.
  ll <- sf::st_transform(L$pts, 4326)
  ll_cast <- sf::st_cast(ll, "GEOMETRY")
  gl <- sf::st_geometry(ll_cast)
  attr(gl, "classes") <- rep("POINT", length(gl))
  ll_classed <- ll_cast
  suppressWarnings(sf::st_geometry(ll_classed) <- gl)
  forms <- list(ll_cast, ll_classed)
  rf_ll <- .gcp_rf(ll)
  rf_utm <- .gcp_rf(L$pts)
  folds <- make_folds(L$pts, k = 4, method = "block_kfold", seed = 1)
  ref <- list(
    ep    = .gcp_quiet(ensure_projected(ll)),
    ep_t  = .gcp_quiet(ensure_projected(ll, target_crs = 32616)),
    hc    = .gcp_quiet(harmonize_crs(ll, L$pts, prefer = "b")),
    mm    = .gcp_quiet(model_metrics(rf_ll, newdata = ll)),
    mm_u  = .gcp_quiet(model_metrics(rf_utm, newdata = ll)),
    pr    = .gcp_quiet(predict(rf_utm, newdata = ll)),
    fs    = .gcp_quiet(fold_separation(folds, ll)))
  for (d in forms) {
    expect_identical(.gcp_quiet(ensure_projected(d)), ref$ep)
    expect_identical(.gcp_quiet(ensure_projected(d, target_crs = 32616)), ref$ep_t)
    expect_identical(.gcp_quiet(harmonize_crs(d, L$pts, prefer = "b")), ref$hc)
    expect_identical(.gcp_quiet(model_metrics(rf_ll, newdata = d)), ref$mm)
    expect_identical(.gcp_quiet(model_metrics(rf_utm, newdata = d)), ref$mm_u)
    expect_identical(.gcp_quiet(predict(rf_utm, newdata = d)), ref$pr)
    expect_identical(.gcp_quiet(fold_separation(folds, d)), ref$fs)
  }
})


# ---- the cells, the levels, the tessellations, the seeds, the surface ------
# summarize_by_cell(), cell_values_at(), determine_optimal_levels(),
# resolution_profile(), create_voronoi_polygons(), build_tessellation()
# (voronoi, triangles), the k-means seeders and predict_surface() stopped on
# such a layer with sf's "not implemented for objects of class
# sfc_GEOMETRY", create_voronoi_polygons() on the column sf::st_cast()
# returns with "Not compatible with STRSXP" in the union of the hull, and
# assign_features_to_polygons() on that column with the same text in the
# join.  Each now reads the layer as the same points in sfc_POINT.

# Any point layer in the two sfc_GEOMETRY forms .gcp_layers() builds.
.gcp_forms_of <- function(x) {
  cast <- sf::st_cast(x, "GEOMETRY")
  g <- sf::st_geometry(cast)
  attr(g, "classes") <- rep("POINT", length(g))
  classed <- cast
  suppressWarnings(sf::st_geometry(classed) <- g)
  stopifnot(inherits(sf::st_geometry(cast), "sfc_GEOMETRY"),
            inherits(sf::st_geometry(classed), "sfc_GEOMETRY"))
  list(cast = cast, classed = classed)
}

# Nine square cells over the layer's box, the layer assigned to them, and a
# 6 x 6 grid of locations over the box.
.gcp_cells <- function() {
  if (is.null(.gcp_store$cells)) {
    L <- .gcp_layers()
    bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(L$pts)))
    cells <- .gcp_quiet(build_tessellation(L$pts, boundary = bnd, method = "square",
                                           approx_n_cells = 9, quiet = TRUE))$cells
    assigned <- .gcp_quiet(assign_features_to_polygons(L$pts, cells, quiet = TRUE))
    locs <- sf::st_as_sf(
      expand.grid(X = seq(500500, 519500, length.out = 6),
                  Y = seq(4000500, 4019500, length.out = 6)),
      coords = c("X", "Y"), crs = 32616)
    .gcp_store$cells <- list(boundary = bnd, cells = cells, assigned = assigned,
                             locs = locs)
  }
  .gcp_store$cells
}


test_that("assign_features_to_polygons() and summarize_by_cell() read the points as in sfc_POINT", {
  L <- .gcp_layers()
  C <- .gcp_cells()
  for (d in .gcp_forms())
    expect_identical(.gcp_quiet(assign_features_to_polygons(d, C$cells, quiet = TRUE)),
                     C$assigned)
  sbc <- function(d, ...) .gcp_quiet(summarize_by_cell(d, "y", "x1", ...))
  ref <- list(sbc(C$assigned, cells_sf = C$cells), sbc(C$assigned))
  for (d in .gcp_forms_of(C$assigned)) {
    expect_identical(sbc(d, cells_sf = C$cells), ref[[1]])
    expect_identical(sbc(d), ref[[2]])
  }
})


test_that("cell_values_at() reads locations, points and seeds as in sfc_POINT", {
  L <- .gcp_layers()
  C <- .gcp_cells()
  vals <- .gcp_quiet(summarize_by_cell(C$assigned, "y", cells_sf = C$cells))
  cva <- function(nd, pts, v = vals, ...) .gcp_quiet(cell_values_at(nd, v, pts, "y", ...))
  ref <- cva(C$locs, C$assigned)
  for (d in .gcp_forms_of(C$locs)) expect_identical(cva(d, C$assigned), ref)
  for (d in .gcp_forms_of(C$assigned)) expect_identical(cva(C$locs, d), ref)
  # The cells of seeds, placed by `seeds`: `values` then carries no polygons.
  seeds <- voronoi_seeds_kmeans(L$pts, k = 6, set_seed = 1)
  tv <- .gcp_quiet(build_tessellation(seeds, boundary = C$boundary, quiet = TRUE))
  seeds$poly_id <- tv$index
  a2 <- .gcp_quiet(assign_features_to_polygons(L$pts, tv$cells, quiet = TRUE))
  v2 <- sf::st_drop_geometry(.gcp_quiet(summarize_by_cell(a2, "y")))
  ref2 <- cva(C$locs, a2, v2, seeds = seeds)
  for (s in .gcp_forms_of(seeds))
    expect_identical(cva(C$locs, a2, v2, seeds = s), ref2)
})


test_that("the level choice and the profile read the points as in sfc_POINT", {
  L <- .gcp_layers()
  dol <- function(d, ...) .gcp_quiet(determine_optimal_levels(d, max_levels = 6L, ...))
  rp  <- function(d, ...) .gcp_quiet(resolution_profile(d, "y", n_levels = 4L, ...))
  ref <- list(dol(L$pts), dol(L$pts, response_var = "y", criterion = "combined"),
              rp(L$pts), rp(L$pts, partition = "kmeans"))
  ref_sel <- .gcp_quiet(select_resolution(ref[[3]]))
  for (d in .gcp_forms()) {
    expect_identical(dol(d), ref[[1]])
    expect_identical(dol(d, response_var = "y", criterion = "combined"), ref[[2]])
    p <- rp(d)
    expect_identical(p, ref[[3]])
    expect_identical(rp(d, partition = "kmeans"), ref[[4]])
    expect_identical(.gcp_quiet(select_resolution(p)), ref_sel)
  }
})


test_that("the Voronoi and triangle builders read the points as in sfc_POINT", {
  L <- .gcp_layers()
  C <- .gcp_cells()
  cvp <- function(d, ...) .gcp_quiet(create_voronoi_polygons(d, quiet = TRUE, ...))
  bt  <- function(d, ...) .gcp_quiet(build_tessellation(d, quiet = TRUE, ...))
  ref <- list(cvp(L$pts), cvp(L$pts, boundary = C$boundary),
              bt(L$pts), bt(L$pts, method = "triangles"),
              bt(L$pts, boundary = C$boundary, method = "triangles"))
  for (d in .gcp_forms()) {
    expect_identical(cvp(d), ref[[1]])
    expect_identical(cvp(d, boundary = C$boundary), ref[[2]])
    expect_identical(bt(d), ref[[3]])
    expect_identical(bt(d, method = "triangles"), ref[[4]])
    expect_identical(bt(d, boundary = C$boundary, method = "triangles"), ref[[5]])
  }
  # Points in that class as the `cover`.
  some <- L$pts[1:20, ]
  ref_c <- list(cvp(some, cover = L$pts), bt(some, cover = L$pts))
  for (d in .gcp_forms()) {
    expect_identical(cvp(some, cover = d), ref_c[[1]])
    expect_identical(bt(some, cover = d), ref_c[[2]])
  }
})


test_that("the k-means seeders read the points as in sfc_POINT", {
  L <- .gcp_layers()
  C <- .gcp_cells()
  gvs <- function(d, ...) .gcp_quiet(get_voronoi_seeds(method = "kmeans", n = 6,
                                                       sample_points = d, set_seed = 1, ...))
  vsk <- function(d, ...) .gcp_quiet(voronoi_seeds_kmeans(d, k = 6, set_seed = 1, ...))
  ref <- list(gvs(L$pts), gvs(L$pts, boundary = C$boundary), vsk(L$pts),
              vsk(L$pts, partition = "kmeans"))
  for (d in .gcp_forms()) {
    expect_identical(gvs(d), ref[[1]])
    expect_identical(gvs(d, boundary = C$boundary), ref[[2]])
    expect_identical(vsk(d), ref[[3]])
    expect_identical(vsk(d, partition = "kmeans"), ref[[4]])
  }
})


test_that("predict_surface() reads `grid` and `covariates` as in sfc_POINT", {
  L <- .gcp_layers()
  C <- .gcp_cells()
  rf  <- .gcp_rf(L$pts)
  gwr <- .gcp_quiet(fit_gwr_model(L$pts, "y", c("x1", "x2"), bandwidth = 20,
                                  adaptive = TRUE))
  ps <- function(fit, g, cv) .gcp_quiet(predict_surface(fit, grid = g, covariates = cv))
  for (fit in list(rf, gwr)) {
    ref <- ps(fit, C$locs, L$pts)
    for (d in .gcp_forms()) expect_identical(ps(fit, C$locs, d), ref)
    # The surface is the grid's rows as given, in their own class; what is
    # predicted on them is the sfc_POINT grid's.
    for (g in .gcp_forms_of(C$locs)) {
      got <- ps(fit, g, L$pts)
      expect_identical(sf::st_geometry(got), sf::st_geometry(g))
      sf::st_geometry(got) <- sf::st_geometry(C$locs)
      expect_identical(got, ref)
    }
  }
})


# ---- duplicates kept, the plots, the cached grid, multipoints ---------------
# build_tessellation(method = "triangles", keep_duplicates = TRUE) reads the
# points without the key of duplicates, which casts them, so it is held here
# on its own, and so is create_voronoi_polygons(keep_duplicates = TRUE).
# plot_tessellation_map(), plot_folds() and plot() of a fit returned a plot
# that stopped with "Not compatible with STRSXP" when it was built, on the
# column sf::st_cast() returns; create_grid_polygons_cached() stopped with
# that text on such points as `boundary`, where create_grid_polygons()
# refuses a boundary that is not polygonal.  A layer of MULTIPOINT features
# in an sfc_GEOMETRY column stopped the Voronoi builders, the join and what
# reads it, and a layer that mixes POINT and MULTIPOINT (sf holds it in an
# sfc_GEOMETRY column) stopped the triangle builder with sf's "not
# implemented for objects of class sfc_GEOMETRY".

# The two sfc_GEOMETRY forms of a layer of one geometry type.
.gcp_forms_typed <- function(x, type) {
  cast <- sf::st_cast(x, "GEOMETRY")
  g <- sf::st_geometry(cast)
  attr(g, "classes") <- rep(type, length(g))
  classed <- cast
  suppressWarnings(sf::st_geometry(classed) <- g)
  list(cast = cast, classed = classed)
}

# What a plot draws: the data of each layer once the plot is built.
.gcp_built <- function(p) ggplot2::ggplot_build(p)$data


test_that("the builders read the points as in sfc_POINT with duplicates kept", {
  L <- .gcp_layers()
  C <- .gcp_cells()
  bt  <- function(d, ...) .gcp_quiet(build_tessellation(d, quiet = TRUE, keep_duplicates = TRUE, ...))
  cvp <- function(d, ...) .gcp_quiet(create_voronoi_polygons(d, quiet = TRUE, keep_duplicates = TRUE, ...))
  ref <- list(bt(L$pts, method = "triangles"),
              bt(L$pts, boundary = C$boundary, method = "triangles"),
              bt(L$pts), cvp(L$pts))
  for (d in .gcp_forms()) {
    expect_identical(bt(d, method = "triangles"), ref[[1]])
    expect_identical(bt(d, boundary = C$boundary, method = "triangles"), ref[[2]])
    expect_identical(bt(d), ref[[3]])
    expect_identical(cvp(d), ref[[4]])
  }
})


test_that("the maps draw the points as in sfc_POINT", {
  skip_if_not_installed("ggplot2")
  L <- .gcp_layers()
  C <- .gcp_cells()
  ptm <- function(...) .gcp_quiet(.gcp_built(plot_tessellation_map(C$cells, ...)))
  folds <- .gcp_quiet(make_folds(L$pts, k = 4, seed = 1))
  pf <- function(...) .gcp_quiet(.gcp_built(plot_folds(folds, ...)))
  ref <- list(ptm(features_sf = L$pts), ptm(seeds_sf = L$pts),
              ptm(features_sf = L$pts, target_crs = 4326),
              ptm(boundary = L$pts), pf(L$pts), pf(L$pts, boundary = L$pts))
  for (d in .gcp_forms()) {
    expect_identical(ptm(features_sf = d), ref[[1]])
    expect_identical(ptm(seeds_sf = d), ref[[2]])
    expect_identical(ptm(features_sf = d, target_crs = 4326), ref[[3]])
    expect_identical(ptm(boundary = d), ref[[4]])
    expect_identical(pf(d), ref[[5]])
    expect_identical(pf(L$pts, boundary = d), ref[[6]])
  }
  # A fit whose layer is in that class, as a custom fit may carry it.
  rf  <- .gcp_rf(L$pts)
  gwr <- .gcp_quiet(fit_gwr_model(L$pts, "y", c("x1", "x2"), bandwidth = 20,
                                  adaptive = TRUE))
  pl <- function(fit, ...) .gcp_quiet(.gcp_built(plot(fit, ...)))
  ref_f <- list(pl(rf), pl(gwr, type = "coefficients"))
  for (d in .gcp_forms()) {
    rf_d <- rf
    rf_d$data_sf <- d
    gwr_d <- gwr
    gwr_d$data_sf <- d
    expect_identical(pl(rf_d), ref_f[[1]])
    expect_identical(pl(gwr_d, type = "coefficients"), ref_f[[2]])
  }
})


test_that("create_grid_polygons_cached() refuses such points as `boundary` as in sfc_POINT", {
  L <- .gcp_layers()
  cache <- new.env(parent = emptyenv())
  expect_error(create_grid_polygons_cached(L$pts, target_cells = 9, cache_env = cache),
               "'boundary' must be polygonal")
  for (d in .gcp_forms())
    expect_error(create_grid_polygons_cached(d, target_cells = 9, cache_env = cache),
                 "'boundary' must be polygonal")
})


test_that("MULTIPOINT features in an sfc_GEOMETRY column are read as in sfc_MULTIPOINT", {
  L <- .gcp_layers()
  C <- .gcp_cells()
  M <- sf::st_cast(L$pts, "MULTIPOINT")
  forms <- .gcp_forms_typed(M, "MULTIPOINT")
  for (d in forms) {
    expect_s3_class(sf::st_geometry(d), "sfc_GEOMETRY")
    expect_true(all(sf::st_geometry_type(d, by_geometry = TRUE) == "MULTIPOINT"))
  }
  cvp <- function(d, ...) .gcp_quiet(create_voronoi_polygons(d, quiet = TRUE, ...))
  bt  <- function(d, ...) .gcp_quiet(build_tessellation(d, quiet = TRUE, ...))
  afp <- function(d) .gcp_quiet(assign_features_to_polygons(d, C$cells, quiet = TRUE))
  sbc <- function(d) .gcp_quiet(summarize_by_cell(afp(d), "y", cells_sf = C$cells))
  ka  <- function(d) .gcp_quiet(kriging_adequacy(afp(d), "y", C$cells))
  rf  <- .gcp_rf(L$pts)
  ps  <- function(d) .gcp_quiet(predict_surface(rf, grid = C$locs, covariates = d))
  ref <- list(cvp(M), cvp(M, boundary = C$boundary), cvp(M, keep_duplicates = TRUE),
              bt(M), bt(M, method = "triangles"),
              bt(M, method = "triangles", keep_duplicates = TRUE),
              afp(M), sbc(M), ka(M), ps(M))
  for (d in forms) {
    expect_identical(cvp(d), ref[[1]])
    expect_identical(cvp(d, boundary = C$boundary), ref[[2]])
    expect_identical(cvp(d, keep_duplicates = TRUE), ref[[3]])
    expect_identical(bt(d), ref[[4]])
    expect_identical(bt(d, method = "triangles"), ref[[5]])
    expect_identical(bt(d, method = "triangles", keep_duplicates = TRUE), ref[[6]])
    expect_identical(afp(d), ref[[7]])
    expect_identical(sbc(d), ref[[8]])
    expect_identical(ka(d), ref[[9]])
    expect_identical(ps(d), ref[[10]])
  }
})


test_that("the triangles read a layer that mixes POINT and MULTIPOINT as the same features in sfc_MULTIPOINT", {
  L <- .gcp_layers()
  C <- .gcp_cells()
  g <- sf::st_geometry(L$pts)
  # Two points held as MULTIPOINT features, the rest as POINT: sf::st_sfc()
  # puts the mix in an sfc_GEOMETRY column.
  mixed <- sf::st_sfc(c(lapply(1:2, function(i) sf::st_multipoint(
                          sf::st_coordinates(g[i])[, 1:2, drop = FALSE])),
                        unclass(g)[-(1:2)]), crs = sf::st_crs(g))
  mix <- L$pts
  sf::st_geometry(mix) <- mixed
  expect_s3_class(sf::st_geometry(mix), "sfc_GEOMETRY")
  as_multi <- mix
  sf::st_geometry(as_multi) <- sf::st_cast(mixed, "MULTIPOINT")
  bt <- function(d, ...) .gcp_quiet(build_tessellation(d, method = "triangles",
                                                      quiet = TRUE, ...))
  expect_identical(bt(mix), bt(as_multi))
  expect_identical(bt(mix, keep_duplicates = TRUE), bt(as_multi, keep_duplicates = TRUE))
  expect_identical(bt(mix, boundary = C$boundary), bt(as_multi, boundary = C$boundary))
  # The triangles are those of the same points in sfc_POINT.
  expect_identical(bt(mix)$cells, bt(L$pts)$cells)
})
