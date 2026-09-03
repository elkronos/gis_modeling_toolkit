# ===========================================================================
# Regression tests for the fifth audit pass.
#
# One block per finding.  Each one FAILED on the tree these fixes landed on,
# and each asserts the number or the condition that was wrong -- not merely
# that the call returns something.
# ===========================================================================

.p5_pts <- function(n = 60, seed = 1, crs = 32632, extent = 1000) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, extent), y = runif(n, 0, extent),
               a = rnorm(n), b = rnorm(n)),
    coords = c("x", "y"), crs = crs
  )
}


# --------------------------------------------------------------------------
# prep_model_data(): names, response overlap, coordinates
# --------------------------------------------------------------------------

test_that("a column name R would parse as an expression is refused", {
  # `reformulate(termlabels = "B5-B4")` builds `~ B5 - B4`: a DIFFERENT model,
  # fitted silently, with $predictor_vars still reporting "B5-B4".  Backticking
  # is not a fix, because the sp coercion GWmodel needs runs the names through
  # make.names() anyway and the formula and the data would then disagree.
  d <- .p5_pts()
  d[["B5-B4"]] <- rnorm(nrow(d))
  d$z <- 2 * d$a + rnorm(nrow(d), 0, 0.1)

  expect_error(prep_model_data(d, "z", c("a", "B5-B4")),
               "not syntactically valid R names")
  expect_error(prep_model_data(d, "z", c("a", "B5-B4")), "B5-B4", fixed = TRUE)
  # A name with a space died in str2lang() AFTER every validation had passed.
  names(d)[names(d) == "b"] <- "soil pH"
  expect_error(prep_model_data(d, "z", c("a", "soil pH")),
               "not syntactically valid R names")
  # Ordinary names are untouched.
  expect_s3_class(prep_model_data(d, "z", "a"), "sf")
})


test_that("the response cannot also be a predictor", {
  d <- .p5_pts(); d$z <- 2 * d$a + rnorm(nrow(d))
  # Leakage in the forest (OOB R2 near 1, response top of the importance
  # table), a silently reduced model in GWR, duplicated rows plus a phantom
  # "<none>" entry in the GWR selection table.  No backend caught it.
  expect_error(prep_model_data(d, "z", c("z", "a")),
               "also listed in 'predictor_vars'")
})


test_that("prep_model_data drops rows whose coordinates are unusable", {
  # st_geometry_type() calls an EMPTY POINT a "POINT", so the POINT checks let
  # it through and nothing ever looked at the coordinates.  One infinite
  # coordinate at predict time makes brms's GP boundary infinite and every
  # basis function NaN.
  d <- .p5_pts(n = 40)
  d$z <- 2 * d$a + rnorm(40)
  g <- sf::st_geometry(d)
  g[[5]] <- sf::st_point()
  sf::st_geometry(d) <- g

  out <- suppressWarnings(prep_model_data(d, "z", "a"))
  expect_equal(nrow(out), 39L)
  expect_false(any(sf::st_is_empty(out)))

  d2 <- .p5_pts(n = 40); d2$z <- rnorm(40)
  g2 <- sf::st_geometry(d2)
  g2[[7]] <- sf::st_point(c(Inf, 0))
  sf::st_geometry(d2) <- g2
  expect_equal(nrow(suppressWarnings(prep_model_data(d2, "z", "a"))), 39L)
})


# --------------------------------------------------------------------------
# CRS selection and the lon/lat heuristic
# --------------------------------------------------------------------------

test_that("a projection is chosen by measurement, not by centroid latitude", {
  # The old rule sent any wide extent with a mid-latitude centroid to an Albers
  # conic whose standard parallels came from the BBOX.  For a trans-equatorial
  # extent those straddle the equator, the cone constant collapses, and the
  # result distorted distances by ~15.6% -- an order of magnitude worse than
  # the single UTM zone it was rejecting (1.65%).  At lat_1 = -lat_2 exactly,
  # PROJ refuses the string and ensure_projected() aborted with an
  # internal-looking "invalid crs".
  set.seed(1)
  d <- sf::st_as_sf(
    data.frame(lon = c(runif(40, -80, -60), -80, -60),
               lat = c(rep(35, 40), -40, 40)),
    coords = c("lon", "lat"), crs = 4326)

  out <- suppressWarnings(ensure_projected(d))
  expect_false(is.na(sf::st_crs(out)))            # used to abort outright

  err_used <- .crs_distance_error(d, sf::st_crs(out))
  err_utm  <- .crs_distance_error(d, sf::st_crs(32619))
  expect_true(is.finite(err_used))
  # Never worse than the zone it was weighing against.
  expect_lte(err_used, err_utm + 1e-9)
  expect_lt(err_used, 0.05)
})


test_that("antimeridian data are not split across the world", {
  # EPSG:3857 is the one answer worse than doing nothing for wrapped data.
  pts <- sf::st_as_sf(data.frame(lon = c(-179.5, 179.5), lat = c(-10, -9)),
                      coords = c("lon", "lat"), crs = 4326)
  out <- suppressWarnings(ensure_projected(pts))

  d_true <- as.numeric(sf::st_distance(pts)[1, 2])
  d_proj <- as.numeric(stats::dist(sf::st_coordinates(out)))
  expect_lt(abs(d_proj / d_true - 1), 0.05)
  expect_lt(d_proj, 5e6)                          # not half a planet
})


test_that(".transform_or_stamp reprojects lon/lat rather than stamping it", {
  ll <- sf::st_as_sf(data.frame(x = c(9.10, 9.15, 9.20), y = c(48.7, 48.75, 48.72)),
                     coords = c("x", "y"))
  expect_warning(out <- .transform_or_stamp(ll, sf::st_crs(32632)),
                 "reprojected")
  ref <- suppressWarnings(ensure_projected(ll, target_crs = sf::st_crs(32632)))
  expect_equal(sf::st_coordinates(out), sf::st_coordinates(ref), tolerance = 1e-9)

  # A planar layer that does NOT look like lon/lat is still stamped, and says so.
  pl <- sf::st_as_sf(data.frame(x = c(1e6, 1.0001e6), y = c(5e6, 5.0001e6)),
                     coords = c("x", "y"))
  expect_warning(out2 <- .transform_or_stamp(pl, sf::st_crs(32632)),
                 "WITHOUT reprojection")
  expect_equal(sf::st_coordinates(out2), sf::st_coordinates(pl))
})


test_that("harmonize_crs and ensure_projected place CRS-less input identically", {
  ll  <- sf::st_as_sf(data.frame(x = c(9.10, 9.15), y = c(48.7, 48.75)),
                      coords = c("x", "y"))
  ref <- sf::st_as_sf(data.frame(x = c(5e5, 5.1e5), y = c(5.39e6, 5.40e6)),
                      coords = c("x", "y"), crs = 32632)

  h <- suppressWarnings(harmonize_crs(ll, ref))
  e <- suppressWarnings(ensure_projected(ll, target_crs = sf::st_crs(32632)))
  expect_equal(sf::st_coordinates(h$a), sf::st_coordinates(e), tolerance = 1e-9)
})


test_that("predict() replays a NEGATIVE crs_assumed decision too", {
  # A CRS-less training layer the heuristic declined to call lon/lat is fitted
  # in its own planar space.  Nothing recorded that, so predict() re-ran the
  # heuristic on newdata ALONE -- and a SUBSET of those same training rows,
  # whose own bbox sits inside the lon/lat envelope, was judged differently
  # from the whole: taken for degrees, reprojected, and predicted ~1e6 m from
  # where it was fitted.
  set.seed(3)
  # A cluster that WOULD pass the heuristic on its own, plus far points that
  # push the layer's bounding box outside the lon/lat envelope.
  near <- data.frame(x = runif(60,   20,  170), y = runif(60, 10, 80))
  far  <- data.frame(x = runif(60, 3000, 5000), y = runif(60, 2000, 4000))
  n <- 120
  d <- sf::st_as_sf(cbind(rbind(near, far), a = rnorm(n)), coords = c("x", "y"))
  d$z <- 2 * d$a + rnorm(n, 0, 0.2)
  expect_false(.looks_like_lonlat(d)$lonlat)       # planar as a whole

  fit <- suppressWarnings(suppressMessages(
    fit_rf_model(d, "z", "a", include_coords = TRUE, num_trees = 60, seed = 1)))
  expect_identical(attr(fit$data_sf, "crs_assumed"), "none")

  # A SUBSET whose own bbox sits inside the lon/lat envelope: the trap.
  keep <- seq_len(60)
  sub  <- d[keep, ]
  expect_true(.looks_like_lonlat(sub)$lonlat)
  p_sub <- suppressWarnings(predict(fit, newdata = sub))
  p_all <- suppressWarnings(predict(fit, newdata = d))
  expect_equal(p_sub, p_all[keep], tolerance = 1e-10)
})


# --------------------------------------------------------------------------
# NNDM
# --------------------------------------------------------------------------

test_that("nndm excludes a held-out point's co-located duplicate", {
  # FNN::get.knn() returns the query point's OWN index in place of one of its
  # exact duplicates.  Dropping the self entry and padding with Inf left the
  # twin missing from the neighbour list altogether, so the sweep could never
  # exclude it: the fold holding out a repeated measurement trained on its
  # exact-location twin while params$realised_distances reported a large
  # buffer for a fold that had none.
  xy   <- rbind(c(0, 0), c(0, 0), c(100, 0), c(200, 0), c(300, 0), c(400, 0))
  pts  <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = 1:6),
                       coords = c("x", "y"), crs = 3857)
  grid <- sf::st_as_sf(data.frame(x = c(-500, 900), y = 0),
                       coords = c("x", "y"), crs = 3857)

  f <- suppressMessages(make_folds(pts, method = "nndm", prediction_points = grid))
  D <- as.matrix(sf::st_distance(pts)); units(D) <- NULL

  expect_false(2L %in% f$folds[[1]]$train)        # the twin is excluded
  expect_false(1L %in% f$folds[[2]]$train)
  # Every fold's REPORTED realised distance is the true one.
  realised <- vapply(seq_along(f$folds), function(i) {
    tr <- f$folds[[i]]$train
    if (!length(tr)) NA_real_ else min(D[i, tr])
  }, numeric(1))
  expect_equal(f$params$realised_distances, realised, tolerance = 1e-8)
})


test_that("nndm honours min_train exactly when n * min_train is fractional", {
  # The sweep permits ceiling(n - 1 - rmin) removals but the neighbour matrix
  # had only floor(n - 1 - rmin) + 1 columns, so every fold that reached the
  # floor kept one training point more than the paper's procedure -- for any
  # odd n at the default min_train = 0.5.
  set.seed(11); n <- 75
  pts <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000)),
    coords = c("x", "y"), crs = 3857)
  # Prediction points far from the training data: the target distances are
  # long, so the sweep removes as many neighbours as the rule allows and the
  # min_train floor is actually reached.
  grid <- sf::st_as_sf(data.frame(x = c(-8000, 9000), y = c(-8000, 9000)),
                       coords = c("x", "y"), crs = 3857)

  f <- suppressMessages(make_folds(pts, method = "nndm", prediction_points = grid))
  sizes <- vapply(f$folds, function(z) length(z$train), integer(1))
  # The rule is (n - 1 - removed) > n * min_train, so the smallest permitted
  # training set is 37 -- floor() left the matrix one column short and stopped
  # the sweep at 38.
  expect_equal(min(sizes), 37L)
  expect_gt(min(sizes), n * 0.5 - 1)
})


test_that("nndm pointizes a polygon prediction layer", {
  # Point-to-POLYGON distance is zero for every cell containing a training
  # point, which pulls the target ECDF to zero, excludes far fewer neighbours,
  # and degenerates the distance-matched CV towards plain LOO.
  set.seed(5); n <- 60
  pts <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000)),
    coords = c("x", "y"), crs = 3857)
  cells <- sf::st_make_grid(pts, cellsize = 250)
  ctrs  <- suppressWarnings(sf::st_centroid(cells))

  f_poly <- suppressMessages(make_folds(pts, method = "nndm",
                                        prediction_points = sf::st_sf(geometry = cells)))
  f_ctr  <- suppressMessages(make_folds(pts, method = "nndm",
                                        prediction_points = sf::st_sf(geometry = ctrs)))
  expect_equal(f_poly$params$n_removed_total, f_ctr$params$n_removed_total)
  expect_equal(lapply(f_poly$folds, `[[`, "train"),
               lapply(f_ctr$folds, `[[`, "train"))
})


# --------------------------------------------------------------------------
# Folds: probe, guards, k
# --------------------------------------------------------------------------

test_that("the fold probe tolerates a row it could not measure", {
  # make_folds() drops an empty geometry and says so, but the probe was taken
  # BEFORE that filter, so it carried NaN coordinates -- and every later cv_*()
  # on the same data died in .check_fold_probe() with R's internal
  # "missing value where TRUE/FALSE needed".
  d <- .p5_pts(n = 40); d$z <- 2 * d$a + rnorm(40)
  g <- sf::st_geometry(d); g[[5]] <- sf::st_point(); sf::st_geometry(d) <- g

  f <- suppressWarnings(suppressMessages(make_folds(d, k = 3, seed = 1)))
  res <- suppressWarnings(suppressMessages(
    cv_spatial(d, "z", "a",
               fit_fn = function(tr) lm_spatial_fit(tr, "z", "a"), folds = f)))
  expect_gt(nrow(res$predictions), 0L)
})


test_that("block_kfold refuses a grid it cannot build", {
  # No cap on nx * ny: block_size in the wrong unit asked for 1e8-1e11 cells
  # and exhausted memory instead of being refused.  create_grid_polygons()
  # guards exactly this and names the CRS units.
  d <- .p5_pts(n = 60, extent = 40000)
  expect_error(make_folds(d, k = 4, method = "block_kfold", block_size = 5),
               "above the")
  expect_error(make_folds(d, k = 4, method = "block_kfold", block_size = 5),
               "metre", fixed = TRUE)
  # A sane block size still works.
  expect_type(suppressMessages(
    make_folds(d, k = 4, method = "block_kfold", block_size = 8000)), "list")
})


test_that("k must be a whole number", {
  d <- .p5_pts(n = 20)
  expect_error(make_folds(d, k = 2.5), "whole number")
  expect_error(make_folds(d, k = NA_real_), "whole number")
  expect_error(make_folds(d, k = c(2, 3)), "whole number")
})


test_that("hand-built folds may not overlap, and unknown IDs are reported", {
  d <- .p5_pts(n = 40); d$z <- 2 * d$a + rnorm(40)
  ids <- seq_len(40)
  bad <- list(list(train = ids, test = ids[1:20]),
              list(train = ids, test = ids[21:40]))
  expect_error(
    suppressMessages(cv_spatial(d, "z", "a",
      fit_fn = function(tr) lm_spatial_fit(tr, "z", "a"), folds = bad)),
    "BOTH")

  good <- list(list(train = ids[21:40], test = ids[1:20]),
               list(train = ids[1:20],  test = ids[21:40]))
  lines <- capture_spatialkit_log(suppressMessages(
    cv_spatial(d, "z", "a", fit_fn = function(tr) lm_spatial_fit(tr, "z", "a"),
               folds = list(list(train = c(ids[21:40], 900:910), test = ids[1:20]),
                            list(train = ids[1:20], test = ids[21:40])))))
  expect_true(log_has(lines, "not in the data"))
  expect_s3_class(suppressMessages(
    cv_spatial(d, "z", "a",
               fit_fn = function(tr) lm_spatial_fit(tr, "z", "a"),
               folds = good))$overall,
    "data.frame")
})


# --------------------------------------------------------------------------
# estimate_sac_range
# --------------------------------------------------------------------------

test_that("the variogram is fitted with a nugget", {
  # A nugget-free model forces the curve through the origin, and gstat's
  # default N/h^2 weights buy that constraint by collapsing the range: with a
  # 50% nugget the fitted range came back at about 0.45 of the truth, so
  # make_folds(auto_range = TRUE) built blocks less than half the correlation
  # length it reported.
  skip_if_not_installed("gstat")
  set.seed(21); n <- 400; a <- 100                # effective range 3a = 300
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  S  <- 0.5 * exp(-D / a) + 0.5 * diag(n)         # half nugget
  z  <- as.numeric(t(chol(S)) %*% rnorm(n))
  d  <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
                     coords = c("x", "y"), crs = 32632)

  r <- suppressWarnings(estimate_sac_range(d, "z", seed = 1))
  expect_true(is.finite(as.numeric(r)))
  expect_gt(as.numeric(r), 150)                   # the nugget-free fit gave ~66
  m <- attr(r, "variogram_model")
  expect_true("Nug" %in% as.character(m$model))
})


test_that("Z is dropped before any distance is measured", {
  skip_if_not_installed("gstat")
  set.seed(31); n <- 200
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  z  <- as.numeric(t(chol(exp(-D / 40) + diag(1e-6, n))) %*% rnorm(n))
  flat <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
                       coords = c("x", "y"), crs = 32632)
  tall <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2],
                                  el = runif(n, 0, 2000), z = z),
                       coords = c("x", "y", "el"), crs = 32632)

  # gstat uses EVERY coordinate dimension, so elevation was folded into each
  # lag and the "range" came back as a length in 3-D.
  expect_equal(as.numeric(suppressWarnings(estimate_sac_range(tall, "z", seed = 1))),
               as.numeric(suppressWarnings(estimate_sac_range(flat, "z", seed = 1))),
               tolerance = 1e-8)
})


test_that("one unusable geometry does not void the whole estimate", {
  skip_if_not_installed("gstat")
  set.seed(41); n <- 200
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  z  <- as.numeric(t(chol(exp(-D / 40) + diag(1e-6, n))) %*% rnorm(n))
  d  <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
                     coords = c("x", "y"), crs = 32632)
  clean <- suppressWarnings(estimate_sac_range(d, "z", seed = 1))

  g <- sf::st_geometry(d); g[[3]] <- sf::st_point(); sf::st_geometry(d) <- g
  holed <- suppressWarnings(estimate_sac_range(d, "z", seed = 1))
  expect_true(is.finite(as.numeric(holed)))       # used to be NA for the layer
  expect_lt(abs(as.numeric(holed) / as.numeric(clean) - 1), 0.15)
})


test_that("the directional maximum is used only when anisotropy is established", {
  skip_if_not_installed("gstat")
  # Isotropic field: each direction sees about a quarter of the pairs, so the
  # max of four is biased upward.  All four must fit, the ratio must be
  # notable, and the widest must stand above the all-pairs estimate.
  set.seed(51); n <- 300
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  z  <- as.numeric(t(chol(exp(-D / 30) + diag(1e-6, n))) %*% rnorm(n))
  d  <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
                     coords = c("x", "y"), crs = 32632)
  r  <- suppressWarnings(estimate_sac_range(d, "z", seed = 1))
  dirs <- attr(r, "directional")
  expect_false(is.null(attr(r, "anisotropy_used")))
  if (!isTRUE(attr(r, "anisotropy_used")))
    expect_lte(as.numeric(r), max(dirs, na.rm = TRUE) + 1e-9)
  if (sum(is.finite(dirs)) < 4L)
    expect_false(isTRUE(attr(r, "anisotropy_used")))
})


# --------------------------------------------------------------------------
# summarize_by_cell
# --------------------------------------------------------------------------

.p5_cells <- function(nc = 6) {
  do.call(rbind, lapply(seq_len(nc), function(i)
    sf::st_sf(poly_id = i, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
      c((i - 1) * 100, 0), c(i * 100, 0), c(i * 100, 100),
      c((i - 1) * 100, 100), c((i - 1) * 100, 0)))), crs = 32632))))
}

.p5_assigned <- function(nc = 6, np = 20, seed = 2, rho_pred = TRUE) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(nc), function(i) {
    sf::st_as_sf(data.frame(
      x = runif(np, (i - 1) * 100 + 5, i * 100 - 5), y = runif(np, 5, 95),
      poly_id = i, resp = rnorm(np),
      pred = if (rho_pred) rnorm(1, 0, 3) + rnorm(np, 0, 0.2) else rnorm(np)),
      coords = c("x", "y"), crs = 32632)
  }))
}

test_that("a variable named `n` is not replaced by the row count", {
  # summarise() makes each new column visible to the expressions after it, so
  # naming the row count `n` put it in scope for the across() calls: a column
  # literally called `n` was summarised as the row count, sd and se NA.
  set.seed(9)
  pts <- .p5_assigned(nc = 3, np = 5)
  pts$n <- rnorm(nrow(pts), 100, 5)
  out <- summarize_by_cell(pts, response_var = "n")

  truth <- tapply(pts$n, pts$poly_id, mean)
  expect_equal(out[["resp_mean_n"]], as.numeric(truth[as.character(out$poly_id)]),
               tolerance = 1e-10)
  expect_false(any(out[["resp_mean_n"]] == out$n))
  expect_false(anyNA(out[["..sd_resp_n"]]))
})


test_that("a predictor's ICC never weights the response", {
  # The fallback to the predictor's ICC is for the case with NO response.
  # Using it whenever the response's own ICC came out <= 0 meant that adding a
  # predictor changed the response's regression weight by a factor of 20.
  pts <- .p5_assigned()
  with_pred <- suppressWarnings(summarize_by_cell(
    pts, response_var = "resp", predictor_vars = "pred", deff = "kish"))
  alone <- suppressWarnings(summarize_by_cell(
    pts, response_var = "resp", deff = "kish"))
  expect_equal(with_pred$cell_weight, alone$cell_weight, tolerance = 1e-10)
})


test_that("a rejected sac cannot size a design effect", {
  skip_if_not_installed("gstat")
  # estimate_sac_range() attaches `variogram_model` even when it REFUSES the
  # range, so the fit can be inspected.  Reading it without checking the value
  # took a refused model as gospel: deff came out equal to n, cell_weight
  # collapsed to 1 and the SEs were inflated 40-50x.
  pts <- .p5_assigned(nc = 4, np = 25)
  fake <- structure(
    NA_real_, class = c("sac_range", "numeric"),
    rejected_reason = "fitted range exceeds the largest lag fitted",
    rejected_range = 3e4,
    variogram_model = gstat::vgm(psill = 1, model = "Exp", range = 1e4,
                                 nugget = 0),
    crs = sf::st_crs(32632))

  # Two: the supplied `sac` is refused, and the internal re-estimate on this
  # small fixture is refused too.
  expect_warning(
    expect_warning(out <- summarize_by_cell(pts, response_var = "resp",
                                            deff = "variogram", sac = fake),
                   "no usable range"),
    "no usable range")
  expect_null(attr(out, "deff_applied"))
  iid <- summarize_by_cell(pts, response_var = "resp", deff = 1)
  expect_equal(out[["..se_resp_resp"]], iid[["..se_resp_resp"]], tolerance = 1e-10)
})


test_that("the ICC guard matches its documentation", {
  # Docs: "at least 2 cells with 2+ observations; falls back to deff = 1
  # otherwise."  All-singleton cells give MSW = 0 and an ICC of exactly 1 -- a
  # design effect of n from data with no within-cell information at all.
  set.seed(13)
  singles <- sf::st_as_sf(
    data.frame(x = seq(5, 795, by = 100), y = 50,
               poly_id = seq_len(8), resp = rnorm(8)),
    coords = c("x", "y"), crs = 32632)
  out <- suppressWarnings(summarize_by_cell(singles, response_var = "resp",
                                            deff = "kish"))
  expect_null(attr(out, "deff_applied"))
})


test_that("duplicate cell IDs are reported", {
  pts <- .p5_assigned(nc = 3, np = 10)
  cells <- .p5_cells(nc = 3)
  cells$poly_id[3] <- 2L
  expect_warning(summarize_by_cell(pts, response_var = "resp", cells_sf = cells),
                 "duplicated value")
})


test_that("every per-cell vector follows the cells_sf join", {
  skip_if_not_installed("gstat")
  # A spatially correlated field, so the variogram actually fits and a
  # per-cell design effect exists to realign.
  set.seed(101)
  nc <- 4; np <- 40
  df <- do.call(rbind, lapply(seq_len(nc), function(i) data.frame(
    x = runif(np, (i - 1) * 100 + 5, i * 100 - 5), y = runif(np, 5, 95),
    poly_id = i)))
  D <- as.matrix(stats::dist(cbind(df$x, df$y)))
  df$resp <- as.numeric(t(chol(exp(-D / 30) + diag(1e-6, nrow(df)))) %*%
                          rnorm(nrow(df)))
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = 32632)

  plain <- suppressWarnings(summarize_by_cell(pts, response_var = "resp",
                                              deff = "variogram"))
  da_p  <- attr(plain, "deff_applied")
  skip_if(is.null(da_p) || is.null(da_p$rbar),
          "the variogram did not fit on this draw")

  cells <- .p5_cells(nc = nc)[c(4, 3, 2, 1), ]     # reversed row order
  out <- suppressWarnings(summarize_by_cell(pts, response_var = "resp",
                                            deff = "variogram", cells_sf = cells))
  da <- attr(out, "deff_applied")
  expect_equal(length(da$rbar), nrow(out))         # used to keep pre-join length

  # deff = 1 + (n - 1) * rbar has to hold ROW-WISE after the join; $rbar used
  # to be left in pre-join order, so deff[i] and rbar[i] described different
  # cells.
  ok <- is.finite(da$deff) & is.finite(da$rbar) & is.finite(out$n)
  expect_true(any(ok))
  expect_equal(da$deff[ok],
               pmin(pmax(1, 1 + (out$n[ok] - 1) * da$rbar[ok]), out$n[ok]),
               tolerance = 1e-6)
  # And the join really did reorder, so this is a discriminating check.
  expect_false(identical(out$poly_id, plain$poly_id))
})


# --------------------------------------------------------------------------
# Tessellation, IDs, seeding
# --------------------------------------------------------------------------

test_that("a point outside every cell gets NA, not the nearest cell", {
  set.seed(17)
  inner <- sf::st_as_sf(data.frame(x = runif(10, 0, 100), y = runif(10, 0, 100)),
                        coords = c("x", "y"), crs = 32632)
  outer <- sf::st_as_sf(data.frame(x = runif(30, 500, 900), y = runif(30, 500, 900)),
                        coords = c("x", "y"), crs = 32632)
  all_pts <- rbind(inner, outer)
  bnd <- sf::st_as_sfc(sf::st_bbox(inner))

  tess <- suppressMessages(build_tessellation(all_pts, boundary = bnd,
                                              method = "square",
                                              approx_n_cells = 4, clip = TRUE))
  expect_equal(sum(is.na(tess$index)), 30L)
  assigned <- suppressWarnings(assign_features_to_polygons(
    all_pts, tess$cells, polygon_id_col = "cell_id", keep_unassigned = TRUE))
  expect_equal(sum(is.na(assigned$cell_id)), sum(is.na(tess$index)))
})


test_that("a square target_cells grid has square cells", {
  bnd <- sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 200, ymax = 100),
                                   crs = sf::st_crs(32632)))
  g <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  bb <- sf::st_bbox(g$geometry[[1]])
  expect_equal(as.numeric(bb["xmax"] - bb["xmin"]),
               as.numeric(bb["ymax"] - bb["ymin"]), tolerance = 1e-8)
})


test_that("a hex grid's max_cells guard measures the grid it will build", {
  bnd <- sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000),
                                   crs = sf::st_crs(32632)))
  # st_make_grid() uses cellsize[1] alone for hexagons, but the guard used
  # BOTH -- so a grid refused via c(10, 10) was built via c(10, 100).
  expect_error(create_grid_polygons(bnd, cellsize = c(10, 10), type = "hex",
                                    max_cells = 200), "max_cells")
  expect_error(suppressMessages(
    create_grid_polygons(bnd, cellsize = c(10, 100), type = "hex",
                         max_cells = 200)), "max_cells")
})


test_that("stable polygon IDs survive a round trip through another CRS", {
  bnd <- sf::st_as_sfc(sf::st_bbox(c(xmin = 9, ymin = 48, xmax = 10, ymax = 49),
                                   crs = sf::st_crs(4326)))
  g0 <- suppressWarnings(create_grid_polygons(sf::st_sf(geometry = bnd),
                                              target_cells = 16, type = "square"))
  base <- ensure_stable_poly_id(g0)
  for (epsg in c(32632, 3857, 2154)) {
    round <- suppressWarnings(sf::st_transform(sf::st_transform(g0, epsg), sf::st_crs(g0)))
    got   <- ensure_stable_poly_id(round)
    # Same cell, same ID, whichever projection the layer arrived in.
    ctr_b <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(base))))
    ctr_g <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(got))))
    ord   <- order(round(ctr_g[, 1], 6), round(ctr_g[, 2], 6))
    ord_b <- order(round(ctr_b[, 1], 6), round(ctr_b[, 2], 6))
    expect_equal(got$poly_id[ord], base$poly_id[ord_b], label = paste("epsg", epsg))
  }
})


test_that("get_voronoi_seeds clusters in 2-D and drops unusable points", {
  set.seed(23)
  cloud <- sf::st_as_sf(
    data.frame(x = c(runif(30, 0, 50), runif(30, 50, 100)),
               y = c(runif(30, 0, 50), runif(30, 50, 100)),
               el = runif(60, 0, 5000)),
    coords = c("x", "y", "el"), crs = 32632)
  bnd <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(sf::st_zm(cloud, drop = TRUE))))
  s3 <- get_voronoi_seeds(bnd, method = "kmeans", n = 2, sample_points = cloud)
  s2 <- voronoi_seeds_kmeans(sf::st_zm(cloud, drop = TRUE), k = 2)
  expect_equal(ncol(sf::st_coordinates(s3)), 2L)
  expect_equal(sort(sf::st_coordinates(s3)[, 1]),
               sort(sf::st_coordinates(s2)[, 1]), tolerance = 1)

  g <- sf::st_geometry(cloud)
  g[[4]] <- sf::st_point(c(NA_real_, NA_real_, NA_real_))
  sf::st_geometry(cloud) <- g
  expect_warning(get_voronoi_seeds(bnd, method = "kmeans", n = 2,
                                   sample_points = cloud),
                 "non-finite coordinates")
})


test_that("assign_features_to_polygons says so when nothing is assigned", {
  pts   <- .p5_pts(n = 20)
  cells <- .p5_cells(nc = 2)
  far   <- sf::st_as_sf(data.frame(x = c(9e5, 9.1e5), y = c(5.9e6, 5.91e6)),
                        coords = c("x", "y"), crs = 32632)
  expect_warning(assign_features_to_polygons(far, cells),
                 "none of the")
})


# --------------------------------------------------------------------------
# GWR and RF
# --------------------------------------------------------------------------

test_that("the local collinearity check sees the intercept and singular windows", {
  skip_if_not_installed("GWmodel")
  # local_xmat omitted the intercept column, so an indicator constant INSIDE a
  # window was not detected as collinear with it; and is.finite(cn) && cn > 1e6
  # discarded kappa = Inf, the worst case there is.  The check also ran before
  # the bandwidth was chosen, using a stand-in window, so it never fired on the
  # default path.
  set.seed(29)
  cl <- rep(1:4, each = 50)
  d <- sf::st_as_sf(data.frame(
    x = c(runif(50, 0, 100), runif(50, 400, 500),
          runif(50, 0, 100), runif(50, 400, 500)),
    y = c(runif(50, 0, 100), runif(50, 0, 100),
          runif(50, 400, 500), runif(50, 400, 500)),
    a = rnorm(200), urban = as.numeric(cl %in% c(2, 4))),
    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + 3 * d$urban + rnorm(200, 0, 0.2)

  # The bandwidth is PINNED rather than selected.  What the check has to see is
  # a window in which `urban` is constant, and whether the backend's own
  # bandwidth search lands on such a window is a property of the backend, not
  # of this fix: 20 nearest neighbours sits well inside one cluster on this
  # layout, so the local design is singular for every sampled location under
  # any GWmodel version.
  # Collect every warning rather than nesting expect_warning(): a fit whose
  # windows are singular also trips the post-fit non-finite count, and whether
  # a given GWmodel version produces NaN coefficients (rather than erroring
  # inside inv()) is a backend property.  Only the spot-check warning is ours
  # to assert here.
  .warns <- function(expr) {
    w <- character(0)
    withCallingHandlers(
      suppressMessages(try(expr, silent = TRUE)),
      warning = function(cond) {
        w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
      })
    w
  }
  expect_true(any(grepl("singular or near-singular",
    .warns(fit_gwr_model(d, "z", c("a", "urban"),
                         adaptive = TRUE, bandwidth = 20)))))

  # And quiet where there is nothing to report: a bandwidth wide enough to span
  # the clusters gives every window both values of `urban`.
  expect_false(any(grepl("singular or near-singular",
    .warns(fit_gwr_model(d, "z", c("a", "urban"),
                         adaptive = TRUE, bandwidth = 199)))))

  # The intercept is what makes the constant indicator collinear, so a check on
  # the predictors alone cannot see it -- assert the arithmetic directly rather
  # than through a backend.
  # Take the window from a cluster where urban == 1: the column is CONSTANT but
  # not zero, so the predictors alone are well conditioned and only the design
  # WITH the intercept is singular.  (A cluster of urban == 0 would be an
  # all-zero column, which the old check happened to catch.)
  xy   <- sf::st_coordinates(d)
  xmat <- as.matrix(sf::st_drop_geometry(d)[, c("a", "urban")])
  i1   <- which(d$urban == 1)[1L]
  nn   <- i1[1L]
  nn   <- order((xy[, 1] - xy[i1, 1])^2 + (xy[, 2] - xy[i1, 2])^2)[1:20]
  expect_true(all(xmat[nn, "urban"] == 1))
  expect_lt(kappa(xmat[nn, ], exact = FALSE), 1e6)          # what it used to check
  expect_gt(kappa(cbind(1, xmat[nn, ]), exact = FALSE), 1e6)
  expect_lt(qr(cbind(1, xmat[nn, ]))$rank, 3L)              # genuinely rank-deficient
})


test_that("fitted.gwr_fit is not fooled by a predictor named yhat", {
  skip_if_not_installed("GWmodel")
  # has_coef_suffix inferred which GWmodel call produced the SDF from the
  # column names, so a PREDICTOR ending in "_coef" made a gwr.basic() SDF look
  # like a gwr.predict() one and switched off the model-term exclusion --
  # after which fitted() returned a coefficient surface.
  set.seed(37); n <- 60
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               runoff_coef = rnorm(n), yhat = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$runoff_coef + 0.5 * d$yhat + rnorm(n, 0, 0.2)
  ref <- d; names(ref)[names(ref) == "yhat"] <- "pred2"
  ref$z <- d$z

  f1 <- suppressMessages(fit_gwr_model(d,   "z", c("runoff_coef", "yhat")))
  f2 <- suppressMessages(fit_gwr_model(ref, "z", c("runoff_coef", "pred2")))
  expect_equal(as.numeric(fitted(f1)), as.numeric(fitted(f2)), tolerance = 1e-8)
})


test_that("an implausible fixed bandwidth is called out", {
  skip_if_not_installed("GWmodel")
  # prep_model_data() projects geographic input before the bandwidth is used,
  # so 0.2 "degrees" became 0.2 metres: every local window empty, all
  # coefficients NaN, summary() reporting n = 0, and nothing raised.
  set.seed(43); n <- 80
  d <- sf::st_as_sf(data.frame(lon = runif(n, 9, 9.3), lat = runif(n, 48.6, 48.9),
                               a = rnorm(n)), coords = c("lon", "lat"), crs = 4326)
  d$z <- 2 * d$a + rnorm(n, 0, 0.2)
  # Only the units warning is asserted.  It is pure R and fires before the fit;
  # what a GWmodel version then DOES with windows that hold no neighbours (all
  # coefficients NaN, an `inv(): matrix is singular` message, an error) varies,
  # and the post-fit non-finite count reports it when it happens rather than
  # promising it will.
  # Collect every warning rather than nesting expect_warning(): how a given
  # GWmodel version reacts to windows holding no neighbours varies (all
  # coefficients NaN, an `inv(): matrix is singular` message, an error), so the
  # number of warnings is a backend property and only this one is ours.
  .warns <- function(expr) {
    w <- character(0)
    withCallingHandlers(
      suppressMessages(try(expr, silent = TRUE)),
      warning = function(cond) {
        w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
      })
    w
  }
  expect_true(any(grepl("ten-thousandth",
                        .warns(fit_gwr_model(d, "z", "a", bandwidth = 0.2,
                                             adaptive = FALSE)))))
  # A bandwidth in the units the fit actually runs in draws no such warning.
  expect_false(any(grepl("ten-thousandth",
                         .warns(fit_gwr_model(d, "z", "a", bandwidth = 5000,
                                              adaptive = FALSE)))))
})


test_that("non-finite local coefficients are counted and reported", {
  skip_if_not_installed("GWmodel")
  # The post-fit half of the same story, asserted without depending on how a
  # given GWmodel version fails: build the count from an SDF directly.  Nothing
  # used to report this at all -- fitted(), residuals(), summary() and
  # model_metrics() all drop non-finite rows, so a fit in which 182 of 200
  # local regressions failed reported n = 18 and R2 = 0.96 as if that were the
  # whole model.
  sdf <- data.frame(Intercept = c(1, NaN, 3, 4), a = c(1, 2, NaN, 4))
  ok  <- vapply(sdf[c("Intercept", "a")],
                function(z) is.finite(suppressWarnings(as.numeric(z))),
                logical(nrow(sdf)))
  expect_equal(sum(!apply(ok, 1L, all)), 2L)
})


test_that("a factor level with no training rows is an error, not a guess", {
  skip_if_not_installed("ranger")
  # .rf_align_levels() checked the level SET, not the levels that occur, so an
  # ordinary subset (or a CV fold holding out a whole class) let a level with
  # zero training rows through and ranger returned another level's prediction.
  set.seed(47); n <- 90
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               g = factor(rep(c("a", "b", "c"), length.out = n))),
                    coords = c("x", "y"), crs = 32632)
  d$z <- ifelse(d$g == "c", 100, 50) + rnorm(n)
  train <- d[d$g != "c", ]                        # level "c" is retained, unused
  expect_true("c" %in% levels(train$g))

  fit <- suppressMessages(fit_rf_model(train, "z", "g", num_trees = 60, seed = 1))
  expect_error(predict(fit, newdata = d[d$g == "c", ]),
               "never grown with")
})


test_that("a logical predictor supplied as text is refused", {
  skip_if_not_installed("ranger")
  # is.numeric / is.factor / is.character are all FALSE for a logical, so the
  # guard skipped it: a "TRUE"/"FALSE" character column (a CSV round trip) was
  # factor-coded 1/2 against splits built on 0/1, sending every row one way.
  set.seed(53); n <- 80
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               flag = rep(c(TRUE, FALSE), length.out = n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- ifelse(d$flag, 1.5, -1.75) + rnorm(n, 0, 0.2)
  fit <- suppressMessages(fit_rf_model(d, "z", "flag", num_trees = 60, seed = 1))

  txt <- d; txt$flag <- as.character(txt$flag)
  expect_error(predict(fit, newdata = txt), "was logical")
  # The honest coercion still works.
  txt2 <- d; txt2$flag <- as.numeric(txt2$flag)
  expect_equal(predict(fit, newdata = txt2), predict(fit, newdata = d),
               tolerance = 1e-10)
})


test_that("cv_rf's seed reaches the forest", {
  skip_if_not_installed("ranger")
  # fit_rf_model() carries seed = 123L and `...` did not override it, so every
  # fold's forest was grown with ranger seed 123 whatever the caller passed --
  # and no argument could change it.
  d <- .p5_pts(n = 90); d$z <- 2 * d$a + rnorm(90)
  f <- suppressMessages(make_folds(d, k = 3, method = "block_kfold", seed = 1))
  r1 <- suppressMessages(cv_rf(d, "z", "a", folds = f, seed = 1, num_trees = 60))
  r2 <- suppressMessages(cv_rf(d, "z", "a", folds = f, seed = 2, num_trees = 60))
  r3 <- suppressMessages(cv_rf(d, "z", "a", folds = f, seed = 1, num_trees = 60))
  expect_false(isTRUE(all.equal(r1$predictions$yhat, r2$predictions$yhat)))
  expect_equal(r1$predictions$yhat, r3$predictions$yhat)     # still reproducible
})


test_that("gwr_model_selection refuses candidates fit_gwr_model would reject", {
  skip_if_not_installed("GWmodel")
  set.seed(59); n <- 60
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               a = rnorm(n),
                               grp = factor(sample(letters[1:3], n, TRUE))),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + rnorm(n, 0, 0.2)
  expect_error(suppressMessages(gwr_model_selection(d, "z", c("a", "grp"))),
               "not numeric")
})


# --------------------------------------------------------------------------
# Metrics, printing, evaluation
# --------------------------------------------------------------------------

test_that("summary() applies the same response guard as model_metrics()", {
  d <- .p5_pts(n = 40)
  d$z <- 2 * d$a + rnorm(40)
  fit <- new_spatial_fit(
    subclass = "lmsurf_fit", engine = stats::lm(z ~ a, sf::st_drop_geometry(d)),
    formula = z ~ a, response_var = "z_chr", predictor_vars = "a",
    data_sf = transform(d, z_chr = as.character(round(d$z, 3))))
  expect_error(summary(fit), "not numeric")
  expect_error(model_metrics(fit), "not numeric")
})


test_that("print() shows one Formula line and the CRS", {
  d <- .p5_pts(n = 40); d$z <- 2 * d$a + rnorm(40)
  fit <- new_spatial_fit(
    subclass = "lmsurf_fit", engine = stats::lm(z ~ a, sf::st_drop_geometry(d)),
    formula = stats::as.formula(paste(
      "z ~ a + gp(..x, ..y, k = 21, c = 1.53561454721031, scale = FALSE,",
      "iso = FALSE)")),
    response_var = "z", predictor_vars = "a", data_sf = d)
  txt <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_equal(length(grep("^  Formula", strsplit(txt, "\n")[[1]])), 1L)
  expect_match(txt, "CRS", fixed = TRUE)
  expect_match(txt, "32632", fixed = TRUE)
})


test_that("weights with a non-zero diagonal are corrected, not accepted", {
  skip_if_not_installed("ranger")
  # Every moment formula assumes w_ii = 0.  A self-inclusive kNN matrix
  # rejected the null on 75% of white-noise residuals at a nominal 5%.
  set.seed(61); n <- 100
  xy <- cbind(runif(n, 0, 100), runif(n, 0, 100))
  d <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = rnorm(n), a = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  fit <- suppressMessages(fit_rf_model(d, "z", "a", num_trees = 60, seed = 1))

  D <- as.matrix(stats::dist(xy)); W <- matrix(0, n, n)
  for (i in seq_len(n)) W[i, order(D[i, ])[1:9]] <- 1   # includes self
  W <- W / rowSums(W)
  expect_warning(r <- residual_morans_i(fit, weights = W), "diagonal")
  expect_gt(r$p_value, 0.01)
})


test_that("compare_models_cv scores every model on identical folds", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("GWmodel")
  # With folds = NULL each cv_*() built its own from k/seed/block_size/..., so
  # a per-model block_size (or seed = NULL) silently gave the models different
  # splits and the table ranked models scored on different data.
  set.seed(67); n <- 150
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + rnorm(n, 0, 0.5)
  res <- suppressWarnings(suppressMessages(compare_models_cv(
    d, "z", "a", models = c("GWR", "RF"), k = 4, seed = 11,
    gwr_args = list(block_size = 150), rf_args = list(num_trees = 60))))

  fg <- res$cv_results$gwr_cv$folds
  fr <- res$cv_results$rf_cv$folds
  key <- function(f) lapply(f, function(z) sort(z$test))
  expect_equal(key(fg), key(fr))
})


test_that("forward selection has a null baseline", {
  skip_if_not_installed("ranger")
  # `best` started at Inf/-Inf with the stopping test gated on
  # is.finite(best), so the first step was accepted whatever its score.
  set.seed(71); n <- 90
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               n1 = rnorm(n), n2 = rnorm(n), n3 = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- rnorm(n)                                  # pure noise
  # lm(), not ranger: the null baseline is an INTERCEPT-ONLY fit, which
  # fit_rf_model() refuses ("a forest has nothing to split on").  The fallback
  # for a backend that cannot fit one is the old unconditional first step, so
  # the baseline has to be tested with a backend that can.
  fitf <- function(tr, vars) lm_spatial_fit(tr, "z", vars)

  sel <- suppressMessages(select_features_forward(
    d, "z", c("n1", "n2", "n3"), fitf, k = 3, seed = 1, quiet = TRUE))
  expect_true(any(sel$history$step == 0L))
  expect_equal(sel$history$variable[sel$history$step == 0L], "<none>")
  # An impossible tolerance now selects nothing at all.
  tight <- suppressMessages(select_features_forward(
    d, "z", c("n1", "n2", "n3"), fitf, k = 3, tol = 1e6, seed = 1, quiet = TRUE))
  expect_length(tight$selected, 0L)
})


test_that("determine_optimal_levels refuses a factor response", {
  set.seed(73); n <- 120
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z  <- 2 * d$a + rnorm(n)
  d$zf <- cut(d$z, 3, labels = c("low", "mid", "high"))
  expect_error(determine_optimal_levels(d, response_var = "zf",
                                        predictor_vars = "a", max_levels = 8),
               "must be numeric or logical")
})


# --------------------------------------------------------------------------
# predict_surface, AOA, plotting, logging
# --------------------------------------------------------------------------

test_that("predict_surface refuses a grid it cannot build", {
  skip_if_not_installed("ranger")
  d <- .p5_pts(n = 60); d$z <- 2 * d$a + rnorm(60)
  fit <- suppressMessages(fit_rf_model(d, "z", "a", num_trees = 40, seed = 1))
  expect_error(predict_surface(fit, cell_size = 0.005), "above the")
  expect_error(predict_surface(fit, n_cells = 1e10), "n_cells", fixed = TRUE)
})


test_that("predict_surface replays the CRS assumption on grid and boundary", {
  skip_if_not_installed("ranger")
  # .replay_crs_assumption() ran on `covariates` only, so a CRS-less grid whose
  # own bbox does not trip the heuristic had the fit's projected CRS STAMPED
  # onto raw degrees; st_nearest_feature() then gave every cell the same
  # covariate row and the surface collapsed to a constant.
  set.seed(79); n <- 80
  train <- sf::st_as_sf(data.frame(lon = runif(n, 9, 10), lat = runif(n, 48, 49)),
                        coords = c("lon", "lat"))       # CRS-less lon/lat
  train$w <- sf::st_coordinates(train)[, 1] * 10
  train$z <- 3 * train$w + rnorm(n, 0, 0.1)
  fit <- suppressWarnings(suppressMessages(
    fit_rf_model(train, "z", "w", num_trees = 60, seed = 1)))
  expect_identical(attr(fit$data_sf, "crs_assumed"), "EPSG:4326")

  grid <- sf::st_as_sf(expand.grid(lon = c(9.2, 9.5), lat = c(48.2, 48.5)),
                       coords = c("lon", "lat"))        # CRS-less, small bbox
  cov  <- train
  surf <- suppressWarnings(suppressMessages(
    predict_surface(fit, grid = grid, covariates = cov)))
  expect_gt(length(unique(round(surf$.pred, 6))), 1L)    # not a constant
})


test_that("the AOA variance test is relative, not absolute", {
  # An absolute floor on the raw sd (1.49e-8) dropped a predictor for being
  # recorded in the wrong UNIT, which flipped an obvious extrapolation into a
  # point comfortably inside the AOA.
  set.seed(83); n <- 60
  tr <- data.frame(a = rnorm(n), b = rnorm(n))
  nw <- data.frame(a = c(0, 0), b = c(0, 25))     # row 2 is far out in b
  sc  <- .aoa_scaling(as.matrix(tr))
  sc2 <- .aoa_scaling(as.matrix(transform(tr, b = b * 1e-9)))
  expect_true(all(sc$keep))
  expect_true(all(sc2$keep))                      # unit change is not "no variance"
  # A genuinely constant column is still dropped.
  sc3 <- .aoa_scaling(as.matrix(transform(tr, b = 0.1)))
  expect_false(sc3$keep[["b"]])
})


test_that("the AOA accepts a polygon newdata for a coordinate model", {
  skip_if_not_installed("ranger")
  set.seed(89); n <- 80
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + rnorm(n, 0, 0.2)
  fit <- suppressMessages(fit_rf_model(d, "z", "a", include_coords = TRUE,
                                       num_trees = 60, seed = 1))
  cells <- sf::st_sf(a = rnorm(9),
                     geometry = sf::st_make_grid(d, n = c(3, 3)))
  out <- suppressWarnings(area_of_applicability(cells, model = fit))
  expect_equal(nrow(out$aoa), 9L)
  expect_equal(length(out$aoa$DI), 9L)
  expect_equal(out$n_inside + out$n_outside + out$n_na, 9L)
})


test_that("plot_tessellation_map builds with a CRS-less overlay", {
  skip_if_not_installed("ggplot2")
  tess <- suppressMessages(build_tessellation(
    .p5_pts(n = 40), method = "voronoi"))
  overlay <- sf::st_as_sf(data.frame(x = c(100, 200), y = c(100, 200)),
                          coords = c("x", "y"))          # no CRS
  p <- suppressWarnings(plot_tessellation_map(tess$cells, features_sf = overlay))
  expect_silent(invisible(ggplot2::ggplot_build(p)))
})


test_that("spatialkit_quiet targets the console appender", {
  before <- logger::log_threshold(namespace = "spatialkit", index = 2)
  on.exit(logger::log_threshold(before, namespace = "spatialkit", index = 2),
          add = TRUE)
  spatialkit_quiet(TRUE)
  expect_identical(logger::log_threshold(namespace = "spatialkit", index = 2),
                   logger::FATAL)
  spatialkit_quiet(FALSE)
  expect_identical(logger::log_threshold(namespace = "spatialkit", index = 2),
                   logger::WARN)
  # Index 1 (the temp-file trace) is untouched.
  expect_identical(logger::log_threshold(namespace = "spatialkit", index = 1),
                   logger::INFO)
})
