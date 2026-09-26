# tests/testthat/test-kriging-adequacy.R
# ---------------------------------------------------------------------------
# kriging_adequacy(): per-cell block-kriging variance from a fitted
# variogram, its ratio to the cell's no-data variance, the comparison with
# s^2/n, the blocked cross-validation statistic of the kriging variance, and
# repeat measurements at one location.
# ---------------------------------------------------------------------------

ka_field <- function(n = 240, seed = 1, psill = 0.8, nugget = 0.2, a = 100) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(psill * exp(-d / a) + diag(nugget + 1e-8, n))) %*% rnorm(n))
  sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
}
ka_bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 32632))
ka_true_sac <- function() {
  vm <- gstat::vgm(psill = 0.8, "Exp", range = 100, nugget = 0.2)
  structure(300, class = "sac_range", variogram_model = vm)
}

test_that("the diagnostics are computed per cell and the CV statistic is near 1 with the true variogram", {
  skip_if_not_installed("gstat")
  pts <- ka_field()
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  ka <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), k = 4, seed = 1)
  expect_s3_class(ka, "kriging_adequacy")
  expect_s3_class(ka, "sf")
  expect_equal(nrow(ka), nrow(cells))
  df <- sf::st_drop_geometry(ka)
  expect_true(all(c("poly_id", "n", "mean", "se", "kr_pred", "kr_var", "kr_ratio",
                    "kr_exceeds_design", "kr_shift") %in% names(df)))
  expect_equal(sum(df$n), nrow(pts))
  expect_true(all(is.finite(df$kr_pred)))
  expect_true(all(df$kr_var >= 0))
  expect_true(all(df$kr_ratio >= 0 & df$kr_ratio <= 1))
  # A well-sampled cell's block variance is a small share of its no-data variance.
  expect_lt(stats::median(df$kr_ratio), 0.2)
  # The plain means agree with summarize_by_cell()'s.
  sm <- summarize_by_cell(asg, response_var = "z")
  expect_equal(df$mean[match(sm$poly_id, df$poly_id)], sm$resp_mean_z)
  expect_equal(df$se[match(sm$poly_id, df$poly_id)], sm$..se_resp_z)
  expect_equal(df$kr_exceeds_design, df$kr_var > df$se^2)
  expect_equal(df$kr_shift, (df$kr_pred - df$mean) / df$se)
  # The kriged estimate is the block-kriging prediction gstat returns.
  bk <- gstat::krige(z ~ 1, locations = pts, newdata = cells,
                     model = attr(ka_true_sac(), "variogram_model"), nmax = 50,
                     debug.level = 0)
  expect_equal(df$kr_pred, as.numeric(bk$var1.pred))
  expect_equal(df$kr_var, as.numeric(bk$var1.var))
  cv <- attr(ka, "cv")
  expect_equal(cv$k, 4L); expect_equal(cv$method, "block_kfold")
  expect_equal(cv$n_pred, nrow(pts))
  expect_gt(cv$zscore_var, 0.6); expect_lt(cv$zscore_var, 1.6)
  expect_equal(attr(ka, "sill"), 1); expect_equal(attr(ka, "nugget"), 0.2)
  expect_true(attr(ka, "range_identified"))
  expect_output(print(ka), "Block-kriging adequacy over 16 cells")
  expect_output(print(ka), "var of standardised error")
})

test_that("empty and one-point cells get a kriged estimate and NA design comparisons", {
  skip_if_not_installed("gstat")
  pts <- ka_field(n = 60, seed = 2)
  cells <- create_grid_polygons(ka_bnd, target_cells = 64, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  ka <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), k = 3, seed = 1)
  df <- sf::st_drop_geometry(ka)
  expect_gt(sum(df$n == 0L), 0L)
  expect_true(all(is.finite(df$kr_pred[df$n == 0L])))
  expect_true(all(is.na(df$kr_exceeds_design[df$n < 2L])))
  expect_true(all(is.na(df$kr_shift[df$n < 2L])))
  expect_true(all(is.na(df$mean[df$n == 0L])))
  # Empty cells are worse determined than full ones.
  expect_gt(mean(df$kr_ratio[df$n == 0L]), mean(df$kr_ratio[df$n >= 3L]))
})

test_that("the variogram is estimated when not supplied, and refused when absent or of an unknown family", {
  skip_if_not_installed("gstat")
  pts <- ka_field(n = 200, seed = 3)
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  ka <- suppressWarnings(kriging_adequacy(asg, "z", cells, k = 3, seed = 1))
  vm <- attr(ka, "variogram")
  expect_s3_class(vm, "data.frame")
  expect_true(all(as.character(vm$model) %in% c("Nug", "Exp", "Sph")))
  expect_equal(attr(ka, "sill"), sum(vm$psill))
  # No model at all.
  expect_error(kriging_adequacy(asg, "z", cells, sac = structure(NA_real_, class = "sac_range")),
               "carries no variogram model")
  # An unsupported family.
  bad <- structure(300, class = "sac_range",
                   variogram_model = gstat::vgm(psill = 1, "Mat", range = 100, nugget = 0.1, kappa = 1.5))
  expect_error(kriging_adequacy(asg, "z", cells, sac = bad), "family .Mat.")
  # An unidentified range is used with a warning.
  unid <- structure(NA_real_, class = "sac_range",
                    variogram_model = attr(ka_true_sac(), "variogram_model"),
                    rejected_reason = "fitted range exceeds the largest lag fitted")
  expect_warning(ka2 <- kriging_adequacy(asg, "z", cells, sac = unid, k = 3), "range was not identified")
  expect_false(attr(ka2, "range_identified"))
})

test_that("supplied folds are honoured and inputs are validated", {
  skip_if_not_installed("gstat")
  pts <- ka_field(n = 200, seed = 4)
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  f <- make_folds(asg, k = 3, method = "random_kfold", seed = 2)
  ka <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = f)
  expect_equal(attr(ka, "cv")$method, "random_kfold")
  expect_equal(attr(ka, "cv")$k, 3L)
  lab <- sample(1:4, nrow(asg), replace = TRUE)
  ka_lab <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = lab)
  expect_equal(attr(ka_lab, "cv")$k, 4L)
  expect_error(kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = 1:5),
               "vector of fold labels with one entry per row")
  expect_error(kriging_adequacy(asg, "nope", cells), "`response_var` must name a column")
  expect_error(kriging_adequacy(sf::st_drop_geometry(asg), "z", cells), "must be an sf object")
  expect_error(kriging_adequacy(asg, "z", pts), "POLYGON, MULTIPOLYGON")
  no_id <- asg; no_id$poly_id <- NULL
  expect_error(kriging_adequacy(no_id, "z", cells), "could not find a cell ID column")
  expect_error(kriging_adequacy(asg, "z", cells, nmax = 0), "`nmax` must be")
  # Rows with a missing response are dropped with a log line, not an error.
  asg2 <- asg; asg2$z[1:5] <- NA
  lines <- capture_spatialkit_log(ka2 <- kriging_adequacy(asg2, "z", cells, sac = ka_true_sac(), k = 3))
  expect_true(log_has(lines, "dropping 5 point"))
  expect_equal(sum(sf::st_drop_geometry(ka2)$n), nrow(asg) - 5L)
  expect_equal(attr(ka2, "n_points"), nrow(asg) - 5L)
})


test_that("print() survives a subset that no longer carries the fitted summary", {
  skip_if_not_installed("gstat")
  pts   <- ka_field(n = 200)
  cells <- create_grid_polygons(ka_bnd, target_cells = 12, type = "square")
  asg   <- assign_features_to_polygons(pts, cells)
  ka    <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), k = 3, seed = 1)

  expect_output(print(ka), "^Block-kriging adequacy over")
  # `[` on a data frame keeps the class and drops the attributes, and knitr
  # reaches print() through knit_print.data.frame() without being asked, so a
  # subset in a knitted document used to abort it with "argument is of length
  # zero" from is.finite(NULL).
  for (sub in list(ka[1:3, ],
                   sf::st_drop_geometry(ka)[, 1:4],
                   sf::st_drop_geometry(ka)[, c("n", "mean")])) {
    expect_error(utils::capture.output(print(sub)), NA)
  }
  expect_output(print(sf::st_drop_geometry(ka)[, 1:4]), "subset")
  # A result whose cross-validation was not computed prints everything else.
  bare <- ka
  attr(bare, "cv") <- list()
  expect_output(print(bare), "blocked CV: not computed")
})


# gstat's own variance of a cell mean with no data, C(B,B) on its own
# discretisation and nugget handling: simple kriging from one datum so far
# away that its covariance with every cell is exactly zero.
ka_gstat_prior <- function(cells, vm) {
  far <- sf::st_sf(z = 0, geometry = sf::st_sfc(sf::st_point(c(1e9, 1e9)),
                                                crs = sf::st_crs(cells)))
  as.numeric(gstat::krige(z ~ 1, far, cells, model = vm, beta = 0, debug.level = 0)$var1.var)
}

test_that("kr_ratio is over each cell's no-data variance, so a cell the data do not reach reads 1", {
  skip_if_not_installed("gstat")
  # kr_var is the variance of a cell MEAN; it was divided by the point sill,
  # which a cell mean never reaches, so empty cells 130-410 m beyond a 90 m
  # effective range read about 0.1 and print() said no cell was above 0.5.
  set.seed(11); n <- 240
  x <- runif(n, 0, 500); y <- runif(n, 0, 1000)        # the western half only
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(0.8 * exp(-d / 30) + diag(0.2 + 1e-8, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  vm <- gstat::vgm(psill = 0.8, "Exp", range = 30, nugget = 0.2)
  ka <- kriging_adequacy(asg, "z", cells, sac = structure(90, class = "sac_range",
                                                          variogram_model = vm),
                         k = 4, seed = 1)
  df <- sf::st_drop_geometry(ka)
  empty <- df$n == 0L
  expect_equal(sum(empty), 8L)
  expect_true(all(df$kr_var[empty] < 0.2))              # far below the point sill of 1
  expect_true(all(df$kr_ratio[empty] > 0.99))
  expect_true(all(df$kr_ratio[!empty] < min(df$kr_ratio[empty])))
  # The denominator is the cell's C(B,B) as gstat block-kriges it.
  prior <- ka_gstat_prior(cells, vm)
  expect_equal(df$kr_ratio, pmin(df$kr_var / prior, 1), tolerance = 1e-5)
  expect_output(print(ka), sprintf("no-data variance of the cell mean .* %d cell\\(s\\) above 0.5",
                                   sum(df$kr_ratio > 0.5)))
  expect_gte(sum(df$kr_ratio > 0.5), 8L)
})

test_that("a large empty cell ranks above small populated ones on kr_ratio", {
  skip_if_not_installed("gstat")
  # With unequal cells, as Voronoi and Delaunay tessellations make them, the
  # ratio over the sill ranked a 600 x 1000 m cell with no data below
  # populated 100 m cells, because a big cell's mean varies little.
  set.seed(4); n <- 300
  x <- runif(n, 0, 400); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(0.8 * exp(-d / 100) + diag(0.2 + 1e-8, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
  sq <- function(x0, x1, y0, y1) sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1),
                                                           c(x0, y1), c(x0, y0))))
  small <- sf::st_make_grid(sf::st_sfc(sq(0, 400, 0, 1000), crs = 32632), cellsize = 100)
  cells <- sf::st_sf(poly_id = seq_len(length(small) + 1L),
                     geometry = c(small, sf::st_sfc(sq(400, 1000, 0, 1000), crs = 32632)))
  asg <- assign_features_to_polygons(pts, cells)
  df <- sf::st_drop_geometry(kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), k = 4, seed = 1))
  big <- df$poly_id == nrow(cells)
  expect_equal(df$n[big], 0L)
  expect_gt(df$kr_ratio[big], 0.99)
  expect_true(all(df$kr_ratio[!big] < 0.5))
  expect_gt(df$kr_ratio[big], max(df$kr_ratio[!big]))
})


# Stations visited three times: a smooth field sampled at 80 sites, each visit
# with its own measurement error of variance `me_var`.
ka_revisits <- function(me_var, seed = 5) {
  st <- ka_field(n = 80, seed = seed, nugget = 0)
  set.seed(seed + 1)
  v <- rbind(st, st, st)
  v$z <- v$z + stats::rnorm(nrow(v), sd = sqrt(me_var))
  v
}

test_that("repeat visits to a station are kriged from their means instead of coming back NA", {
  skip_if_not_installed("gstat")
  # Two observations at one location get the full sill as their covariance in
  # gstat, so every kriging system holding a pair was singular: 80 stations x
  # 3 visits gave NA in all 16 cells and no CV prediction, with no warning.
  rv <- ka_revisits(me_var = 0.4)        # visits differ by more than the 0.2 nugget
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(rv, cells)
  vm <- attr(ka_true_sac(), "variogram_model")
  expect_warning(ka <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), k = 4,
                                        seed = 1, nmax = 1000),
                 "240 point\\(s\\) share a location .* 80 distinct locations")
  df <- sf::st_drop_geometry(ka)
  expect_true(all(is.finite(df$kr_pred)) && all(is.finite(df$kr_var)))
  expect_equal(sum(df$n), 240L)                           # plain means use every visit
  expect_equal(attr(ka, "n_points"), 240L)
  expect_equal(attr(ka, "n_locations"), 80L)
  expect_equal(attr(ka, "cv")$n_pred, 80L)
  expect_true(is.finite(attr(ka, "cv")$zscore_var))
  # With the whole nugget differing between visits, kriging the 80 means is
  # kriging all 240 visits with the nugget as their measurement error.
  bk <- gstat::krige(z ~ 1, rv, cells, model = vm[vm$model != "Nug", ],
                     weights = rep(1 / 0.2, nrow(rv)), debug.level = 0)
  expect_equal(df$kr_pred, as.numeric(bk$var1.pred), tolerance = 1e-8)
  expect_equal(df$kr_var, as.numeric(bk$var1.var), tolerance = 1e-8)
  expect_output(print(ka), "240 points at 80 distinct locations")
  expect_output(print(ka), "4 folds, 80 locations")
})

test_that("identical repeat records change nothing the kriging reports", {
  skip_if_not_installed("gstat")
  # A mean of identical replicates is one observation, nugget and all.
  pts <- ka_field(n = 120, seed = 7)
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  lab <- rep(1:4, length.out = nrow(pts))
  one <- kriging_adequacy(assign_features_to_polygons(pts, cells), "z", cells,
                          sac = ka_true_sac(), folds = lab)
  expect_warning(two <- kriging_adequacy(assign_features_to_polygons(rbind(pts, pts), cells),
                                         "z", cells, sac = ka_true_sac(), folds = c(lab, lab)),
                 "share a location")
  a <- sf::st_drop_geometry(one); b <- sf::st_drop_geometry(two)
  expect_equal(b$n, 2L * a$n)
  # To 1e-6: gstat's block variance from a model nugget and from the same
  # variance passed as a measurement error agree to about 1e-7, not to 1e-15.
  expect_equal(b[c("mean", "kr_pred", "kr_var", "kr_ratio")], a[c("mean", "kr_pred", "kr_var", "kr_ratio")],
               tolerance = 1e-6)
  expect_equal(attr(two, "cv")[c("zscore_var", "zscore_mean", "rmse", "n_pred")],
               attr(one, "cv")[c("zscore_var", "zscore_mean", "rmse", "n_pred")], tolerance = 1e-6)
})

test_that("a kriging system gstat cannot solve raises a warning and print() counts it", {
  skip_if_not_installed("gstat")
  # gstat answers a singular system with NA, and at debug.level 0 says
  # nothing; print() then claimed an estimate for every empty cell.  Points a
  # micrometre from ten others under a nugget-free Gaussian model are singular
  # without sharing a location.
  pts <- ka_field(n = 60, seed = 3)
  xy <- sf::st_coordinates(pts)
  twin <- sf::st_as_sf(data.frame(x = xy[1:10, 1] + 1e-6, y = xy[1:10, 2], z = pts$z[1:10] + 0.1),
                       coords = c("x", "y"), crs = 32632)
  cells <- create_grid_polygons(ka_bnd, target_cells = 64, type = "square")
  asg <- assign_features_to_polygons(rbind(pts, twin), cells)
  sac <- structure(170, class = "sac_range", variogram_model = gstat::vgm(1, "Gau", 100))
  w <- character()
  ka <- withCallingHandlers(
    kriging_adequacy(asg, "z", cells, sac = sac, k = 4, seed = 1),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  df <- sf::st_drop_geometry(ka)
  n_na <- sum(!is.finite(df$kr_pred))
  expect_gt(n_na, 0L)
  expect_true(any(grepl(sprintf("returned no estimate for %d of 64 cell", n_na), w)))
  expect_true(any(grepl("no cross-validation prediction it could use", w)))
  expect_output(print(ka), sprintf("no kriged estimate for %d of the 64 cells", n_na))
  n_empty_kr <- sum(df$n == 0L & is.finite(df$kr_pred))
  expect_lt(n_empty_kr, sum(df$n == 0L))
  expect_output(print(ka), sprintf("available for %d of them", n_empty_kr))
})



# ---------------------------------------------------------------------------
# Round 2 of the review.
# ---------------------------------------------------------------------------

# The standardised errors of kriging each split's held-out points from that
# split's own training rows, by hand.
ka_manual_cv <- function(pts, splits, vm, nmax = 50) {
  unlist(lapply(splits, function(s) {
    kf <- gstat::krige(z ~ 1, pts[s$train, ], pts[s$test, ], model = vm,
                       nmax = nmax, debug.level = 0)
    (pts$z[s$test] - kf$var1.pred) / sqrt(kf$var1.var)
  }))
}

test_that("buffered leave-one-out and NNDM folds keep the points they exclude out of the kriging", {
  skip_if_not_installed("gstat")
  # Their exclusion zones live only in each split's train set, and the fold
  # labels passed on were 1:n, so both ran as plain leave-one-out and
  # print() called it blocked CV.
  vm <- attr(ka_true_sac(), "variogram_model")
  pts <- ka_field(n = 120, seed = 9)
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  fb <- suppressWarnings(make_folds(asg, k = 1, method = "buffered_loo", buffer = 250, seed = 1))
  expect_lt(mean(lengths(lapply(fb$folds, `[[`, "train"))), nrow(asg) - 1)
  ka <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = fb)
  cv <- attr(ka, "cv")
  zs <- ka_manual_cv(asg, fb$folds, vm)
  expect_equal(cv$zscore_var, stats::var(zs), tolerance = 1e-10)
  expect_equal(cv$n_pred, nrow(asg))
  expect_equal(cv$k, nrow(asg))
  expect_identical(cv$method, "buffered_loo")
  loo <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = seq_len(nrow(asg)))
  expect_gt(abs(cv$rmse - attr(loo, "cv")$rmse), 0.01)
  expect_output(print(ka), "buffered leave-one-out CV (buffered_loo, 120 folds", fixed = TRUE)
  expect_output(print(loo), "  CV (supplied labels", fixed = TRUE)

  # NNDM on a clustered sample predicted over the whole square.
  set.seed(13)
  cx <- runif(8, 100, 900); cy <- runif(8, 100, 900)
  x <- pmin(pmax(rep(cx, each = 15) + rnorm(120, sd = 30), 1), 999)
  y <- pmin(pmax(rep(cy, each = 15) + rnorm(120, sd = 30), 1), 999)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(0.8 * exp(-d / 100) + diag(0.2 + 1e-8, 120))) %*% rnorm(120))
  cl <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
  acl <- assign_features_to_polygons(cl, cells)
  pp <- sf::st_as_sf(sf::st_make_grid(ka_bnd, n = 12, what = "centers"))
  fn <- suppressWarnings(make_folds(acl, k = 1, method = "nndm", prediction_points = pp, seed = 1))
  skip_if(sum(lengths(lapply(fn$folds, `[[`, "train"))) == nrow(acl) * (nrow(acl) - 1L),
          "NNDM excluded no neighbours on this draw")
  kn <- kriging_adequacy(acl, "z", cells, sac = ka_true_sac(), folds = fn)
  expect_equal(attr(kn, "cv")$zscore_var, stats::var(ka_manual_cv(acl, fn$folds, vm)),
               tolerance = 1e-10)
  expect_output(print(kn), "NNDM leave-one-out CV (nndm", fixed = TRUE)
})

test_that("the points are kriged in the CRS the variogram was fitted in", {
  skip_if_not_installed("gstat")
  # The range is a length in attr(sac, "crs"); a sac fitted in metres met
  # points in km and was read as 100 km.
  pts <- ka_field(n = 150, seed = 8)
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  sac <- ka_true_sac(); attr(sac, "crs") <- sf::st_crs(32632)
  lab <- rep(1:4, length.out = nrow(asg))
  m <- kriging_adequacy(asg, "z", cells, sac = sac, folds = lab)
  km <- sf::st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")
  k <- kriging_adequacy(sf::st_transform(asg, km), "z", sf::st_transform(cells, km),
                        sac = sac, folds = lab)
  a <- sf::st_drop_geometry(m); b <- sf::st_drop_geometry(k)
  expect_equal(b$kr_pred, a$kr_pred, tolerance = 1e-6)
  expect_equal(b$kr_var, a$kr_var, tolerance = 1e-6)
  expect_equal(attr(k, "cv")$zscore_var, attr(m, "cv")$zscore_var, tolerance = 1e-6)
  expect_equal(sf::st_crs(k), sf::st_crs(32632))
})

test_that("a variogram of residuals is used with a warning that it understates the variance", {
  skip_if_not_installed("gstat")
  pts <- ka_field(n = 120, seed = 10)
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  sac <- ka_true_sac()
  attr(sac, "detrended") <- TRUE; attr(sac, "detrend_method") <- "ols"
  lab <- rep(1:4, length.out = nrow(asg))
  expect_warning(kriging_adequacy(asg, "z", cells, sac = sac, folds = lab),
                 "variogram of the residuals on predictors .*understate")
  attr(sac, "detrended") <- FALSE
  expect_no_warning(kriging_adequacy(asg, "z", cells, sac = sac, folds = lab))
})

test_that("why a range was not identified is kept on the result and said as it is", {
  skip_if_not_installed("gstat")
  # Every refusal was explained as a sill never reached, and only a TRUE/FALSE
  # was kept, so a non-converged fit could not be told from the others.
  pts <- ka_field(n = 120, seed = 10)
  cells <- create_grid_polygons(ka_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  lab <- rep(1:4, length.out = nrow(asg))
  unid <- structure(NA_real_, class = "sac_range",
                    variogram_model = attr(ka_true_sac(), "variogram_model"),
                    rejected_reason = "variogram model did not converge")
  w <- character()
  ka <- withCallingHandlers(kriging_adequacy(asg, "z", cells, sac = unid, folds = lab),
                            warning = function(cnd) {
                              w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
                            })
  expect_true(any(grepl("range was not identified \\(variogram model did not converge\\)", w)))
  expect_false(any(grepl("sill was never reached", w)))
  expect_identical(attr(ka, "rejected_reason"), "variogram model did not converge")
  expect_output(print(ka), "range not identified (variogram model did not converge)", fixed = TRUE)
  ok <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = lab)
  expect_identical(attr(ok, "rejected_reason"), NA_character_)
})

test_that("a cell holding more locations than nmax is kriged from all of them, not its middle", {
  skip_if_not_installed("gstat")
  # gstat takes the nmax locations nearest a block's centre, so a cell with
  # more than nmax points was kriged from its central few.
  vm <- attr(ka_true_sac(), "variogram_model")
  pts <- ka_field(n = 300, seed = 12)
  sq <- function(x0, x1, y0, y1) sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1),
                                                           c(x0, y1), c(x0, y0))))
  small <- sf::st_make_grid(sf::st_sfc(sq(600, 1000, 0, 1000), crs = 32632), cellsize = 200)
  cells <- sf::st_sf(poly_id = seq_len(length(small) + 1L),
                     geometry = c(sf::st_sfc(sq(0, 600, 0, 1000), crs = 32632), small))
  asg <- assign_features_to_polygons(pts, cells)
  lab <- rep(1:3, length.out = nrow(asg))
  ka <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = lab, nmax = 30)
  df <- sf::st_drop_geometry(ka)
  expect_gt(df$n[1], 30L)
  expect_equal(df$kr_n_used[1], df$n[1] + 30L)
  # Its own points plus the 30 nearest its centre outside it, all used.
  xy <- sf::st_coordinates(asg)
  cen <- sp::coordinates(sf::as_Spatial(sf::st_geometry(cells)[1]))
  own <- asg$poly_id == 1L
  d <- sqrt((xy[, 1] - cen[1, 1])^2 + (xy[, 2] - cen[1, 2])^2); d[own] <- Inf
  sel <- c(which(own), order(d)[1:30])
  bk <- gstat::krige(z ~ 1, asg[sel, ], cells[1, ], model = vm, debug.level = 0)
  expect_equal(df$kr_pred[1], as.numeric(bk$var1.pred), tolerance = 1e-10)
  expect_equal(df$kr_var[1], as.numeric(bk$var1.var), tolerance = 1e-10)
  mid <- gstat::krige(z ~ 1, asg, cells[1, ], model = vm, nmax = 30, debug.level = 0)
  expect_gt(abs(df$kr_pred[1] - mid$var1.pred), 1e-3)
  # Cells whose own points all lie among the 30 nearest their centre keep
  # gstat's own neighbourhood.
  keep <- df$kr_n_used == 30L
  expect_gt(sum(keep), 5L)
  gs <- gstat::krige(z ~ 1, asg, cells[which(keep), ], model = vm, nmax = 30, debug.level = 0)
  expect_equal(df$kr_pred[keep], as.numeric(gs$var1.pred), tolerance = 1e-10)
  # A system above max_neighbours is left out with a warning, and counted.
  expect_warning(small_cap <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = lab,
                                               nmax = 30, max_neighbours = 100),
                 "1 cell\\(s\\) left out .*max_neighbours = 100")
  sc <- sf::st_drop_geometry(small_cap)
  expect_true(is.na(sc$kr_pred[1]) && is.na(sc$kr_n_used[1]) && is.na(sc$kr_ratio[1]))
  expect_equal(sc$kr_pred[-1], df$kr_pred[-1])
  expect_equal(attr(small_cap, "cells_left_out")[["size"]], 1L)
  expect_output(print(small_cap), "no kriged estimate for 1 of the 11 cells \\(1 left out as holding")
  expect_error(kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), max_neighbours = 0),
               "`max_neighbours` must be")
})

test_that("a cell whose bounding box dwarfs its area is left out instead of discretised", {
  skip_if_not_installed("gstat")
  # gstat lays its 500-point grid over a cell's whole bounding box, so a
  # sliver's memory grows with box / area: +592 MB at 7,072.
  sq <- sf::st_polygon(list(rbind(c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0))))
  strip <- sf::st_buffer(sf::st_linestring(rbind(c(0, 0), c(1000, 1000))), 5)
  rest <- sf::st_cast(sf::st_sfc(sf::st_difference(sq, strip)), "POLYGON")
  cells <- sf::st_sf(poly_id = seq_len(length(rest) + 1L),
                     geometry = sf::st_sfc(c(list(sf::st_intersection(strip, sq)), as.list(rest)),
                                           crs = 32632))
  ratio <- spatialkit:::.cell_box_ratio(cells)
  expect_gt(ratio[1], 50); expect_lt(ratio[1], 1000)
  expect_true(all(ratio[-1] < 3))
  pts <- ka_field(n = 150, seed = 14)
  asg <- assign_features_to_polygons(pts, cells)
  lab <- rep(1:3, length.out = nrow(asg))
  # Under the default bound the strip is kriged like any other cell ...
  all_in <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = lab)
  expect_true(all(is.finite(sf::st_drop_geometry(all_in)$kr_pred)))
  expect_equal(attr(all_in, "cells_left_out"), c(shape = 0L, size = 0L))
  # ... and above a tighter one it is left out, with a warning and a count.
  expect_warning(ka <- kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), folds = lab,
                                        max_box_ratio = 50),
                 "1 cell\\(s\\) left out .*max_box_ratio = 50")
  df <- sf::st_drop_geometry(ka)
  expect_true(is.na(df$kr_pred[1]) && is.na(df$kr_var[1]) && is.na(df$kr_ratio[1]))
  expect_equal(df$kr_pred[-1], sf::st_drop_geometry(all_in)$kr_pred[-1])
  expect_equal(attr(ka, "cells_left_out")[["shape"]], 1L)
  expect_output(print(ka), "1 left out as too thin or scattered")
  # Parts far apart and a cell that is mostly hole score the same way.
  far <- sf::st_multipolygon(list(
    list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0))),
    list(rbind(c(990, 0), c(1000, 0), c(1000, 10), c(990, 10), c(990, 0)))))
  holed <- sf::st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)),
                               rbind(c(5, 5), c(5, 95), c(95, 95), c(95, 5), c(5, 5))))
  r2 <- spatialkit:::.cell_box_ratio(sf::st_sf(geometry = sf::st_sfc(far, holed, crs = 32632)))
  expect_equal(r2, c(1000 * 10 / 200, 100^2 / (100^2 - 90^2)))
  expect_error(kriging_adequacy(asg, "z", cells, sac = ka_true_sac(), max_box_ratio = 0.5),
               "`max_box_ratio` must be")
})
