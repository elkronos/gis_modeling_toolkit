# tests/testthat/test-kriging-adequacy.R
# ---------------------------------------------------------------------------
# kriging_adequacy(): per-cell block-kriging variance from a fitted
# variogram, its ratio to the cell's no-data variance, the comparison with
# s^2/n, and the blocked cross-validation statistic of the kriging variance.
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

