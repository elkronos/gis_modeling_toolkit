# tests/testthat/test-plotting-diagnostics.R
# ---------------------------------------------------------------------------
# The curve behind the chosen point: plot_cv_metrics(), plot.aoa(),
# plot_calibration(), the response overlay on plot.spatial_fit(type =
# "variogram"), and the sweep plots plot.resolution_profile(),
# plot.feature_selection(), plot.gwr_model_selection().
#
# ggplot objects are built lazily; ggplot2::ggplot_build() forces the layer
# computation, and the assertions read the built data or the labels.
# ---------------------------------------------------------------------------

diag_points <- function(n = 240, seed = 11) {
  set.seed(seed)
  pts <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n), c = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  xy <- sf::st_coordinates(pts)
  pts$z <- 2 * pts$a - pts$b + sin(xy[, 1] / 200) + cos(xy[, 2] / 250) + rnorm(n, 0, 0.5)
  pts
}

layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], character(1))

# A cv result from the lm stand-in: no optional backend needed.
diag_cv <- function(pts, k = 4, metrics = NULL, p = NULL) {
  fit_fn <- function(train_sf) lm_spatial_fit(train_sf, "z", c("a", "b"))
  suppressMessages(cv_spatial(pts, "z", c("a", "b"), fit_fn = fit_fn, k = k,
                              seed = 1, metrics = metrics, p = p))
}


# ---- 9.1 -----------------------------------------------------------------

test_that("plot_cv_metrics draws one point per fold and the pooled line", {
  skip_if_not_installed("ggplot2")
  pts <- diag_points()
  cv  <- diag_cv(pts, k = 4)
  p   <- plot_cv_metrics(cv, "RMSE")
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  geoms <- layer_geoms(p)
  expect_true("GeomHline" %in% geoms)
  pt <- b$data[[which(geoms == "GeomPoint")]]
  expect_equal(nrow(pt), nrow(cv$fold_metrics))
  expect_equal(sort(pt$y), sort(cv$fold_metrics$RMSE))
  hl <- b$data[[which(geoms == "GeomHline")]]
  expect_equal(hl$yintercept, cv$overall$RMSE)
  expect_match(p$labels$caption, "pooled value")
  # A user metric column is drawn like any other.
  cv2 <- diag_cv(pts, k = 4, metrics = function(y, yhat) c(MedAE = stats::median(abs(y - yhat))))
  p2  <- plot_cv_metrics(cv2, "MedAE")
  b2  <- ggplot2::ggplot_build(p2)
  expect_equal(b2$data[[which(layer_geoms(p2) == "GeomHline")]]$yintercept, cv2$overall$MedAE)
})

test_that("plot_cv_metrics refuses an all-NA column with the reason, and names unknown columns", {
  skip_if_not_installed("ggplot2")
  pts <- diag_points()
  cv  <- diag_cv(pts, k = 3)
  expect_error(plot_cv_metrics(cv, "Adj_R2"), "NA in every fold")
  expect_error(plot_cv_metrics(cv, "Adj_R2"), "unless `p` was passed")
  expect_error(plot_cv_metrics(cv, "nope"), "not a column of the per-fold metrics")
  expect_error(plot_cv_metrics(cv, "nope"), "Available: .*RMSE")
  expect_error(plot_cv_metrics(list(a = 1), "RMSE"), "must be the list returned by")
  expect_error(plot_cv_metrics(cv, c("RMSE", "MAE")), "single column name")
  # With p supplied Adj_R2 is finite per fold and draws.
  cv_p <- diag_cv(pts, k = 3, p = 2L)
  expect_no_error(ggplot2::ggplot_build(plot_cv_metrics(cv_p, "Adj_R2")))
  # A per-fold extra with no pooled counterpart draws without the line and
  # says so.
  cv$fold_metrics$bandwidth <- c(40, 42, 41)
  pb <- plot_cv_metrics(cv, "bandwidth")
  expect_false("GeomHline" %in% layer_geoms(pb))
  expect_match(pb$labels$caption, "No pooled value")
})

test_that("plot_cv_metrics facets a compare_models_cv() result by model", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ranger")
  pts <- diag_points()
  cmp <- suppressMessages(suppressWarnings(
    compare_models_cv(pts, "z", c("a", "b"), models = "RF", k = 3,
                      rf_args = list(num_trees = 60), quiet = TRUE)))
  p <- plot_cv_metrics(cmp, "MAE")
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  hl <- b$data[[which(layer_geoms(p) == "GeomHline")]]
  expect_equal(hl$yintercept, cmp$overall$MAE)
  # A two-model comparison, faked from one, gets two panels.
  cmp2 <- cmp
  cmp2$by_fold <- rbind(cmp$by_fold, transform(cmp$by_fold, model = "GWR", MAE = MAE * 2))
  cmp2$overall <- rbind(cmp$overall, transform(cmp$overall, model = "GWR", MAE = MAE * 2))
  p2 <- plot_cv_metrics(cmp2, "MAE")
  expect_s3_class(p2$facet, "FacetWrap")
  b2 <- ggplot2::ggplot_build(p2)
  expect_equal(length(unique(b2$data[[which(layer_geoms(p2) == "GeomPoint")]]$PANEL)), 2L)
})

test_that("plot_cv_metrics refuses a result with no surviving fold", {
  skip_if_not_installed("ggplot2")
  pts <- diag_points(n = 60)
  cv <- suppressWarnings(suppressMessages(
    cv_spatial(pts, "z", "a", fit_fn = function(train_sf) stop("nope"), k = 3, seed = 1)))
  expect_error(plot_cv_metrics(cv, "RMSE"), "no per-fold metrics")
})


# ---- 9.2 -----------------------------------------------------------------

test_that("plot.aoa draws both distributions and the threshold", {
  skip_if_not_installed("ggplot2")
  pts <- diag_points()
  set.seed(12)
  new <- sf::st_as_sf(
    data.frame(x = runif(100, 0, 1000), y = runif(100, 0, 1000),
               a = rnorm(100, 1.5), b = rnorm(100)),
    coords = c("x", "y"), crs = 32632)
  aoa <- area_of_applicability(new, train_sf = pts, predictor_vars = c("a", "b"),
                               folds = make_folds(pts, k = 4, method = "block_kfold", seed = 1))
  p <- plot(aoa)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  geoms <- layer_geoms(p)
  expect_true("GeomVline" %in% geoms)
  expect_equal(b$data[[which(geoms == "GeomVline")]]$xintercept, aoa$threshold)
  # Two ECDF curves: one per set.
  ecdf_layer <- b$data[[which(geoms == "GeomStep")]]
  expect_equal(length(unique(ecdf_layer$group)), 2L)
  expect_match(p$labels$subtitle, sprintf("^%d of %d prediction locations outside", aoa$n_outside, aoa$n_new))
  expect_match(p$labels$caption, "cross-validated over block_kfold folds")

  ph <- plot(aoa, type = "histogram")
  expect_no_error(ggplot2::ggplot_build(ph))
  expect_true(all(c("GeomBar", "GeomPath", "GeomVline") %in% layer_geoms(ph)))

  # No folds: the caption says the threshold is optimistic.
  aoa0 <- area_of_applicability(new, train_sf = pts, predictor_vars = c("a", "b"))
  expect_match(plot(aoa0)$labels$caption, "not cross-validated")
  expect_error(plot.aoa(list(a = 1)), "must be the object returned by area_of_applicability")
})


# ---- 9.3 -----------------------------------------------------------------

test_that("the response variogram is overlaid on the residual variogram", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  # Response with strong short-range structure; the predictor explains part.
  set.seed(5); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  field <- as.numeric(t(chol(exp(-d / 50) + diag(1e-4, n))) %*% rnorm(n))
  w <- rnorm(n)
  pts <- sf::st_as_sf(data.frame(x = x, y = y, w = w, z = field + 2 * w),
                      coords = c("x", "y"), crs = 3857)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  p_on  <- plot(fit, type = "variogram")
  p_off <- plot(fit, type = "variogram", response = FALSE)
  expect_no_error(ggplot2::ggplot_build(p_on))
  g_on  <- layer_geoms(p_on); g_off <- layer_geoms(p_off)
  # Two more layers with the overlay: the hollow points and the dashed line.
  expect_equal(sum(g_on == "GeomPoint"), sum(g_off == "GeomPoint") + 1L)
  expect_gte(sum(g_on == "GeomLine"), sum(g_off == "GeomLine"))
  expect_null(p_off$labels$caption)
  expect_match(p_on$labels$caption, "Hollow points, dashed line: Response \\(z\\)")
  # The overlay is the response's own variogram: same lags, larger
  # semivariance than the residuals' (the model absorbed the 2 * w part).
  b <- ggplot2::ggplot_build(p_on)
  pt_layers <- which(g_on == "GeomPoint")
  resid_pts <- b$data[[pt_layers[1L]]]; resp_pts <- b$data[[pt_layers[2L]]]
  expect_equal(nrow(resid_pts), nrow(resp_pts))
  expect_gt(mean(resp_pts$y), mean(resid_pts$y))
  # Both ranges identified here, so the caption compares the sills.
  sac_res <- estimate_sac_range(sf::st_sf(r = residuals(fit), geometry = sf::st_geometry(pts)), "r")
  sac_z   <- estimate_sac_range(pts, "z")
  if (is.finite(sac_res) && is.finite(sac_z))
    expect_match(p_on$labels$caption, "Residual sill is [0-9]+% of the response sill")
  else
    expect_match(p_on$labels$caption, "Sills not compared")
  # The subtitle (the residual range) is untouched by the overlay.
  expect_identical(p_on$labels$subtitle, p_off$labels$subtitle)
})


# ---- 9.4 -----------------------------------------------------------------

test_that("plot_calibration reads the nominal levels off the column names", {
  skip_if_not_installed("ggplot2")
  # A cv_bayes()-shaped result, built by hand: coverage columns per fold and
  # the pooled summary, as cv_bayes() returns them.
  fm <- data.frame(fold = 1:3, n_train = 80, n_test = 40, n_pred = c(40, 38, 42),
                   RMSE = 1, MAE = 1, MAPE = 1, SMAPE = 1, R2 = 0.5, Adj_R2 = NA,
                   n_MAPE = 40L, n_SMAPE = 40L,
                   coverage_50 = c(0.40, 0.45, 0.38), coverage_80 = c(0.70, 0.72, 0.69),
                   coverage_95 = c(0.90, 0.88, 0.91), CRPS = 0.3)
  cv <- list(overall = data.frame(RMSE = 1), fold_metrics = fm,
             predictive_coverage = list(coverage_50 = 0.41, coverage_80 = 0.703,
                                        coverage_95 = 0.897, mean_CRPS = 0.3))
  p <- plot_calibration(cv)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  geoms <- layer_geoms(p)
  pooled <- b$data[[which(geoms == "GeomPoint")[2L]]]
  expect_equal(pooled$x, c(0.5, 0.8, 0.95))
  expect_equal(pooled$y, c(0.41, 0.703, 0.897))
  per_fold <- b$data[[which(geoms == "GeomPoint")[1L]]]
  expect_equal(nrow(per_fold), 9L)
  expect_match(p$labels$subtitle, "too narrow at every level")
  expect_match(p$labels$caption, "3 nominal levels")
  expect_match(p$labels$caption, "seq\\(0.1, 0.9, by = 0.1\\)")

  # Without the pooled summary the fold-weighted mean is used.
  cv$predictive_coverage <- NULL
  b2 <- ggplot2::ggplot_build(plot_calibration(cv))
  pooled2 <- b2$data[[which(geoms == "GeomPoint")[2L]]]
  expect_equal(pooled2$y[1L], sum(fm$coverage_50 * fm$n_pred) / sum(fm$n_pred))
  # A compare_models_cv() result is read through $bayes_cv.
  expect_no_error(plot_calibration(list(overall = data.frame(model = "Bayesian"), bayes_cv = cv)))
  # No coverage columns: refused with the reason.
  expect_error(plot_calibration(list(fold_metrics = fm[, 1:12], overall = data.frame())),
               "no coverage_\\* columns")
  expect_error(plot_calibration(list(a = 1)), "must be the list returned by cv_bayes")
})


# ---- 9.5 -----------------------------------------------------------------

test_that("plot.resolution_profile draws one panel per criterion with the selected level marked", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  pts  <- diag_points(n = 300, seed = 3)
  prof <- suppressWarnings(resolution_profile(pts, response_var = "z", n_levels = 8, seed = 1))
  p <- plot(prof)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  n_panels <- length(unique(b$layout$layout$PANEL))
  drawn <- c("cp", "reliability", "elbow", "moran_z")
  drawn <- drawn[vapply(drawn, function(cn) any(is.finite(prof[[cn]])), logical(1))]
  expect_equal(n_panels, length(drawn))
  # The chosen points agree with select_resolution().
  geoms <- layer_geoms(p)
  chosen <- b$data[[which(geoms == "GeomPoint")[length(which(geoms == "GeomPoint"))]]]
  expect_equal(nrow(chosen), length(drawn))
  best <- vapply(drawn, function(cn) select_resolution(prof, criterion = cn)$best, integer(1))
  # The x axis may be log10-scaled (it is when the ladder spans a factor of
  # 8 or more), in which case the built data carry log10(levels).
  expect_true(setequal(round(chosen$x, 6), round(best, 6)) ||
              setequal(round(chosen$x, 6), round(log10(best), 6)))
  expect_match(p$labels$caption, "Red: the level each criterion selects")

  p2 <- plot(prof, criteria = c("cp", "wss"))
  expect_equal(length(unique(ggplot2::ggplot_build(p2)$layout$layout$PANEL)), 2L)
  expect_error(plot(prof, criteria = "zzz"), "unknown criteria: zzz")
  expect_error(plot.resolution_profile(data.frame(levels = 1)), "must come from resolution_profile")
})

test_that("plot.feature_selection draws the path, the candidates and the stop", {
  skip_if_not_installed("ggplot2")
  pts <- diag_points(n = 200, seed = 7)
  fit_fn <- function(train_sf, vars) lm_spatial_fit(train_sf, "z", vars)
  sel <- suppressMessages(select_features_forward(pts, "z", c("a", "b", "c"),
                                                  fit_fn = fit_fn, k = 3, quiet = TRUE))
  expect_s3_class(sel, "feature_selection")
  expect_true(is.list(sel))
  p <- plot(sel)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  geoms <- layer_geoms(p)
  # Every scored candidate appears once, as a faint or a path point.
  n_pts <- sum(vapply(b$data[which(geoms == "GeomPoint")], nrow, integer(1)))
  n_hist <- sum(is.finite(sel$history$score))
  expect_equal(n_pts, n_hist + length(sel$selected) * 0 + 1L)  # + the chosen marker
  chosen <- b$data[[which(geoms == "GeomPoint")[length(which(geoms == "GeomPoint"))]]]
  expect_equal(chosen$x, length(sel$selected))
  expect_match(p$labels$caption, sprintf("Selected: %s", paste(sel$selected, collapse = ", ")))
  # The lm stand-in scores the null model, so step 0 is on the path.
  expect_true(0L %in% sel$history$step)
  expect_match(p$labels$caption, "Step 0 is the intercept-only model")

  # A hold-out score is drawn as a separate mark and named.
  sel2 <- suppressMessages(suppressWarnings(
    select_features_forward(pts, "z", c("a", "b", "c"), fit_fn = fit_fn, k = 3,
                            quiet = TRUE, select_on = "split")))
  p2 <- plot(sel2)
  if (is.finite(sel2$score_holdout)) {
    expect_match(p2$labels$caption, "hold-out score")
    expect_equal(sum(layer_geoms(p2) == "GeomPoint"), sum(geoms == "GeomPoint") + 1L)
  }
  expect_error(plot.feature_selection(list(selected = "a")), "with its `history` frame")
})

test_that("plot.gwr_model_selection marks the winner and reports its lead", {
  skip_if_not_installed("ggplot2")
  # Hand-built object in the shape gwr_model_selection() returns, so the
  # plot is tested without the GWmodel backend.
  tab <- data.frame(rank = 1:4, n_vars = c(2L, 3L, 1L, 1L),
                    variables = c("a + b", "a + b + c", "a", "b"),
                    criterion = c(500, 510, 700, 800), stringsAsFactors = FALSE)
  x <- structure(list(best = c("a", "b"), table = tab, criterion = "AICc",
                      bandwidth = 60, adaptive = TRUE), class = "gwr_model_selection")
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  geoms <- layer_geoms(p)
  chosen <- b$data[[which(geoms == "GeomPoint")[length(which(geoms == "GeomPoint"))]]]
  expect_equal(chosen$x, 2); expect_equal(chosen$y, 500)
  expect_match(p$labels$subtitle, "4 models evaluated at bandwidth 60 \\(adaptive\\)")
  expect_match(p$labels$subtitle, "leads the runner-up by 10.00")
  expect_match(p$labels$caption, "Criterion: AICc")
  x$table$criterion <- NA_real_
  expect_error(plot(x), "no model has a finite criterion")
})
