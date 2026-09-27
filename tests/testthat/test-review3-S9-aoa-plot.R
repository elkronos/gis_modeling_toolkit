# tests/testthat/test-review3-S9-aoa-plot.R
# ---------------------------------------------------------------------------
# Third review pass over the area of applicability, the prediction grid and
# the plots: rows outside on a dropped predictor that the AOA plot left out,
# a grid that did not cover the training extent, fold labels that crashed
# area_of_applicability(), a weights hint that could not work, and captions
# and subtitles that said something the data did not.
# ---------------------------------------------------------------------------

r3_layer_geoms <- function(p) unname(vapply(p$layers, function(l) class(l$geom)[1L], character(1)))

# 60 training rows with a zero-variance land-cover dummy, and 40 prediction
# rows of which the last 10 take a value the training data never has on it.
r3_aoa_dropped <- function(lc_new = rep(0:1, c(30, 10)), na_first = TRUE) {
  set.seed(1)
  train <- sf::st_as_sf(data.frame(x = runif(60), y = runif(60), a = rnorm(60),
                                   b = rnorm(60), lc = 0),
                        coords = c("x", "y"), crs = 32632)
  new <- sf::st_as_sf(data.frame(x = runif(40), y = runif(40), a = rnorm(40),
                                 b = rnorm(40), lc = lc_new),
                      coords = c("x", "y"), crs = 32632)
  if (na_first) new$a[1] <- NA
  suppressWarnings(area_of_applicability(new, train_sf = train,
                                         predictor_vars = c("a", "b", "lc")))
}


# ---- plot.aoa(): DI = Inf and DI = NA rows ---------------------------------

test_that("plot.aoa counts DI = Inf rows in the prediction curve and names them", {
  skip_if_not_installed("ggplot2")
  res <- r3_aoa_dropped()
  n_inf <- sum(is.infinite(res$aoa$DI))
  expect_identical(n_inf, 10L)
  expect_identical(res$n_na, 1L)
  p <- plot(res)
  b <- ggplot2::ggplot_build(p)
  # The step layer stays first and grouped by set.
  expect_identical(r3_layer_geoms(p)[1L], "GeomStep")
  expect_equal(length(unique(b$data[[1]]$group)), 2L)
  pred <- b$data[[1]][b$data[[1]]$colour == "#2166AC", ]
  # Its height at the threshold is the share inside among the rows with a DI;
  # the Inf rows used to be dropped, and the curve read 0.97 inside where
  # the subtitle counted 11 of 40 outside.
  at_thr <- max(pred$y[pred$x <= res$threshold])
  expect_equal(at_thr, res$n_inside / (res$n_new - res$n_na))
  # It tops out at the finite share, not at 1.
  n_fin <- res$n_new - res$n_na - n_inf
  expect_equal(max(pred$y), n_fin / (n_fin + n_inf))
  # The training curve still reaches 1.
  tr <- b$data[[1]][b$data[[1]]$colour == "grey45", ]
  expect_equal(max(tr$y), 1)
  cap <- p$labels$caption
  expect_match(cap, "10 prediction locations outside on a dropped predictor \\(DI = Inf\\) are off the axis")
  expect_match(cap, "1 with a missing predictor \\(DI = NA\\) is not drawn")
  expect_no_error(ggplot2::ggplot_build(plot(res, type = "histogram")))
  expect_match(plot(res, type = "histogram")$labels$caption, "DI = Inf")
})

test_that("plot.aoa gives the true reason when every row is outside on a dropped predictor", {
  skip_if_not_installed("ggplot2")
  res <- r3_aoa_dropped(lc_new = rep(1, 40), na_first = FALSE)
  expect_true(all(is.infinite(res$aoa$DI)))
  expect_error(plot(res), "40 of 40 prediction locations are outside on a predictor dropped")
  expect_error(plot(res), "DI = Inf")
  expect_false(tryCatch(plot(res), error = function(e)
    grepl("missing or non-finite predictor", conditionMessage(e))))
})

test_that("plot.aoa draws no DI = Inf or NA caption lines when there are none", {
  skip_if_not_installed("ggplot2")
  res <- r3_aoa_dropped(lc_new = rep(0, 40), na_first = FALSE)
  p <- plot(res)
  expect_false(grepl("DI = Inf|DI = NA", p$labels$caption))
  pred <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(max(pred$y[pred$colour == "#2166AC"]), 1)
})


# ---- plot.aoa(): the fold line of the caption --------------------------------

test_that("plot.aoa's caption reads correctly for every kind of fold input", {
  skip_if_not_installed("ggplot2")
  set.seed(2); n <- 120
  train <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                                   a = rnorm(n), b = rnorm(n)),
                        coords = c("x", "y"), crs = 32632)
  new <- train[1:20, ]
  aoa_with <- function(folds) area_of_applicability(new, train_sf = train,
                                                    predictor_vars = c("a", "b"),
                                                    folds = folds)
  # No folds: the threshold is small and the AOA conservative, as
  # ?area_of_applicability says; the caption called it optimistic.
  cap0 <- plot(aoa_with(NULL))$labels$caption
  expect_match(cap0, "not cross-validated")
  expect_match(cap0, "conservative")
  expect_false(grepl("optimistic", cap0))
  # Labels and bare splits read "over from fold labels; method unknown folds".
  cap_l <- plot(aoa_with(rep(1:4, 30)))$labels$caption
  expect_match(cap_l, "cross-validated over folds given as labels \\(method unknown\\)")
  f <- make_folds(train, k = 4, method = "block_kfold", block_size = 300, seed = 1)
  cap_s <- plot(aoa_with(f$folds))$labels$caption
  expect_match(cap_s, "cross-validated over folds given as train/test splits \\(method unknown\\)")
  for (cap in c(cap_l, cap_s)) expect_false(grepl("over from|unknown folds", cap))
  expect_match(plot(aoa_with(f))$labels$caption, "cross-validated over block_kfold folds")
})


# ---- area_of_applicability(): fold labels and weights ------------------------

test_that("numeric fold labels that print alike are grouped as cv_*() groups them", {
  set.seed(3); n <- 40
  tr <- sf::st_as_sf(data.frame(x = runif(n), y = runif(n), a = rnorm(n), b = rnorm(n)),
                     coords = c("x", "y"), crs = 32632)
  # 0.3 and 0.1 + 0.2 differ in the last bits and print alike: one fold in
  # cv_*(); area_of_applicability() stopped with "factor level [2] is
  # duplicated".
  lab <- rep(c(0.3, 0.1 + 0.2, 0.5, 0.7), each = 10)
  expect_false(identical(0.3, 0.1 + 0.2))
  expect_no_error(res <- area_of_applicability(tr[1:5, ], train_sf = tr,
                                               predictor_vars = c("a", "b"),
                                               folds = lab))
  expect_identical(res$params$n_folds, 3L)
  tests_of <- function(sp) lapply(sp, `[[`, "test")
  aoa_sp <- spatialkit:::.aoa_fold_splits(lab, n)
  cv_sp  <- spatialkit:::.folds_from_labels(lab, data.frame(..row_id = seq_len(n)),
                                            "cv_spatial")
  expect_identical(tests_of(aoa_sp), tests_of(cv_sp))
  expect_identical(tests_of(aoa_sp), list(1:20, 21:30, 31:40))
  # The refusals are unchanged.
  expect_error(spatialkit:::.aoa_fold_splits(c(1, NA, 2, 2), 4),
               "area_of_applicability\\(\\): `folds` contains missing labels")
  expect_error(spatialkit:::.aoa_fold_splits(rep(1, 4), 4),
               "area_of_applicability\\(\\): `folds` must define at least two")
  expect_error(spatialkit:::.aoa_fold_splits(1:3, 4), "has 3 labels but the training data has 4 rows")
})

test_that("a NaN weight is refused with its cause, not with the pmax() advice", {
  set.seed(4); n <- 40
  tr <- sf::st_as_sf(data.frame(x = runif(n), y = runif(n), a = rnorm(n), b = rnorm(n)),
                     coords = c("x", "y"), crs = 32632)
  aoa_w <- function(w) area_of_applicability(tr[1:5, ], train_sf = tr,
                                             predictor_vars = c("a", "b"), weights = w)
  # pmax(NaN, 0) is NaN, so the advice reproduced the error it came with.
  expect_error(aoa_w(c(a = NaN, b = 1)), "`weights` is NaN for .a.\\. .*no row is out of bag")
  expect_error(aoa_w(c(a = NaN, b = 1)), "weights = NULL")
  err <- tryCatch(aoa_w(c(a = NaN, b = 1)), error = conditionMessage)
  expect_false(grepl("so pass pmax\\(importance, 0\\)", err))
  # Inf and NA are named without the out-of-bag story.
  err_inf <- tryCatch(aoa_w(c(a = 1, b = Inf)), error = conditionMessage)
  expect_match(err_inf, "`weights` is not finite for .b.\\.$")
  expect_match(tryCatch(aoa_w(c(a = NA, b = 1)), error = conditionMessage),
               "`weights` is not finite for .a.")
  # A negative weight keeps the advice that works for it.
  expect_error(aoa_w(c(a = -0.1, b = 1)), "pmax\\(importance, 0\\)")
})

test_that("the NaN importance of a forest with no out-of-bag rows gets the NaN message", {
  skip_if_not_installed("ranger")
  set.seed(1); n <- 60
  tr <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                                a = rnorm(n), b = rnorm(n)),
                     coords = c("x", "y"), crs = 32632)
  tr$z <- tr$a + rnorm(n)
  fit <- suppressWarnings(suppressMessages(
    fit_rf_model(tr, "z", c("a", "b"), num_trees = 20, seed = 1,
                 replace = FALSE, sample_fraction = 1)))
  w <- pmax(fit$info$importance, 0)
  expect_true(all(is.nan(w)))
  expect_error(area_of_applicability(tr, model = fit, weights = w),
               "is NaN for .a., .b.\\. Permutation importance is NaN when no row is out of bag")
})


# ---- predict_surface(): the grid covers the training extent ------------------

test_that("the automatic grid covers every training point and is centred on the box", {
  pts <- surf_test_points(n = 120)
  fit <- lm_spatial_fit(pts)
  bb  <- sf::st_bbox(pts)
  P   <- sf::st_coordinates(pts)
  for (cs in c(100, 334, 77.7)) {
    g  <- predict_surface(fit, cell_size = cs, covariates = pts)
    xy <- sf::st_coordinates(g)
    ux <- sort(unique(xy[, 1])); uy <- sort(unique(xy[, 2]))
    # floor() cells from xmin left up to a cell uncovered on the east and
    # north: 14 of 120 training points in no cell at cell_size = 100.
    inside <- P[, 1] >= min(ux) - cs / 2 & P[, 1] <= max(ux) + cs / 2 &
              P[, 2] >= min(uy) - cs / 2 & P[, 2] <= max(uy) + cs / 2
    expect_true(all(inside), label = sprintf("cell_size %s", cs))
    # Symmetric about the box, every centre inside it, spacing exact.
    expect_equal(min(ux) - bb[["xmin"]], bb[["xmax"]] - max(ux), tolerance = 1e-8)
    expect_equal(min(uy) - bb[["ymin"]], bb[["ymax"]] - max(uy), tolerance = 1e-8)
    expect_gt(min(ux), bb[["xmin"]]); expect_lt(max(ux), bb[["xmax"]])
    expect_true(all(abs(diff(ux) - cs) < 1e-6))
    # The overhang is less than one cell.
    expect_lt(length(ux) * cs - (bb[["xmax"]] - bb[["xmin"]]), cs)
  }
  # An exact multiple gains no column: 1000 / 100 is 10, 0.3 / 0.1 is 3.
  grid_fn <- spatialkit:::.make_prediction_grid
  crs <- sf::st_crs(32632)
  bb1 <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 500), crs = crs)
  g1 <- grid_fn(bb1, crs, cell_size = 100)
  expect_equal(sort(unique(g1$..grid_x)), seq(50, 950, by = 100))
  expect_equal(sort(unique(g1$..grid_y)), seq(50, 450, by = 100))
  bb2 <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 0.3, ymax = 0.3), crs = crs)
  expect_identical(length(unique(grid_fn(bb2, crs, cell_size = 0.1)$..grid_x)), 3L)
})


# ---- plot_folds(): the unit in the subtitle ---------------------------------

test_that("plot_folds' subtitle gives the CRS's linear unit, not its identifier", {
  skip_if_not_installed("ggplot2")
  set.seed(1); n <- 80
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  p4 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs"
  pts_p4  <- sf::st_as_sf(xy, coords = c("x", "y"), crs = p4)
  pts_wkt <- sf::st_set_crs(sf::st_set_crs(pts_p4, NA), sf::st_crs(sf::st_crs(pts_p4)$wkt))
  pts_epsg <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  sub_of <- function(pts, ...) {
    f <- suppressWarnings(make_folds(pts, seed = 1, ...))
    plot_folds(f, pts)$labels$subtitle
  }
  # The proj string, the whole ~1300-character WKT, or "EPSG:32632 units".
  for (pts in list(pts_p4, pts_wkt, pts_epsg)) {
    sub <- sub_of(pts, k = 5, method = "block_kfold", block_size = 300)
    expect_match(sub, "^Block size 300 \\(metre\\)\n")
    expect_lt(nchar(sub), 120)
  }
  expect_identical(sub_of(pts_epsg, method = "buffered_loo", buffer = 100),
                   "Leave-one-out with a 100 buffer (metre)")
})


# ---- plot_tessellation_map(): labels ----------------------------------------

test_that("plot_tessellation_map labels lon/lat cells quietly and draws a units label", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  pts <- sf::st_as_sf(data.frame(x = 5e5 + runif(40, 0, 1000), y = 5e6 + runif(40, 0, 1000)),
                      coords = c("x", "y"), crs = 32632)
  cells <- build_tessellation(pts, method = "voronoi", quiet = TRUE)$cells
  # geom_sf_text() ran st_point_on_surface() again on the label points at
  # print, which warned for every lon/lat layer.
  p_ll <- plot_tessellation_map(sf::st_transform(cells, 4326), labels = TRUE)
  expect_no_warning(b_ll <- ggplot2::ggplot_build(p_ll))
  txt <- b_ll$data[[which(r3_layer_geoms(p_ll) == "GeomText")]]
  expect_equal(nrow(txt), nrow(cells))
  # An st_area() label failed at print: "units package is not attached".
  cells$area <- sf::st_area(cells)
  p_u <- plot_tessellation_map(cells, labels = TRUE, label_col = "area")
  expect_no_error(b_u <- ggplot2::ggplot_build(p_u))
  lab <- b_u$data[[which(r3_layer_geoms(p_u) == "GeomText")]]$label
  expect_true(all(grepl("\\[m\\^2\\]$", lab)))
  cells$lag <- as.difftime(seq_len(nrow(cells)), units = "days")
  p_d <- plot_tessellation_map(cells, labels = TRUE, label_col = "lag")
  expect_no_error(ggplot2::ggplot_build(p_d))
  # A vector of names is refused by name, not with R's coercion error.
  expect_error(plot_tessellation_map(cells, fill_col = c("area", "cell_id")),
               "`fill_col` must be a single column name")
  expect_error(plot_tessellation_map(cells, labels = TRUE, label_col = c("area", "cell_id")),
               "`label_col` must be a single column name")
  expect_error(plot_tessellation_map(cells, fill_col = NA_character_),
               "`fill_col` must be a single column name")
})


# ---- plot_cv_metrics(): counts and unnamed models ------------------------------

test_that("plot_cv_metrics draws no pooled line for a count column", {
  skip_if_not_installed("ggplot2")
  # overall$n_pred is the total over the folds: the line sat at 150 against
  # folds of 30.
  cv <- list(
    fold_metrics = data.frame(fold = 1:5, n_pred = c(30, 30, 32, 28, 30),
                              RMSE = c(1, 1.2, 0.9, 1.1, 1), n_MAPE = c(30, 30, 32, 28, 30)),
    overall = data.frame(RMSE = 1.04, n_pred = 150L, n_MAPE = 150L))
  for (m in c("n_pred", "n_MAPE")) {
    p <- plot_cv_metrics(cv, m)
    expect_false("GeomHline" %in% r3_layer_geoms(p))
    expect_match(p$labels$caption, sprintf("No pooled value: `%s` is a count per fold; `overall` holds the total", m))
    expect_no_error(ggplot2::ggplot_build(p))
  }
  # RMSE is unchanged.
  b <- ggplot2::ggplot_build(p <- plot_cv_metrics(cv, "RMSE"))
  expect_equal(b$data[[which(r3_layer_geoms(p) == "GeomHline")]]$yintercept, 1.04)
})

test_that("plot_cv_metrics names a model with no per-fold values and keeps the survivor's strip", {
  skip_if_not_installed("ggplot2")
  # RF's bandwidth is NA in every fold and in `overall`; it vanished
  # unmentioned, and with one model left no strip said which one was drawn.
  cmp <- list(
    by_fold = data.frame(model = rep(c("GWR", "RF"), each = 3), fold = rep(1:3, 2),
                         bandwidth = c(40, 42, 41, NA, NA, NA), n_pred = 25,
                         stringsAsFactors = FALSE),
    overall = data.frame(model = c("GWR", "RF"), RMSE = c(1.8, 1.7),
                         stringsAsFactors = FALSE))
  p <- plot_cv_metrics(cmp, "bandwidth")
  expect_match(p$labels$caption, "No pooled value")
  expect_match(p$labels$caption, "Not drawn: RF \\(no finite per-fold `bandwidth`\\)")
  expect_s3_class(p$facet, "FacetWrap")
  b <- ggplot2::ggplot_build(p)
  expect_identical(as.character(b$layout$layout$model), "GWR")
  # A single cv_*() result still has no strip.
  cv <- list(fold_metrics = data.frame(fold = 1:3, RMSE = c(1, 2, 3)),
             overall = data.frame(RMSE = 2))
  expect_false(inherits(plot_cv_metrics(cv, "RMSE")$facet, "FacetWrap"))
})


# ---- plot.spatial_fit(type = "variogram"): the no-range clause ----------------

test_that("the overlay caption says no range was identified, and why for the response", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  set.seed(5); n <- 120
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 80) + diag(1e-4, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 3857)
  sac <- estimate_sac_range(pts, "z")
  skip_if_not(is.finite(sac))
  no_range <- function(s, reason) {
    s[1L] <- NA_real_
    attr(s, "rejected_reason") <- reason
    s
  }
  draw <- function(main, ov) .draw_sac_variogram(main, what = "Residual variogram",
                                                  overlay = ov, overlay_label = "Response (z)")
  cap <- function(p) gsub("\n", " ", p$labels$caption)
  # A non-converged response fit was captioned as "reached no identified
  # sill", as if the curve were still rising.
  c1 <- cap(draw(sac, no_range(sac, "variogram model did not converge")))
  expect_match(c1, "Sills not compared: the response variogram has no identified range \\(variogram model did not converge\\)")
  expect_false(grepl("reached no", c1))
  c2 <- cap(draw(no_range(sac, NULL), sac))
  expect_match(c2, "Sills not compared: the residual variogram has no identified range\\.")
  c3 <- cap(draw(no_range(sac, NULL),
                 no_range(sac, "empirical variogram decreases with distance")))
  expect_match(c3, "neither variogram has an identified range \\(response: empirical variogram decreases with distance\\)")
})
