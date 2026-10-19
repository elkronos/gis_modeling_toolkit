# tests/testthat/test-fix-integration-callers.R
# ---------------------------------------------------------------------------
# Where one function's change meets another's.
#
# 1. estimate_sac_range() refuses a layer with no spatial structure and
#    attaches a pure nugget.  Its callers ask one predicate for that verdict
#    (.sac_no_structure(), R/utils.R) and act on it the same way.
# 2. predict() on a GWR counts, per call, the locations it cannot estimate
#    and those beyond the reach of the fit.  predict_surface(), which
#    predicts in chunks, and cv_gwr(), which predicts fold by fold, raise
#    each warning once with totals.
# 3. fit_gwr_model() fits the model with no predictors, the kernel-weighted
#    local mean, which select_features_forward() takes as a GWR's baseline.
# 4. model_metrics() refuses the predict() arguments that return something
#    other than one value per row.
# 5. kNNDM folds are named by the functions that label a fold scheme.
# ---------------------------------------------------------------------------

fic_conditions <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = val, conds = w,
       messages = vapply(w, conditionMessage, character(1)),
       classes = vapply(w, function(c) class(c)[1L], character(1)))
}

# 240 points on a 1000 m square: an exponential field (effective range 150,
# nugget 0.2 on a sill of 1) and a white-noise column, in a 4 x 4 grid.
fic_layer <- local({
  val <- NULL
  function() {
    if (is.null(val)) {
      set.seed(4401); n <- 240
      x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
      D <- as.matrix(stats::dist(cbind(x, y)))
      z <- drop(crossprod(chol(0.8 * exp(-3 * D / 150) + diag(0.2, n)), rnorm(n)))
      pts <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, z = z, noise = rnorm(n)),
                          coords = c("x", "y"), crs = 32632)
      bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
      cells <- create_grid_polygons(bnd, target_cells = 16, type = "square")
      val <<- list(pts = pts, cells = cells,
                   asg = suppressWarnings(assign_features_to_polygons(pts, cells)))
    }
    val
  }
})


# --- 1. one verdict for "no spatial structure" --------------------------------

test_that("the two no-structure refusals are the shared constants, and one predicate reads them", {
  skip_if_not_installed("gstat")
  s <- fic_layer()
  no <- suppressWarnings(estimate_sac_range(s$pts, "noise"))
  expect_identical(attr(no, "rejected_reason"), .sac_reason_no_structure)
  expect_true(.sac_no_structure(no))
  expect_true(.vgm_is_pure_nugget(attr(no, "variogram_model")))
  # The flat variogram: white noise let past the test.
  set.seed(70302)
  wn <- sf::st_as_sf(data.frame(x = 5e5 + runif(300, 0, 1000), y = 5e6 + runif(300, 0, 1000),
                                z = rnorm(300)), coords = c("x", "y"), crs = 32632)
  flat <- suppressWarnings(estimate_sac_range(wn, "z", structure_alpha = 1))
  expect_identical(attr(flat, "rejected_reason"), .sac_reason_flat)
  expect_true(.sac_no_structure(flat))
  expect_true(.vgm_is_pure_nugget(attr(flat, "variogram_model")))
  # Not the other refusals, an accepted range, a bare NA or nothing.
  ok  <- estimate_sac_range(s$pts, "z")
  far <- suppressWarnings(estimate_sac_range(s$pts, "z", range_frac = 0.05))
  expect_true(is.finite(ok)); expect_true(is.na(far))
  expect_false(.sac_no_structure(ok))
  expect_false(.sac_no_structure(far))
  expect_false(.sac_no_structure(NA_real_))
  expect_false(.sac_no_structure(NULL))
  expect_false(.sac_no_structure(250))
  # A reason written by hand counts when it says so; a finite range never does.
  by_hand <- structure(NA_real_, rejected_reason = "no spatial structure in this layer")
  expect_true(.sac_no_structure(by_hand))
  expect_false(.sac_no_structure(structure(120, rejected_reason = .sac_reason_no_structure)))
})

test_that("every caller acts on the no-structure verdict, and none treats it as a failed fit", {
  skip_if_not_installed("gstat")
  s <- fic_layer()
  no <- suppressWarnings(estimate_sac_range(s$pts, "noise"))

  # resolution_profile(): no nugget for Cp, no warning about the model, and
  # select_resolution() says the count is not identified.
  rp <- fic_conditions(resolution_profile(s$pts, "noise", sac = no, n_levels = 4, nstart = 2,
                                          estimate = NULL, obs_per_cell = NULL,
                                          partition = "kmeans"))
  expect_identical(attr(rp$value, "cp_noise")$source, "finest-level residual mean square")
  expect_false(any(grepl("not a fitted value|optimiser|did not converge", rp$messages)))
  sel <- suppressWarnings(select_resolution(rp$value, criterion = "cp"))
  expect_false(sel$structured)
  expect_true(any(grepl(.sac_reason_no_structure, sel$note, fixed = TRUE)))

  # summarize_by_cell(): a design effect of 1 by the variogram route.
  sm <- fic_conditions(summarize_by_cell(s$asg, "noise", deff = "variogram", sac = no,
                                         estimate = "mean"))
  expect_length(sm$conds, 0L)
  expect_true(all(sm$value$deff_method == "variogram"))
  expect_true(all(attr(sm$value, "deff_applied")$deff == 1, na.rm = TRUE))
  # The flat variogram takes the same route.
  set.seed(70302)
  wn <- sf::st_as_sf(data.frame(x = 5e5 + runif(300, 0, 1000), y = 5e6 + runif(300, 0, 1000),
                                z = rnorm(300)), coords = c("x", "y"), crs = 32632)
  flat <- suppressWarnings(estimate_sac_range(wn, "z", structure_alpha = 1))
  wn_cells <- create_grid_polygons(sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(wn))),
                                   target_cells = 9, type = "square")
  wn_asg <- suppressWarnings(assign_features_to_polygons(wn, wn_cells))
  sf2 <- fic_conditions(summarize_by_cell(wn_asg, "z", deff = "variogram", sac = flat,
                                          estimate = "mean"))
  expect_false(any(vapply(sf2$conds, inherits, logical(1), "spatialkit_deff_fallback")))
  expect_true(all(sf2$value$deff_method == "variogram"))
  # A refusal that says "no spatial structure" but carries a model with a
  # range is not that verdict's model: it is set aside like any refusal.
  odd <- no
  attr(odd, "variogram_model") <- attr(estimate_sac_range(s$pts, "z"), "variogram_model")
  so <- fic_conditions(summarize_by_cell(s$asg, NULL, predictor_vars = "noise",
                                         deff = "variogram", sac = odd, estimate = "mean"))
  expect_identical(sum(vapply(so$conds, inherits, logical(1), "spatialkit_deff_fallback")), 1L)

  # kriging_adequacy(): the model is used, and the warning says what it is.
  ka <- fic_conditions(kriging_adequacy(s$asg, "noise", s$cells, sac = no, k = 3))
  hit <- grepl("range was not identified", ka$messages, fixed = TRUE)
  expect_identical(sum(hit), 1L)
  expect_match(ka$messages[hit], .sac_reason_no_structure, fixed = TRUE)
  expect_match(ka$messages[hit], "is a pure nugget: no spatial correlation was found")
  expect_false(grepl("did not pin down", ka$messages[hit], fixed = TRUE))
  # Another refusal keeps the sentence written for it.
  far <- suppressWarnings(estimate_sac_range(s$pts, "z", range_frac = 0.05))
  kf <- fic_conditions(kriging_adequacy(s$asg, "z", s$cells, sac = far, k = 3))
  expect_true(any(grepl("is one the data did not pin down", kf$messages, fixed = TRUE)))

  # plot(): the reason is the subtitle.
  skip_if_not_installed("ggplot2")
  expect_match(plot(no)$labels$subtitle, "^No effective range: no spatial structure detected")
})

test_that("a selection's centres carry the profile's hull as their cover", {
  skip_if_not_installed("gstat")
  s <- fic_layer()
  prof <- suppressWarnings(resolution_profile(s$pts, "z", levels = c(4, 9, 16), nstart = 3))
  sel <- suppressWarnings(select_resolution(prof))
  expect_s3_class(attr(sel$seeds, "cover"), "sfc")
  expect_equal(as.numeric(sf::st_area(attr(sel$seeds, "cover"))),
               as.numeric(sf::st_area(attr(prof, "hull"))))
})


# --- 2. the GWR's counted warnings, once per caller ---------------------------

fic_gwr_layer <- function(n = 150, seed = 4) {
  set.seed(seed)
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                               a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + rnorm(n, 0, 0.3)
  d
}

test_that("predict_surface() raises a GWR's extrapolation warning once, with the grid's totals", {
  skip_if_not_installed("sp")
  d <- fic_gwr_layer()
  fit <- suppressMessages(fit_gwr_model(d, "z", "a", bandwidth = 30))
  g <- expand.grid(x = 5e5 + seq(-1500, 2500, length.out = 20),
                   y = 5e6 + seq(-1500, 2500, length.out = 15))
  g$a <- 0
  grid <- sf::st_as_sf(g, coords = c("x", "y"), crs = 32632)
  one  <- fic_conditions(predict_surface(fit, grid, chunk_size = 5000L))
  many <- fic_conditions(predict_surface(fit, grid, chunk_size = 50L))
  for (r in list(one, many)) {
    expect_identical(sum(r$classes == "spatialkit_gwr_extrapolation"), 1L)
    w <- r$conds[[which(r$classes == "spatialkit_gwr_extrapolation")]]
    expect_match(conditionMessage(w), "^predict_surface\\(\\): [0-9]+ of 300 location\\(s\\)")
    expect_identical(w$n, 300L)
  }
  w1 <- one$conds[[which(one$classes == "spatialkit_gwr_extrapolation")]]
  wm <- many$conds[[which(many$classes == "spatialkit_gwr_extrapolation")]]
  expect_identical(wm$n_far, w1$n_far)
  expect_gt(wm$n_far, 50L)
  expect_identical(wm$calls, 6L)
  expect_identical(w1$calls, 1L)
  # One fit, one reach.
  expect_length(wm$reach, 1L)
  expect_equal(wm$reach, w1$reach)
  expect_equal(many$value$.pred, one$value$.pred)
  # A backend that raises neither class is not touched.
  lm_fit <- lm_spatial_fit(d, "z", "a")
  expect_length(fic_conditions(predict_surface(lm_fit, grid, chunk_size = 50L))$conds, 0L)
})

test_that("cv_gwr() raises each counted warning once for all its folds", {
  skip_if_not_installed("sp")
  d <- fic_gwr_layer(n = 160, seed = 9)
  # Four strips as folds, and a fixed bandwidth that does not reach across a
  # strip: every fold has held-out points no training point carries weight
  # for.
  x <- sf::st_coordinates(d)[, 1]
  strips <- cut(x, 4, labels = FALSE)
  per_fold <- fic_conditions(lapply(1:4, function(k) {
    f <- suppressMessages(fit_gwr_model(d[strips != k, ], "z", "a", adaptive = FALSE,
                                        bandwidth = 120))
    predict(f, d[strips == k, ])
  }))
  n_na_folds <- sum(per_fold$classes == "spatialkit_gwr_not_estimable")
  expect_gte(n_na_folds, 2L)                      # one warning per fold that has any
  total_na <- sum(vapply(per_fold$conds[per_fold$classes == "spatialkit_gwr_not_estimable"],
                         function(w) w$n_na, integer(1)))
  cv <- fic_conditions(cv_gwr(d, "z", "a", folds = strips, adaptive = FALSE, bandwidth = 120))
  expect_identical(sum(cv$classes == "spatialkit_gwr_not_estimable"), 1L)
  w <- cv$conds[[which(cv$classes == "spatialkit_gwr_not_estimable")]]
  expect_match(conditionMessage(w), "^cv_gwr\\(\\): ")
  expect_identical(w$n_na, total_na)
  expect_identical(w$calls, n_na_folds)
  expect_identical(sum(is.na(cv$value$predictions$yhat)), total_na)

  # Adaptive folds are separate fits, each with its own reach: the total
  # gives the smallest and the largest.
  far <- d
  sf::st_geometry(far)[strips == 4] <- sf::st_geometry(d)[strips == 4] + c(6000, 0)
  sf::st_crs(far) <- 32632
  strips2 <- cut(sf::st_coordinates(far)[, 1], c(-Inf, 5e5 + 250, 5e5 + 500, 5e5 + 800, Inf),
                 labels = FALSE)
  cv2 <- fic_conditions(cv_gwr(far, "z", "a", folds = strips2, bandwidth = 20))
  expect_identical(sum(cv2$classes == "spatialkit_gwr_extrapolation"), 1L)
  w2 <- cv2$conds[[which(cv2$classes == "spatialkit_gwr_extrapolation")]]
  expect_match(conditionMessage(w2), "^cv_gwr\\(\\): ")
  if (w2$calls > 1L) {
    expect_length(w2$reach, 2L)
    expect_match(conditionMessage(w2), "each call predicted from a fit of its own")
  } else expect_length(w2$reach, 1L)
})


# --- 3. the GWR with no predictors --------------------------------------------

test_that("fit_gwr_model() with no predictors is the kernel-weighted local mean", {
  skip_if_not_installed("sp")
  set.seed(7); n <- 150
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                               a = rnorm(n), n1 = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  xy <- sf::st_coordinates(d)
  d$z <- 1 + sin(xy[, 1] / 300) + 2 * d$a + rnorm(n, 0, 0.3)
  local_mean_at <- function(at, bw = 30, kernel = "bisquare", adaptive = TRUE)
    vapply(seq_len(nrow(at)), function(i) {
      dd <- sqrt((xy[, 1] - at[i, 1])^2 + (xy[, 2] - at[i, 2])^2)
      w  <- .gw_kernel_weights(dd, bw, kernel, adaptive = adaptive)
      sum(w * d$z) / sum(w)
    }, numeric(1))

  r <- fic_conditions(fit_gwr_model(d, "z", character(0), bandwidth = 30))
  f0 <- r$value
  expect_length(r$conds, 0L)
  expect_s3_class(f0, "gwr_fit")
  expect_identical(f0$predictor_vars, character(0))
  expect_identical(deparse(f0$formula), "z ~ 1")
  expect_equal(as.numeric(fitted(f0)), local_mean_at(xy), tolerance = 1e-10)
  expect_identical(names(coef(f0)), "Intercept")
  expect_equal(as.numeric(coef(f0)$Intercept), as.numeric(fitted(f0)), tolerance = 1e-10)
  # At new locations, with standard errors.
  nd <- d[1:20, ]
  sf::st_geometry(nd) <- sf::st_geometry(nd) + c(13, -7)
  sf::st_crs(nd) <- 32632
  expect_equal(as.numeric(predict(f0, nd)), local_mean_at(sf::st_coordinates(nd)),
               tolerance = 1e-10)
  ps <- predict(f0, nd, se.fit = TRUE)
  expect_true(all(is.finite(ps$se.fit)) && all(ps$se.pred > ps$se.fit))
  # A fixed Gaussian kernel, and a bandwidth the fit searches for.
  fg <- suppressMessages(fit_gwr_model(d, "z", character(0), adaptive = FALSE,
                                       bandwidth = 200, kernel = "gaussian"))
  expect_equal(as.numeric(fitted(fg)),
               local_mean_at(xy, bw = 200, kernel = "gaussian", adaptive = FALSE),
               tolerance = 1e-10)
  for (ad in c(TRUE, FALSE)) {
    fs <- fic_conditions(fit_gwr_model(d, "z", character(0), adaptive = ad))
    expect_true(is.finite(fs$value$info$bandwidth))
    expect_true(is.finite(fs$value$info$AICc))
    expect_false(isTRUE(fs$value$info$bandwidth_is_fallback))
  }
  # It prints, is scored, cross-validates and is compared like any GWR.
  expect_output(print(f0), "Formula : z ~ 1", fixed = TRUE)
  expect_true(is.finite(model_metrics(f0)$RMSE))
  cv <- suppressMessages(cv_gwr(d, "z", character(0), k = 3, bandwidth = 30))
  expect_identical(cv$n_folds_succeeded, 3L)
  fa <- suppressMessages(fit_gwr_model(d, "z", "a", bandwidth = 30))
  cmp <- suppressWarnings(compare_models(list(none = f0, a = fa)))
  expect_lt(cmp$AICc[cmp$model == "a"], cmp$AICc[cmp$model == "none"])
  # The map of its one coefficient is drawn, with nothing masked.
  skip_if_not_installed("ggplot2")
  p <- plot(f0, type = "coefficients")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "Intercept")
  # What is not a character vector is still refused by name.
  expect_error(fit_gwr_model(d, "z", NULL, bandwidth = 30),
               "`predictor_vars` must be a character vector")
  expect_error(fit_gwr_model(d, "z", 3, bandwidth = 30),
               "`predictor_vars` must be a character vector")
})

test_that("a GWR sweep takes the local mean as its null model, the baseline `null_fn` used to supply", {
  skip_if_not_installed("sp")
  set.seed(4007); n <- 120
  d <- sf::st_as_sf(data.frame(x = stats::runif(n, 0, 2000), y = stats::runif(n, 0, 2000),
                               a1 = stats::rnorm(n), n1 = stats::rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- d$a1 + stats::rnorm(n, 0, 0.5)
  gwr_fn <- function(tr, v) fit_gwr_model(tr, "z", v, bandwidth = 30)
  lines <- capture_spatialkit_log(
    sel <- suppressWarnings(suppressMessages(
      select_features_forward(d, "z", c("a1", "n1"), gwr_fn, k = 3, seed = 7, quiet = TRUE))))
  expect_identical(sel$baseline$type, "null model")
  expect_null(sel$baseline$first)
  expect_identical(sel$history$variable[1], "<none>")
  expect_true("a1" %in% sel$selected)
  expect_false(log_has(lines, "the baseline is the training-fold mean"))
  # Step 0 is the local mean under the sweep's kernel, on the sweep's folds.
  fo <- suppressWarnings(suppressMessages(
    make_folds(d, k = 3, method = "block_kfold", seed = 7, response_var = "z",
               predictor_vars = c("a1", "n1"))))
  e <- unlist(lapply(fo$folds, function(sp) {
    tr <- d[sp$train, ]; te <- d[sp$test, ]
    txy <- sf::st_coordinates(tr); at <- sf::st_coordinates(te)
    te$z - vapply(seq_len(nrow(at)), function(i) {
      dd <- sqrt((txy[, 1] - at[i, 1])^2 + (txy[, 2] - at[i, 2])^2)
      w  <- .gw_kernel_weights(dd, 30, "bisquare", adaptive = TRUE)
      sum(w * tr$z) / sum(w)
    }, numeric(1))
  }))
  expect_equal(sel$baseline$score, sqrt(mean(e^2)), tolerance = 1e-8)
  # A learner that refuses the empty set still falls back to the fold mean,
  # and the caution says how to get the local mean.
  refusing <- function(tr, v) {
    if (!length(v)) stop("needs a predictor", call. = FALSE)
    fit_gwr_model(tr, "z", v, bandwidth = 30)
  }
  lines2 <- capture_spatialkit_log(
    sel2 <- suppressWarnings(suppressMessages(
      select_features_forward(d, "z", c("a1", "n1"), refusing, k = 3, seed = 7, quiet = TRUE))))
  expect_identical(sel2$baseline$type, "fold mean")
  expect_true(log_has(lines2, "pass character\\(0\\) on to fit_gwr_model\\(\\)"))
})


# --- 4. model_metrics() and predict() arguments that change the shape ---------

test_that("model_metrics() refuses draws = TRUE and se.fit = TRUE by name", {
  d <- surf_test_points(80)
  fit <- lm_spatial_fit(d, "z", "w")
  expect_error(model_metrics(fit, newdata = d, draws = TRUE),
               "model_metrics\\(\\): `draws = TRUE` cannot be passed through `...`.*matrix of posterior draws")
  expect_error(model_metrics(fit, newdata = d, se.fit = TRUE),
               "model_metrics\\(\\): `se.fit = TRUE` cannot be passed.*list with the standard errors")
  # FALSE asks for the default shape and is passed on.
  expect_equal(model_metrics(fit, newdata = d, draws = FALSE)$RMSE,
               model_metrics(fit, newdata = d)$RMSE)
  skip_if_not_installed("sp")
  g <- suppressMessages(fit_gwr_model(fic_gwr_layer(), "z", "a", bandwidth = 30))
  expect_error(model_metrics(g, newdata = fic_gwr_layer(), se.fit = TRUE),
               "`se.fit = TRUE` cannot be passed")
})


# --- 5. kNNDM folds are named --------------------------------------------------

test_that("kNNDM folds get their own label in plot_folds() and print.kriging_adequacy()", {
  s <- fic_layer()
  grid <- sf::st_as_sf(expand.grid(x = 5e5 + seq(50, 950, by = 100),
                                   y = 5e6 + seq(50, 950, by = 100)),
                       coords = c("x", "y"), crs = 32632)
  f <- suppressWarnings(suppressMessages(
    make_folds(s$pts, k = 3, method = "knndm", prediction_points = grid, seed = 1)))
  expect_identical(f$method, "knndm")
  expect_match(.fold_subtitle(f), "kNNDM")
  expect_null(.fold_subtitle(list(method = "something else", params = list())))
  skip_if_not_installed("gstat")
  ok <- estimate_sac_range(s$pts, "z")
  ka <- suppressWarnings(suppressMessages(
    kriging_adequacy(s$asg, "z", s$cells, sac = ok, folds = f)))
  expect_output(print(ka), "kNNDM k-fold CV (knndm, 3 folds", fixed = TRUE)
  skip_if_not_installed("ggplot2")
  expect_match(plot_folds(f, s$pts)$labels$subtitle, "kNNDM")
})


# --- 6. the projection chosen for a lon/lat layer does not depend on row order -

test_that("ensure_projected() and make_folds(auto_range = TRUE) give one answer in any row order", {
  skip_if_not_installed("gstat")
  # 400 points over the conterminous United States, a field of effective
  # range 600 km (the layer of test-fix-variogram-invariance.R).  The
  # candidates for a lon/lat layer were scored on rows taken at even spacing
  # down the layer: as generated it got a Lambert azimuthal CRS, in four of
  # six other row orders an Albers one, and make_folds(), which projects the
  # layer before the range is estimated, sized its blocks at 505,530 or
  # 506,838 accordingly.
  set.seed(14); n <- 400
  lon <- runif(n, -120, -70); lat <- runif(n, 25, 50)
  ll <- sf::st_as_sf(data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326)
  xy <- sf::st_coordinates(sf::st_transform(ll, 5070))
  ll$z <- as.numeric(crossprod(chol(0.8 * exp(-3 * as.matrix(stats::dist(xy)) / 600e3) +
                                      diag(0.2, n)), rnorm(n)))
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  crs0 <- sf::st_crs(quiet(ensure_projected(ll)))
  f0 <- quiet(make_folds(ll, k = 3, method = "block_kfold", auto_range = TRUE,
                         response_var = "z", seed = 1))
  expect_true(is.finite(f0$params$sac_range))
  for (s in 1:6) {
    set.seed(s); p <- sample(n)
    expect_true(sf::st_crs(quiet(ensure_projected(ll[p, ]))) == crs0, info = paste("order", s))
    f <- quiet(make_folds(ll[p, ], k = 3, method = "block_kfold", auto_range = TRUE,
                          response_var = "z", seed = 1))
    expect_identical(as.numeric(f$params$sac_range), as.numeric(f0$params$sac_range),
                     info = paste("order", s))
    expect_identical(f$params$crs, f0$params$crs, info = paste("order", s))
  }
  # The score itself: the same number for the same points in any order.
  set.seed(3); p <- sample(n)
  expect_identical(.crs_distance_error(ll[p, ], sf::st_crs(5070)),
                   .crs_distance_error(ll, sf::st_crs(5070)))
})
