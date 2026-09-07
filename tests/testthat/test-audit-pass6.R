# ===========================================================================
# Regression tests for the sixth audit pass.
#
# One block per finding.  Each one FAILED on the tree these fixes landed on,
# and each asserts the number or the condition that was wrong -- not merely
# that the call returns something.  The pass-6 findings that live naturally
# beside an existing test (the k-NN tie rule, the rotation invariance of
# estimate_sac_range(), the elbow ordering, the mutation-testing gaps) are in
# the topical files; this file holds the rest.
# ===========================================================================

.p6_pts <- function(n = 90, seed = 1, crs = 32632, extent = 1000) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, extent), y = runif(n, 0, extent),
               w = rnorm(n), a = rnorm(n)),
    coords = c("x", "y"), crs = crs
  )
}

# A deterministic lm-backed fit_fn for cv_spatial(): (train_sf) -> spatial_fit.
.p6_lm_fit <- function(train_sf) lm_spatial_fit(train_sf, "z", "w")


# --------------------------------------------------------------------------
# A-01: the inverse-gamma length-scale prior is solved exactly, at any ratio
# --------------------------------------------------------------------------

test_that(".lscale_invgamma pins both tails whatever the ratio or scale of the bounds", {
  # The two tail conditions have one exact solution -- a 1-D root in the
  # shape -- yet the 2-D Nelder-Mead used before stopped in a flat basin from
  # upper/lower ~ 70 upward (every clustered layout) and depended on the
  # absolute scale of the bounds, so the fit silently fell back to the
  # half-normal the same file calls the worse prior.
  ig <- spatialkit:::.lscale_invgamma
  p_below <- function(q, s) 1 - stats::pgamma(1 / q, shape = s$shape, rate = s$scale)
  bounds <- list(c(0.01, 0.65), c(0.01, 0.7), c(0.01, 1), c(0.01, 10),
                 c(0.01305, 1.2013), c(0.00482, 2.0173), c(1, 1e4), c(1e-4, 5))
  for (b in bounds) {
    s <- ig(b[1], b[2])
    expect_true(isTRUE(s$ok), info = paste("bounds", b[1], b[2]))
    expect_equal(p_below(b[1], s), 0.01, tolerance = 1e-6,
                 info = paste("lower tail at", b[1], b[2]))
    expect_equal(1 - p_below(b[2], s), 0.01, tolerance = 1e-6,
                 info = paste("upper tail at", b[1], b[2]))
  }
  # The solution depends on the bounds only through their ratio, up to scale.
  s1 <- ig(0.01, 1); s2 <- ig(1, 100)
  expect_equal(s1$shape, s2$shape, tolerance = 1e-8)
  expect_equal(s2$scale / s1$scale, 100, tolerance = 1e-8)
  # A different tail target is honoured too.
  s <- ig(0.01, 10, tail = 0.05)
  expect_equal(p_below(0.01, s), 0.05, tolerance = 1e-6)
  expect_equal(1 - p_below(10, s), 0.05, tolerance = 1e-6)
})


# --------------------------------------------------------------------------
# A-03: projection candidates are scored against the WGS84 ellipsoid
# --------------------------------------------------------------------------

test_that(".geod_distance is Vincenty's ellipsoidal distance, not a spherical one", {
  gd <- spatialkit:::.geod_distance
  # Vincenty (1975) / Geoscience Australia test line: Flinders Peak to
  # Buninyong, 54 972.271 m.
  dms <- function(d, m, s) sign(d) * (abs(d) + m / 60 + s / 3600)
  expect_equal(gd(dms(144, 25, 29.52440), dms(-37, 57, 3.72030),
                  dms(143, 55, 35.38390), dms(-37, 39, 10.15610)),
               54972.271, tolerance = 1e-8)
  # One degree of longitude on the equator is a * pi / 180 exactly.
  expect_equal(gd(0, 0, 1, 0), 6378137 * pi / 180, tolerance = 1e-9)
  # One degree of latitude from the equator: 110 574.39 m on WGS84 (a sphere
  # of radius 6371 km says 111 194.9 m -- the 0.56% the old measure was off).
  expect_equal(gd(0, 0, 0, 1), 110574.39, tolerance = 1e-6)
  expect_gt(abs(gd(0, 0, 0, 1) / (6371008.8 * pi / 180) - 1), 5e-3)
  # Vectorised, symmetric, zero on identical points.
  lon <- c(13.4, 77.6, -44, 129); lat <- c(52.5, 13, -31, -47)
  d_ab <- gd(lon, lat, rev(lon), rev(lat))
  expect_length(d_ab, 4L)
  expect_equal(d_ab, rev(d_ab))
  expect_equal(gd(lon, lat, lon, lat), rep(0, 4))
})

test_that(".crs_distance_error measures the projection, not the sphere", {
  # sf::st_distance() on lon/lat is an s2 great-circle distance on a sphere
  # of radius 6371 km, which is itself 0.24-0.56% off WGS84 -- the size of
  # the projection errors being ranked.  On a ~1-degree layer the old
  # measure reported UTM's error at ~0.3% where the true worst case is
  # ~0.03%, and in 16 of 40 random wide extents it picked the wrong zone.
  set.seed(3)
  n <- 40
  ll <- sf::st_as_sf(data.frame(lon = 13.4 + runif(n, -0.6, 0.6),
                                lat = 52.5 + runif(n, -0.4, 0.4)),
                     coords = c("lon", "lat"), crs = 4326)
  err <- spatialkit:::.crs_distance_error(ll, sf::st_crs(32633))
  expect_true(is.finite(err))
  expect_lt(err, 1e-3)
  # What the sphere would have said for the same layer, computed here from
  # the great-circle formula: several times larger, and not the projection's.
  xy  <- sf::st_coordinates(ll); rad <- pi / 180
  ij  <- which(lower.tri(matrix(0, n, n)), arr.ind = TRUE)
  la1 <- xy[ij[, 2], 2] * rad; la2 <- xy[ij[, 1], 2] * rad
  dl  <- (xy[ij[, 1], 1] - xy[ij[, 2], 1]) * rad
  d_sph <- 6371008.8 * acos(pmin(1, sin(la1) * sin(la2) + cos(la1) * cos(la2) * cos(dl)))
  d_prj <- as.numeric(stats::dist(sf::st_coordinates(sf::st_transform(ll, 32633))))
  expect_gt(max(abs(d_prj / d_sph - 1)), 2e-3)
  expect_gt(max(abs(d_prj / d_sph - 1)) / err, 3)
})


# --------------------------------------------------------------------------
# B-03: the null-model probe in select_features_forward() is silent
# --------------------------------------------------------------------------

test_that("select_features_forward() does not echo the fold failures of its null probe", {
  # Backends that refuse an empty predictor set make the intercept-only probe
  # fail in every fold; suppressWarnings() silenced the R condition but not
  # the LOGGER lines, so every successful RF/GWR run printed k "fold i fit
  # failed" lines and an "all k folds failed" line -- word for word what a
  # broken run prints -- even with quiet = TRUE.
  d <- .p6_pts(n = 80, seed = 4)
  d$z <- 2 * d$w + rnorm(80, 0, 0.2)
  refusing_fit <- function(train_sf, predictor_vars) {
    if (!length(predictor_vars))
      stop("this backend needs at least one predictor", call. = FALSE)
    lm_spatial_fit(train_sf, "z", predictor_vars)
  }
  lines <- capture_spatialkit_log(
    sel <- suppressWarnings(suppressMessages(
      select_features_forward(d, "z", c("w", "a"), fit_fn = refusing_fit,
                              k = 3, method = "random_kfold", seed = 1,
                              quiet = TRUE))),
    level = logger::INFO)
  expect_true("w" %in% sel$selected)
  expect_false(log_has(lines, "fit failed; skipping"))
  expect_false(log_has(lines, "all 3 folds failed"))
  # The probe's fallback is still announced (as a message, not a log line)
  # when quiet = FALSE.
  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(select_features_forward(d, "z", c("w", "a"),
                                             fit_fn = refusing_fit, k = 3,
                                             method = "random_kfold", seed = 1)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("null \\(intercept-only\\) model could not be scored", msgs)))
})


# --------------------------------------------------------------------------
# C-01 / C-05: conditions name the function the user called
# --------------------------------------------------------------------------

test_that("cv_rf() conditions name cv_rf(), and the fold remapper names no internal", {
  skip_if_not_installed("ranger")
  d <- .p6_pts(n = 80, seed = 5)
  d$z <- 2 * d$w + rnorm(80, 0, 0.2)
  expect_message(
    cv <- cv_rf(d, "z", "w", k = 3, num_trees = 30, seed = 1),
    "^cv_rf\\(\\): no folds supplied")
  expect_equal(cv$n_folds_succeeded, 3L)
  # A fit_fn that is not a function is refused in cv_spatial()'s own name,
  # and in cv_rf()'s when reached through the wrapper's plumbing.
  expect_error(cv_spatial(d, "z", "w", fit_fn = 1), "^cv_spatial\\(\\): `fit_fn`")
  # Folds naming rows that are not in the data: the log line used to open
  # with ".remap_folds():", a function no user ever called.
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  f$folds[[1]]$test <- c(f$folds[[1]]$test, 900L)
  lines <- capture_spatialkit_log(
    suppressMessages(suppressWarnings(
      cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit, folds = f, seed = 1))),
    level = logger::INFO)
  expect_true(log_has(lines, "cross-validation: 1 fold entr"))
  expect_false(log_has(lines, "\\.remap_folds\\(\\)"))
  expect_false(log_has(lines, "\\.cv_run_folds\\(\\)"))
})

test_that("coerce_to_points(mode = 'line_midpoint') names itself when it refuses", {
  ml <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_multilinestring(list(rbind(c(0, 0), c(1, 1)), rbind(c(2, 2), c(3, 3)))),
    crs = 32632))
  expect_error(coerce_to_points(ml, mode = "line_midpoint"),
               "^coerce_to_points\\(\\): method \"line_midpoint\" only supports LINESTRING")
})


# --------------------------------------------------------------------------
# C-02: cv_bayes(seed = ) reaches the sampler, once per fold
# --------------------------------------------------------------------------

test_that("cv_bayes() hands each fold its own sampler seed, drawn from `seed`", {
  # fit_bayesian_spatial_model() carries seed = 123 and cv_bayes() never set
  # it, so every fold of every run sampled from Stan seed 123 and
  # cv_bayes(seed = ) changed nothing on fixed folds.  The backend is mocked:
  # what is under test is the seed cv_bayes() passes, not Stan.
  d <- .p6_pts(n = 60, seed = 6)
  d$z <- 2 * d$w + rnorm(60, 0, 0.2)
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  seen <- list()
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, response_var, predictor_vars,
                                          ..., seed = 123) {
      seen[[length(seen) + 1L]] <<- seed
      lm_spatial_fit(data_sf, response_var, predictor_vars)
    },
    .package = "spatialkit")
  run <- function(seed, fit_args = list()) {
    seen <<- list()
    suppressWarnings(cv_bayes(d, "z", "w", folds = f, seed = seed,
                              fit_args = fit_args,
                              compute_pred_intervals = FALSE))
    unlist(seen)
  }
  s1  <- run(1)
  s1b <- run(1)
  s2  <- run(2)
  expect_length(s1, 3L)
  expect_false(any(s1 == 123L))                 # not the backend default
  expect_equal(length(unique(s1)), 3L)          # one per fold
  expect_identical(s1, s1b)                     # reproducible from `seed`
  expect_false(identical(s1, s2))               # and a function of it
  # A seed given in fit_args is fixed for every fold, as documented.
  expect_identical(run(1, fit_args = list(seed = 7L)), rep(7L, 3L))
})


# --------------------------------------------------------------------------
# C-03: a misspelt `newdata` is an error, not an in-sample answer
# --------------------------------------------------------------------------

test_that("the evaluation functions refuse arguments in `...` when `newdata` is absent", {
  # `...` is forwarded to predict(), which checks it only on the
  # out-of-sample branch.  model_metrics(gwr, newdta = hold) reported the
  # in-sample RMSE of 1.086 where the held-out answer was 25.24, with the
  # same return shape.
  d <- .p6_pts(n = 80, seed = 7)
  d$z <- 2 * d$w + rnorm(80, 0, 0.2)
  fit  <- lm_spatial_fit(d, "z", "w")
  hold <- d[1:20, ]; hold$z <- hold$z + 25
  ok <- model_metrics(fit, newdata = hold)
  expect_gt(ok$RMSE, 20)
  expect_error(model_metrics(fit, newdta = hold),
               "^model_metrics\\(\\): argument\\(s\\) `newdta` were supplied but `newdata` was not")
  expect_error(evaluate_insample(fit, newdta = hold),
               "^evaluate_insample\\(\\): argument\\(s\\) `newdta`")
  expect_error(compare_models(list(lm = fit), newdta = hold),
               "^compare_models\\(\\): argument\\(s\\) `newdta`")
  # The in-sample call with nothing extra still works, and so does the
  # correctly spelt out-of-sample one.
  expect_lt(model_metrics(fit)$RMSE, 1)
  expect_equal(compare_models(list(lm = fit), newdata = hold)$RMSE[1], ok$RMSE)
})


# --------------------------------------------------------------------------
# C-06: one `newdata` contract for every predict() method
# --------------------------------------------------------------------------

test_that("predict() on the three fit classes refuses the same inputs with the same words", {
  # A bare sfc died inside two methods with R's "argument must be coercible
  # to non-negative integer"; a character predictor was refused by name by
  # two and returned all-NA with a generic backend warning from the third; a
  # missing column was reported by two different functions in two wordings.
  chk <- spatialkit:::.check_predict_newdata
  d <- .p6_pts(n = 40, seed = 8)
  d$z <- 2 * d$w + rnorm(40, 0, 0.2)
  fit <- lm_spatial_fit(d, "z", c("w", "a"))
  nd  <- d[1:5, ]
  expect_error(chk(sf::st_geometry(nd), fit, "predict.x_fit"),
               "^predict.x_fit\\(\\): `newdata` is a bare geometry column \\(sfc\\).*'w', 'a'")
  expect_error(chk(sf::st_drop_geometry(nd), fit, "predict.x_fit"),
               "^predict.x_fit\\(\\): `newdata` must be an sf object.*st_as_sf")
  expect_error(chk(nd[, "w"], fit, "predict.x_fit"),
               "^predict.x_fit\\(\\): `newdata` is missing predictor column\\(s\\) 'a'\\.$")
  txt <- nd; txt$a <- as.character(txt$a)
  expect_error(chk(txt, fit, "predict.x_fit"),
               "^predict.x_fit\\(\\): predictor 'a' was numeric at fit time but is character in `newdata`")
  expect_identical(chk(nd, fit, "predict.x_fit"), nd)
  # A logical predictor at fit time is refused if it arrives as character
  # too, and numeric-for-logical is accepted.
  d2 <- d; d2$a <- d2$a > 0
  fit2 <- lm_spatial_fit(d2, "z", c("w", "a"))
  nd2 <- d2[1:5, ]; nd2$a <- as.character(nd2$a)
  expect_error(chk(nd2, fit2, "f"), "predictor 'a' was logical at fit time but is character")
  nd3 <- d2[1:5, ]; nd3$a <- as.numeric(nd3$a)
  expect_identical(chk(nd3, fit2, "f"), nd3)

  # ... and the three methods route through it, so the message is identical
  # apart from the method name.
  same_words <- function(msg, method) {
    expect_match(msg, sprintf("^%s\\(\\): `newdata` is a bare geometry column", method))
  }
  fake_bayes <- structure(
    list(engine = structure(list(), class = "brmsfit"), data_sf = d,
         response_var = "z", predictor_vars = c("w", "a")),
    class = c("bayesian_fit", "spatial_fit"))
  same_words(conditionMessage(tryCatch(predict(fake_bayes, newdata = sf::st_geometry(nd)),
                                       error = identity)), "predict.bayesian_fit")
  if (requireNamespace("ranger", quietly = TRUE)) {
    rf <- fit_rf_model(d, "z", c("w", "a"), num_trees = 20, seed = 1)
    same_words(conditionMessage(tryCatch(predict(rf, newdata = sf::st_geometry(nd)),
                                         error = identity)), "predict.rf_fit")
    expect_error(predict(rf, newdata = txt),
                 "^predict.rf_fit\\(\\): predictor 'a' was numeric at fit time but is character")
    expect_error(predict(rf, newdata = nd[, "w"]),
                 "^predict.rf_fit\\(\\): `newdata` is missing predictor column\\(s\\) 'a'")
  }
  if (requireNamespace("GWmodel", quietly = TRUE) &&
      requireNamespace("sp", quietly = TRUE)) {
    gwr <- suppressWarnings(fit_gwr_model(d, "z", c("w", "a"), bandwidth = 20,
                                          adaptive = TRUE))
    same_words(conditionMessage(tryCatch(predict(gwr, newdata = sf::st_geometry(nd)),
                                         error = identity)), "predict.gwr_fit")
    expect_error(predict(gwr, newdata = txt),
                 "^predict.gwr_fit\\(\\): predictor 'a' was numeric at fit time but is character")
    expect_error(predict(gwr, newdata = nd[, "w"]),
                 "^predict.gwr_fit\\(\\): `newdata` is missing predictor column\\(s\\) 'a'")
    # Unknown arguments in `...` are refused rather than swallowed.
    expect_error(predict(gwr, newdata = nd, newdta = nd),
                 "^predict.gwr_fit\\(\\): unused argument\\(s\\) `newdta`")
  }
})


# --------------------------------------------------------------------------
# C-07: an all-folds-failed run still returns a typed fold_metrics frame
# --------------------------------------------------------------------------

test_that("cv_spatial() returns the metric columns even when every fold failed", {
  d <- .p6_pts(n = 60, seed = 9)
  d$z <- rnorm(60)
  always_fails <- function(train_sf) stop("no fit today", call. = FALSE)
  expect_warning(
    cv <- suppressMessages(cv_spatial(d, "z", "w", fit_fn = always_fails, k = 3, seed = 1)),
    "all folds failed")
  fm <- cv$fold_metrics
  expect_s3_class(fm, "data.frame")
  expect_equal(nrow(fm), 0L)
  expect_true(all(c("fold", "n_train", "n_test", "n_pred", "RMSE", "MAE", "MAPE",
                    "SMAPE", "R2", "Adj_R2") %in% names(fm)))
  # The thing a caller does with it works instead of erroring on a missing
  # column: cv_gwr() and cv_bayes() already returned this shape.
  expect_equal(nrow(subset(fm, RMSE < 5)), 0L)
  expect_equal(cv$n_folds_succeeded, 0L)
  expect_equal(cv$n_folds_attempted, 3L)
})


# --------------------------------------------------------------------------
# C-08: rows that no fold names are reported, not silently dropped
# --------------------------------------------------------------------------

test_that("rows named by no fold raise a warning that counts them", {
  # A folds object built on site[1:45, ] and applied to all 90 rows scored
  # 45 of them and reported attempted = succeeded = 3 with nothing said.
  d <- .p6_pts(n = 90, seed = 10)
  d$z <- 2 * d$w + rnorm(90, 0, 0.2)
  d$..row_id <- seq_len(90)
  half <- make_folds(d[1:45, ], k = 3, method = "random_kfold", seed = 1)
  expect_warning(
    cv <- suppressMessages(cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit,
                                      folds = half, seed = 1)),
    "cross-validation: 45 of 90 rows in the data are named by no fold")
  expect_equal(nrow(cv$predictions), 45L)
  # And nothing is said when every row is covered.
  full <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  expect_no_warning(
    suppressMessages(cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit, folds = full, seed = 1)))
})


# --------------------------------------------------------------------------
# C-10: a vector of fold labels is accepted by every cv_*()
# --------------------------------------------------------------------------

test_that("cv_spatial() accepts a vector of fold labels, one per row", {
  # area_of_applicability() documented and accepted the third shape; the
  # cv_*() functions died on it with "$ operator is invalid for atomic
  # vectors".
  d <- .p6_pts(n = 60, seed = 11)
  d$z <- 2 * d$w + rnorm(60, 0, 0.2)
  d$..row_id <- seq_len(60)
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  labels <- f$assignment$fold[match(d$..row_id, f$assignment$row_id)]
  cv_lab <- suppressMessages(cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit,
                                        folds = labels, seed = 1))
  cv_obj <- suppressMessages(cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit,
                                        folds = f, seed = 1))
  expect_equal(nrow(cv_lab$predictions), 60L)
  # Same train/test partition, so the same held-out predictions per row.
  o1 <- order(cv_lab$predictions$..row_id); o2 <- order(cv_obj$predictions$..row_id)
  expect_equal(cv_lab$predictions$yhat[o1], cv_obj$predictions$yhat[o2])
  expect_equal(cv_lab$overall$RMSE, cv_obj$overall$RMSE)
  # Character labels work the same way; wrong length and a single level are
  # refused by name.
  expect_equal(suppressMessages(cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit,
                                           folds = paste0("f", labels), seed = 1))$overall$RMSE,
               cv_obj$overall$RMSE)
  expect_error(cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit, folds = labels[-1]),
               "^cv_spatial\\(\\): `folds` has 59 labels but the data has 60 rows")
  expect_error(cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit, folds = rep(1L, 60)),
               "^cv_spatial\\(\\): `folds` must define at least two non-empty folds")
  expect_error(cv_spatial(d, "z", "w", fit_fn = .p6_lm_fit, folds = c(NA, labels[-1])),
               "^cv_spatial\\(\\): `folds` contains missing labels")
})


# --------------------------------------------------------------------------
# E-05: block_size is validated, and the grid guard formats its own message
# --------------------------------------------------------------------------

test_that("make_folds(block_kfold) validates block_size and reports a huge grid cleanly", {
  d <- .p6_pts(n = 60, seed = 12)
  for (bad in list(NA, c(100, 200), -250, 0, "250", Inf)) {
    expect_error(make_folds(d, k = 4, method = "block_kfold", block_size = bad, seed = 1),
                 "^make_folds\\(\\): `block_size` must be a single positive number",
                 info = paste("block_size =", paste(bad, collapse = ",")))
  }
  # A unit mistake asks for a grid past 2^31 cells on a side; %d used to die
  # with "invalid format" instead of the documented refusal.
  expect_error(make_folds(d, k = 4, method = "block_kfold", block_size = 1e-9, seed = 1),
               "the requested grid is [0-9]{12} x [0-9]{12} = .* cells, above the 1,000,000")
  expect_error(make_folds(d, k = 4, method = "block_kfold", block_size = 1e-6, seed = 1),
               "above the 1,000,000 this function will build")
})


# --------------------------------------------------------------------------
# E-06: estimate_sac_range() says which variable it modelled
# --------------------------------------------------------------------------

test_that("estimate_sac_range() refuses an unknown predictor and reports a failed detrending", {
  skip_if_not_installed("gstat")
  d <- .p6_pts(n = 150, seed = 13)
  xy <- sf::st_coordinates(d)
  D  <- as.matrix(stats::dist(xy))
  d$trend <- 0.003 * xy[, 1]
  d$z <- d$trend + as.numeric(t(chol(exp(-D / 60) + diag(1e-6, 150))) %*% rnorm(150))
  # intersect() used to drop the unknown names silently and model the RAW
  # response, so make_folds(auto_range = TRUE) sized blocks from the wrong
  # variogram with nothing said.
  expect_error(estimate_sac_range(d, "z", predictor_vars = c("trnd", "a")),
               "^estimate_sac_range\\(\\): predictor_vars 'trnd' not found in the data")
  expect_error(make_folds(d, k = 4, method = "block_kfold", auto_range = TRUE,
                          response_var = "z", predictor_vars = "trnd", seed = 1),
               "predictor_vars 'trnd' not found")
  ok <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "trend"))
  expect_true(attr(ok, "detrended"))
  raw <- suppressWarnings(estimate_sac_range(d, "z"))
  expect_false(attr(raw, "detrended"))
  # An all-NA predictor makes lm() fail: the fallback to the raw response is
  # a different estimand, so it is an R warning, and the attribute says so.
  d$bad <- NA_real_
  expect_warning(
    fb <- estimate_sac_range(d, "z", predictor_vars = "bad"),
    "OLS detrending on bad failed .* fitted to the RAW response")
  expect_false(attr(fb, "detrended"))
  expect_equal(as.numeric(fb), as.numeric(raw))
})


# --------------------------------------------------------------------------
# E-07: a boundary that misses the points is refused, a partial one reported
# --------------------------------------------------------------------------

test_that("make_folds(block_kfold, boundary = ) says when the boundary does not cover the points", {
  d <- .p6_pts(n = 80, seed = 14)
  sq <- function(x0, y0, side) sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + side, y0), c(x0 + side, y0 + side), c(x0, y0 + side), c(x0, y0)))),
    crs = 32632))
  far <- sq(5000, 5000, 1000)
  expect_error(make_folds(d, k = 4, method = "block_kfold", boundary = far, seed = 1),
               "^make_folds\\(block_kfold\\): none of the 80 points fall inside `boundary`")
  # A boundary covering part of the extent: the points outside are counted
  # and the region is extended, and every point still gets a block.
  part <- sq(0, 0, 600)
  n_out <- sum(!lengths(sf::st_intersects(d, part)))
  expect_gt(n_out, 0L)
  expect_warning(
    f <- make_folds(d, k = 4, method = "block_kfold", boundary = part, seed = 1),
    sprintf("%d of 80 points fall outside `boundary`; the block region has been extended", n_out))
  expect_equal(sort(unlist(lapply(f$folds, `[[`, "test"))), seq_len(80))
  # The single-block refusal names what produced the grid: an automatic grid
  # over a zero-extent layer, not "block_nx/block_ny" the caller never passed.
  one_place <- sf::st_as_sf(data.frame(x = rep(500, 20), y = rep(500, 20), w = rnorm(20)),
                            coords = c("x", "y"), crs = 32632)
  expect_error(suppressWarnings(make_folds(one_place, k = 4, method = "block_kfold", seed = 1)),
               "the automatic grid \\(block_multiplier x k blocks over the extent\\) produces a single block")
  expect_error(suppressWarnings(make_folds(one_place, k = 4, method = "block_kfold",
                                           block_nx = 1, block_ny = 1, seed = 1)),
               "the requested block_nx/block_ny produces a single block")
})


# --------------------------------------------------------------------------
# F-01: the package's log formatter is pinned, not inherited
# --------------------------------------------------------------------------

test_that("spatialkit's log lines are immune to the user's global logger formatter", {
  # Every helper hands logger an ALREADY-formatted string.  The namespace
  # inherited the user's global formatter, so under formatter_sprintf every
  # message with a literal "%" -- the CRS distortion figures, the GWR
  # collinearity percentage -- hard-errored with "too few arguments", and
  # because .warn_and_log() logs before it warns, the R warning died with
  # it.  Under the default formatter_glue a "{...}" in a fold error was
  # re-evaluated.
  # logger returns the formatter as the expression that generated it.
  old <- eval(logger::log_formatter(namespace = "global"), asNamespace("logger"))
  on.exit(logger::log_formatter(old, namespace = "global"), add = TRUE)
  logger::log_formatter(logger::formatter_sprintf, namespace = "global")
  # Re-run the package's own initialisation with the hostile global config
  # in place, as if the user had configured logger BEFORE library(spatialkit).
  spatialkit:::.onLoad(NULL, "spatialkit")
  on.exit(spatialkit:::.onLoad(NULL, "spatialkit"), add = TRUE)

  lines <- capture_spatialkit_log({
    expect_no_error(spatialkit:::.log_warn("distance error %.2f%% against %.2f%%", 0.49, 2.37))
    expect_warning(spatialkit:::.warn_and_log("solver diverged at {iter=3}; %d%% done", 50),
                   "solver diverged at \\{iter=3\\}; 50% done")
  })
  expect_true(log_has(lines, "distance error 0.49% against 2.37%"))
  expect_true(log_has(lines, "solver diverged at \\{iter=3\\}; 50% done"))
  expect_false(log_has(lines, "diverged at 3"))
})


# --------------------------------------------------------------------------
# F-02 / H-01: ranger runs on one thread unless the session opts in
# --------------------------------------------------------------------------

test_that("fit_rf_model() and predict.rf_fit() default to one thread, or mc.cores", {
  skip_if_not_installed("ranger")
  # NULL reached ranger as "use every core" (hardware_concurrency()), and
  # cv_rf(parallel = TRUE) multiplied that by the worker count.  CRAN's
  # limit is two.  ranger keeps no record of the thread count, so the call
  # is intercepted and what it received is written to a file -- which also
  # survives the fork that cv_rf(parallel = 2) makes.
  seen <- tempfile(fileext = ".txt")
  real_ranger <- ranger::ranger
  local_mocked_bindings(
    ranger = function(..., num.threads = NULL) {
      cat(sprintf("fit %s\n", format(num.threads)), file = seen, append = TRUE)
      real_ranger(..., num.threads = num.threads)
    },
    .package = "ranger")
  # predict() dispatches on the engine's class through the S3 registry, so a
  # probe subclass in front of "ranger" sees what predict.rf_fit() passes.
  registerS3method("predict", "p6_ranger_probe",
                   function(object, data, ..., num.threads = NULL) {
                     cat(sprintf("predict %s\n", format(num.threads)), file = seen, append = TRUE)
                     NextMethod()
                   })
  read_seen <- function() {
    x <- if (file.exists(seen)) readLines(seen) else character(0)
    unlink(seen)
    x
  }
  d <- .p6_pts(n = 80, seed = 15)
  d$z <- 2 * d$w + rnorm(80, 0, 0.2)

  fit <- fit_rf_model(d, "z", "w", num_trees = 20, seed = 1)
  expect_identical(read_seen(), "fit 1")
  class(fit$engine) <- c("p6_ranger_probe", class(fit$engine))
  invisible(predict(fit, newdata = d[1:10, ]))
  expect_identical(read_seen(), "predict 1")
  # The session's mc.cores opt-in is honoured by both ...
  withr_opt <- options(mc.cores = 2L); on.exit(options(withr_opt), add = TRUE)
  fit2 <- fit_rf_model(d, "z", "w", num_trees = 20, seed = 1)
  class(fit2$engine) <- c("p6_ranger_probe", class(fit2$engine))
  invisible(predict(fit2, newdata = d[1:10, ]))
  expect_identical(read_seen(), c("fit 2", "predict 2"))
  # ... an explicit value wins in both places, and the forest does not
  # depend on the count.
  fit3 <- fit_rf_model(d, "z", "w", num_trees = 20, seed = 1, num_threads = 1)
  expect_identical(read_seen(), "fit 1")
  class(fit3$engine) <- c("p6_ranger_probe", class(fit3$engine))
  invisible(predict(fit3, newdata = d[1:10, ], num.threads = 1))
  expect_identical(read_seen(), "predict 1")
  expect_equal(fitted(fit3), fitted(fit2))
  options(mc.cores = NULL)
})

test_that("cv_rf() runs each forked fold on one thread", {
  skip_if_not_installed("ranger")
  skip_on_cran()
  skip_on_os("windows")
  seen <- tempfile(fileext = ".txt")
  real_ranger <- ranger::ranger
  local_mocked_bindings(
    ranger = function(..., num.threads = NULL) {
      cat(sprintf("%s\n", format(num.threads)), file = seen, append = TRUE)
      real_ranger(..., num.threads = num.threads)
    },
    .package = "ranger")
  d <- .p6_pts(n = 80, seed = 16)
  d$z <- 2 * d$w + rnorm(80, 0, 0.2)
  withr_opt <- options(mc.cores = 2L); on.exit(options(withr_opt), add = TRUE)
  n_workers <- spatialkit:::.resolve_n_cores(2L)
  suppressMessages(cv_rf(d, "z", "w", k = 3, num_trees = 20, seed = 1, parallel = 2L))
  got <- readLines(seen); unlink(seen)
  expect_length(got, 3L)
  # With forks, one thread each; sequentially the session opt-in applies.
  expect_true(all(got == if (n_workers > 1L) "1" else "2"))
  # An explicit num_threads still passes through under forks.
  suppressMessages(cv_rf(d, "z", "w", k = 3, num_trees = 20, seed = 1, parallel = 2L,
                         num_threads = 2L))
  got <- readLines(seen); unlink(seen)
  expect_true(all(got == "2"))
  options(mc.cores = NULL)
})


# --------------------------------------------------------------------------
# F-04: seed = NULL still makes the parallel path reproducible from set.seed()
# --------------------------------------------------------------------------

test_that("cv_spatial(seed = NULL) is reproducible from set.seed(), in parallel too", {
  d <- make_cv_test_points(n = 120, seed = 1)
  run <- function(parallel) {
    set.seed(777)
    suppressMessages(cv_spatial(d, "z", "w", fit_fn = boot_fit_fn, k = 4,
                                seed = NULL, parallel = parallel))$overall$RMSE
  }
  expect_equal(run(FALSE), run(FALSE))
  # The caller's stream is consumed, as any RNG-using call would consume it:
  # a different state gives a different answer, and the state moves on.
  set.seed(778)
  other <- suppressMessages(cv_spatial(d, "z", "w", fit_fn = boot_fit_fn, k = 4,
                                       seed = NULL))$overall$RMSE
  expect_false(isTRUE(all.equal(other, run(FALSE))))
  set.seed(777); invisible(suppressMessages(cv_spatial(d, "z", "w", fit_fn = boot_fit_fn,
                                                       k = 4, seed = NULL)))
  after <- runif(1)
  set.seed(777); fresh <- runif(1)
  expect_false(after == fresh)
})

test_that("cv_spatial(seed = NULL, parallel = 2) reproduces the sequential run", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(spatialkit:::.resolve_n_cores(2L) < 2L, "single-core machine")
  # Unseeded forked workers were seeded by mclapply() from the clock and the
  # process ID: set.seed(777); cv_*(seed = NULL, parallel = 2) gave 0.5687 /
  # 0.5649 / 0.5623 on three runs while the sequential call was reproducible.
  d <- make_cv_test_points(n = 120, seed = 1)
  run <- function(parallel) {
    set.seed(777)
    suppressMessages(cv_spatial(d, "z", "w", fit_fn = boot_fit_fn, k = 4,
                                seed = NULL, parallel = parallel))$overall$RMSE
  }
  expect_equal(run(2L), run(2L))
  expect_equal(run(2L), run(FALSE))
})


# --------------------------------------------------------------------------
# F-05: clear_grid_cache() removes only its own entries
# --------------------------------------------------------------------------

test_that("clear_grid_cache() leaves the user's other objects in the environment alone", {
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))), crs = 32632))
  e <- new.env(parent = emptyenv())
  assign("important", 42, envir = e)
  assign("my_results", list(a = 1), envir = e)
  create_grid_polygons_cached(bnd, target_cells = 4, cache_env = e)
  create_grid_polygons_cached(bnd, target_cells = 9, cache_env = e)
  keys <- ls(e, all.names = TRUE)
  expect_equal(sum(startsWith(keys, "spatialkit_grid::")), 2L)
  # rm(ls()) used to wipe every binding and count the user's objects among
  # the "entries removed".
  expect_equal(clear_grid_cache(e), 2L)
  expect_identical(sort(ls(e, all.names = TRUE)), c("important", "my_results"))
  expect_equal(get("important", envir = e), 42)
  expect_equal(clear_grid_cache(e), 0L)
})


# --------------------------------------------------------------------------
# F-06: warnings raised inside a forked fold reach the caller
# --------------------------------------------------------------------------

test_that("warnings raised inside a fold are relayed from the parallel path", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(spatialkit:::.resolve_n_cores(2L) < 2L, "single-core machine")
  # R conditions do not cross the fork, so fit_gwr_model()'s documented
  # integer-response warning, raised in every fold, reached nobody under
  # parallel = 2 while the sequential run showed all of them.
  d <- .p6_pts(n = 80, seed = 17)
  d$z <- 2 * d$w + rnorm(80, 0, 0.2)
  warning_fit <- function(train_sf) {
    warning("this fold has something to say", call. = FALSE)
    lm_spatial_fit(train_sf, "z", "w")
  }
  w <- character(0)
  res <- withCallingHandlers(
    suppressMessages(cv_spatial(d, "z", "w", fit_fn = warning_fit, k = 4, seed = 1,
                                parallel = 2L)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_equal(res$n_folds_succeeded, 4L)
  expect_true(any(grepl("this fold has something to say", w)))
  # Once per distinct message, not once per fold.
  expect_equal(sum(grepl("this fold has something to say", w)), 1L)
  # And the numbers are the sequential ones.
  seq_res <- suppressWarnings(suppressMessages(
    cv_spatial(d, "z", "w", fit_fn = warning_fit, k = 4, seed = 1)))
  expect_equal(res$overall$RMSE, seq_res$overall$RMSE)
})


# --------------------------------------------------------------------------
# G-04: the response's design effect comes from the RESPONSE variogram
# --------------------------------------------------------------------------

test_that("summarize_by_cell(deff = 'variogram') corrects the response SE with the response variogram", {
  skip_if_not_installed("gstat")
  # With predictor_vars supplied the internal variogram was a RESIDUAL one --
  # the correlation of the part the predictors do not explain, which is
  # weaker -- and grand-mean coverage fell from 0.93 to 0.51 the moment a
  # predictor was listed, which is the normal thing to do.
  set.seed(18); n <- 240
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  D <- as.matrix(stats::dist(cbind(x, y)))
  p <- as.numeric(t(chol(exp(-D / 100) + diag(1e-6, n))) %*% rnorm(n))   # smooth predictor
  z <- 2 * p + as.numeric(t(chol(exp(-D / 15) + diag(1e-6, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, p = p,
                                 poly_id = 1L + (x > 500) + 2L * (y > 500)),
                      coords = c("x", "y"), crs = 32632)
  resp_only <- suppressWarnings(summarize_by_cell(pts, "z", deff = "variogram"))
  with_pred <- suppressWarnings(summarize_by_cell(pts, "z", predictor_vars = "p",
                                                  deff = "variogram"))
  expect_equal(with_pred[["..se_resp_z"]], resp_only[["..se_resp_z"]], tolerance = 1e-10)
  expect_equal(attr(with_pred, "deff_applied")$deff, attr(resp_only, "deff_applied")$deff,
               tolerance = 1e-10)
  # ... which is what passing the response variogram explicitly gives, and
  # NOT what the residual variogram gives.
  sac_resp  <- suppressWarnings(estimate_sac_range(pts, "z"))
  sac_resid <- suppressWarnings(estimate_sac_range(pts, "z", predictor_vars = "p"))
  skip_if(is.na(sac_resp) || is.na(sac_resid), "the variograms did not fit on this draw")
  explicit <- summarize_by_cell(pts, "z", predictor_vars = "p", deff = "variogram",
                                sac = sac_resp)
  expect_equal(with_pred[["..se_resp_z"]], explicit[["..se_resp_z"]], tolerance = 1e-10)
  residual <- summarize_by_cell(pts, "z", predictor_vars = "p", deff = "variogram",
                                sac = sac_resid)
  expect_false(isTRUE(all.equal(with_pred[["..se_resp_z"]], residual[["..se_resp_z"]],
                                tolerance = 1e-6)))
})


# --------------------------------------------------------------------------
# H-02: core counts default to the session's mc.cores opt-in, capped
# --------------------------------------------------------------------------

test_that("parallel defaults follow getOption('mc.cores') and never exceed the machine", {
  # `cores = max(1L, parallel::detectCores() - 1L)` was the documented
  # default for the Bayesian backend and the auto-detect path of every
  # cv_*(): 63 workers on a 64-core host, and a hard error under
  # _R_CHECK_LIMIT_CORES_.
  expect_identical(formals(fit_bayesian_spatial_model)$cores,
                   quote(getOption("mc.cores", 1L)))
  withr_opt <- options(mc.cores = 1L); on.exit(options(withr_opt), add = TRUE)
  expect_identical(spatialkit:::.resolve_n_cores(TRUE), 1L)
  options(mc.cores = NULL)
  auto <- suppressMessages(spatialkit:::.resolve_n_cores(TRUE))
  expect_lte(auto, max(1L, parallel::detectCores(logical = TRUE)))
  detected <- parallel::detectCores(logical = FALSE)
  if (!is.na(detected) && detected > 1L) expect_lte(auto, detected - 1L)
})
