# tests/testthat/test-review3-final-A.R
# ---------------------------------------------------------------------------
# Third review, final pass: folds for area_of_applicability() after
# prep_model_data() removed rows, and their provenance; a boundary without a
# CRS in make_folds(), cv_*() and predict_surface(); a `sac` that is a bare
# range in resolution_profile() and summarize_by_cell(); and cross-validation
# cautions raised with .warn_and_log() instead of a log line and a warning.
# ---------------------------------------------------------------------------

# 120 points on a 1 km square with a response, two predictors and a spatial
# trend, in a projected CRS.
.fa_pts <- function(n = 120, seed = 3) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  a <- rnorm(n); b <- rnorm(n)
  sf::st_as_sf(data.frame(x = x, y = y, a = a, b = b,
                          z = a + 0.002 * x + rnorm(n, 0, 0.3)),
               coords = c("x", "y"), crs = 32632)
}

# Conditions an expression raises, muffled: warning messages and classes, and
# the number of messages.
.fa_conditions <- function(expr) {
  w_msg <- character(0); w_cls <- list(); n_msg <- 0L
  val <- withCallingHandlers(expr,
    warning = function(w) {
      w_msg <<- c(w_msg, conditionMessage(w))
      w_cls <<- c(w_cls, list(class(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      n_msg <<- n_msg + 1L
      invokeRestart("muffleMessage")
    })
  list(value = val, warnings = w_msg, classes = w_cls, messages = n_msg)
}

.fa_quiet <- function(expr) {
  logger::with_log_threshold(expr, threshold = logger::FATAL,
                             namespace = "spatialkit", index = 2)
}


# --- area_of_applicability(): folds after prep_model_data() dropped rows ----

test_that("area_of_applicability() takes the folds cv_*() takes when prep_model_data() dropped a row", {
  pts <- .fa_pts()
  pts$z[50] <- NA
  folds <- make_folds(pts, k = 4, method = "block_kfold", seed = 1)
  prep  <- .fa_quiet(prep_model_data(pts, "z", c("a", "b")))
  expect_identical(nrow(prep), 119L)
  fit <- lm_spatial_fit(prep, "z", c("a", "b"))
  new <- pts[1:5, ]

  # It stopped with "fold 1 refers to rows outside 1:119" -- for the folds
  # cv_*() accepts on the same layer.
  lines <- capture_spatialkit_log(
    res <- area_of_applicability(new, model = fit, folds = folds))
  expect_true(log_has(lines, "name 1 row\\(s\\) that prep_model_data\\(\\) removed"))

  # The same threshold as the folds with that row taken out by hand, as
  # positions in the model's training data.
  kept <- setdiff(seq_len(nrow(pts)), 50L)
  by_hand <- lapply(folds$folds, function(s)
    list(train = match(setdiff(s$train, 50L), kept),
         test  = match(setdiff(s$test, 50L), kept)))
  ref <- area_of_applicability(new, model = fit, folds = by_hand)
  expect_equal(res$threshold, ref$threshold)
  expect_equal(res$train_DI, ref$train_DI)

  # A label vector with one label per row of the layer fitted from.
  lab <- area_of_applicability(new, model = fit, folds = folds$assignment$fold)
  expect_equal(lab$threshold, ref$threshold)
  # The bare splits, which name row 120, past the 119 training rows.
  bare <- area_of_applicability(new, model = fit, folds = folds$folds)
  expect_equal(bare$threshold, ref$threshold)

  # Folds built on the model's own training data are positions in it, and
  # are still read that way: the same threshold as on the layer without its
  # record, which `[` removes.
  f_own <- make_folds(prep, k = 4, method = "block_kfold", seed = 1)
  own   <- area_of_applicability(new, model = fit, folds = f_own)
  plain <- prep[seq_len(nrow(prep)), ]
  expect_null(attr(plain, "dropped"))
  expect_equal(own$threshold,
               area_of_applicability(new, train_sf = plain,
                                     predictor_vars = c("a", "b"),
                                     folds = f_own)$threshold)

  # The layer's own ..row_id is honoured the same way.
  p2 <- pts
  p2$..row_id <- 1000L + seq_len(nrow(p2))
  f2 <- make_folds(p2, k = 4, method = "block_kfold", seed = 1)
  fit2 <- lm_spatial_fit(.fa_quiet(prep_model_data(p2, "z", c("a", "b"))),
                         "z", c("a", "b"))
  expect_equal(area_of_applicability(new, model = fit2, folds = f2)$threshold,
               ref$threshold)
})

test_that("a fold ID naming a row the data never had is still refused", {
  pts <- .fa_pts()
  pts$z[50] <- NA
  folds <- make_folds(pts, k = 4, method = "block_kfold", seed = 1)
  fit <- lm_spatial_fit(.fa_quiet(prep_model_data(pts, "z", c("a", "b"))),
                        "z", c("a", "b"))
  folds$folds[[1]]$test <- c(folds$folds[[1]]$test, 999L)
  expect_error(area_of_applicability(pts[1:5, ], model = fit, folds = folds),
               "names 1 ..row_id value\\(s\\) that are not in the training data \\(first: 999\\)")
})

test_that("area_of_applicability() refuses folds built on other rows, as cv_*() do", {
  pts  <- .fa_pts()
  prep <- .fa_quiet(prep_model_data(pts, "z", c("a", "b")))
  fb   <- make_folds(prep, k = 4, method = "block_kfold", seed = 1)
  ok   <- area_of_applicability(prep[1:5, ],
                                model = lm_spatial_fit(prep, "z", c("a", "b")),
                                folds = fb)
  expect_true(is.finite(ok$threshold))
  # The same rows in another order: every fold ID now names another point.
  sorted <- prep[order(prep$z), ]
  expect_error(
    area_of_applicability(prep[1:5, ],
                          model = lm_spatial_fit(sorted, "z", c("a", "b")),
                          folds = fb),
    "built from different data")
  # A bare list of splits carries no probe and is read as positions, as before.
  expect_no_error(
    area_of_applicability(prep[1:5, ],
                          model = lm_spatial_fit(sorted, "z", c("a", "b")),
                          folds = fb$folds))
})

test_that("folds built on polygons still apply to a fit on the points they were reduced to", {
  pts <- .fa_pts()
  polys <- sf::st_sf(sf::st_drop_geometry(pts),
                     geometry = sf::st_geometry(sf::st_buffer(pts, 10)))
  f_poly <- make_folds(polys, k = 4, method = "block_kfold", seed = 1)
  prep   <- .fa_quiet(prep_model_data(polys, "z", c("a", "b")))
  expect_true(all(sf::st_geometry_type(prep) == "POINT"))
  fit <- lm_spatial_fit(prep, "z", c("a", "b"))
  lines <- capture_spatialkit_log(
    res <- area_of_applicability(prep[1:5, ], model = fit, folds = f_poly))
  expect_true(is.finite(res$threshold))
  expect_true(log_has(lines, "skipping the provenance check"))
})


# --- a boundary without a CRS ------------------------------------------------

test_that("make_folds() warns about a CRS-less boundary, prediction_points and blocks, naming them", {
  pts <- .fa_pts()
  bnd <- sf::st_set_crs(sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(pts))), NA)
  # Named by the function and the argument, and the CRS stamped by its label:
  # "the supplied `crs`" named an argument make_folds() does not have.
  expect_warning(
    f <- .fa_quiet(make_folds(pts, k = 4, method = "block_kfold", seed = 1,
                              boundary = bnd)),
    "make_folds\\(\\): `boundary` has no CRS.*stamping the target CRS \\('EPSG:32632'\\) WITHOUT reprojection")
  ref <- make_folds(pts, k = 4, method = "block_kfold", seed = 1,
                    boundary = sf::st_set_crs(bnd, 32632))
  expect_identical(f$assignment$fold, ref$assignment$fold)

  blk <- sf::st_set_crs(sf::st_make_grid(bnd, n = c(3, 3)), NA)
  expect_warning(
    .fa_quiet(make_folds(pts, k = 4, method = "block_kfold", seed = 1,
                         blocks = sf::st_sf(geometry = blk))),
    "make_folds\\(\\): `blocks` has no CRS")

  grid <- sf::st_set_crs(sf::st_as_sf(sf::st_sample(sf::st_as_sfc(sf::st_bbox(pts)),
                                                    64, type = "regular")), NA)
  expect_warning(
    .fa_quiet(make_folds(pts[1:60, ], method = "nndm", prediction_points = grid)),
    "make_folds\\(\\): `prediction_points` has no CRS")
})

test_that("cv_*() warn once about a CRS-less lon/lat boundary, naming the caller", {
  pts <- .fa_pts()
  bll <- sf::st_set_crs(sf::st_transform(
    sf::st_as_sf(sf::st_buffer(sf::st_as_sfc(sf::st_bbox(pts)), 50)), 4326), NA)
  # prep_model_data() and make_folds() each warned about it, and neither
  # named cv_spatial() or `boundary`.
  res <- .fa_conditions(.fa_quiet(
    cv_spatial(pts, "z", "a", fit_fn = function(tr) lm_spatial_fit(tr, "z", "a"),
               k = 3, boundary = bll)))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "^cv_spatial\\(\\): `boundary` has no CRS; its coordinates look like lon/lat")
  expect_identical(res$value$overall$n_pred, nrow(pts))

  # A planar one is stamped with the data's CRS, with one warning naming the
  # caller rather than make_folds().
  bnd <- sf::st_set_crs(sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(pts))), NA)
  res <- .fa_conditions(.fa_quiet(
    cv_spatial(pts, "z", "a", fit_fn = function(tr) lm_spatial_fit(tr, "z", "a"),
               k = 3, boundary = bnd)))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "^cv_spatial\\(\\): `boundary` has no CRS and its coordinates do not look like lon/lat")
})

test_that("predict_surface() warns about a CRS-less boundary, naming it", {
  pts <- surf_test_points(n = 60)
  fit <- lm_spatial_fit(pts, "z", "w")
  bnd <- sf::st_set_crs(sf::st_as_sfc(sf::st_bbox(pts)), NA)
  expect_warning(
    g <- .fa_quiet(predict_surface(fit, n_cells = 100, covariates = pts,
                                   boundary = bnd)),
    "predict_surface\\(\\): `boundary` has no CRS")
  expect_gt(nrow(g), 0L)
})


# --- `sac` as a bare range ---------------------------------------------------

test_that("resolution_profile() refuses a units `sac` and a character one by name", {
  skip_if_not_installed("units")
  pts <- .fa_pts()
  # set_units(1.5, "km") was read as a range of 1.5 m: a floor of millions of
  # cells, with nothing said.
  expect_error(resolution_profile(pts, "z", sac = units::set_units(1.5, "km"),
                                  n_levels = 4),
               "`sac` must be an estimate_sac_range\\(\\) result.*got 1\\.5 \\[km\\]")
  expect_error(resolution_profile(pts, "z", sac = "300", n_levels = 4),
               "got an object of class character")
})

test_that("summarize_by_cell() says so when a `sac` without a model is set aside", {
  pts <- .fa_pts()
  pts$poly_id <- rep(1:6, each = 20)
  # With no response to estimate a variogram from, the fallback warning names
  # the value that could not be used.
  res <- .fa_conditions(.fa_quiet(
    summarize_by_cell(pts, deff = "variogram", sac = 300)))
  expect_length(res$warnings, 1L)
  expect_true("spatialkit_deff_fallback" %in% res$classes[[1L]])
  expect_match(res$warnings, "supplied `sac` \\(300\\) carries no fitted variogram model")

  skip_if_not_installed("gstat")
  # With one, the value was ignored and a variogram estimated without a word.
  res <- .fa_conditions(.fa_quiet(
    summarize_by_cell(pts, "z", deff = "variogram", sac = 300)))
  hit <- grepl("supplied `sac` \\(300\\) carries no fitted variogram model", res$warnings)
  expect_identical(sum(hit), 1L)
  expect_match(res$warnings[hit], "It was set aside, and the design effect uses the variogram estimated from `response_var` instead")
  expect_false("spatialkit_deff_fallback" %in% res$classes[[which(hit)]])
})


# --- cautions raised once ----------------------------------------------------

test_that("make_folds() cautions are logged and raised with the same text, once under knitr", {
  old <- spatialkit_quiet(FALSE)
  withr::defer(spatialkit_quiet(old))
  pts <- .fa_pts()
  lines <- capture_spatialkit_log(
    expect_warning(make_folds(pts, k = 4, method = "block_kfold", seed = 1,
                              auto_range = TRUE),
                   "make_folds(): auto_range requires response_var; ignoring.",
                   fixed = TRUE))
  expect_true(log_has(lines, "auto_range requires response_var; ignoring"))

  withr::local_options(knitr.in.progress = TRUE)
  res <- .fa_conditions(utils::capture.output(
    f <- make_folds(pts, k = 4, method = "block_kfold", seed = 1,
                    auto_range = TRUE),
    type = "message"))
  expect_identical(res$messages, 0L)
  expect_length(res$warnings, 1L)
})

test_that("an all-failed cross-validation logs the text it warns with", {
  pts <- .fa_pts()
  fails <- function(tr) stop("no fit here")
  lines <- capture_spatialkit_log(
    expect_warning(suppressMessages(cv_spatial(pts, "z", "a", fit_fn = fails, k = 3)),
                   paste0("cv_spatial(): all folds failed (all 3 folds failed to ",
                          "produce predictions); cross-validation results ",
                          "contain no predictions. First error:"),
                   fixed = TRUE))
  expect_true(log_has(lines, paste0("cv_spatial\\(\\): all folds failed \\(all 3 folds ",
                                    "failed to produce predictions\\); cross-validation ",
                                    "results contain no predictions\\. First error:")))
})
