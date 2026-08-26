# =============================================================================
# Random forest backend
#
# .rf_frame() runs without ranger; everything else is gated.  The frame builder
# is separated out precisely because fitting and prediction must produce
# identical columns in identical order -- ranger matches by name, and a silent
# mismatch there would produce plausible-looking wrong predictions rather than
# an error.
# =============================================================================

mk_rf_pts <- function(n = 200, seed = 21, crs = 32632) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
    coords = c("x", "y"), crs = crs
  )
  d$z <- 2 * d$a - d$b + rnorm(n, 0, 0.3)
  d
}


# --- .rf_frame ---------------------------------------------------------------

test_that(".rf_frame selects the predictors in order and drops geometry", {
  pts <- mk_rf_pts(10)
  f <- .rf_frame(pts, c("b", "a"))
  expect_s3_class(f, "data.frame")
  expect_identical(names(f), c("b", "a"))
  expect_identical(nrow(f), 10L)
  expect_false("geometry" %in% names(f))
  expect_false("z" %in% names(f))
})

test_that(".rf_frame appends coordinates only when asked", {
  pts <- mk_rf_pts(10)
  expect_identical(names(.rf_frame(pts, "a", include_coords = FALSE)), "a")
  f <- .rf_frame(pts, "a", include_coords = TRUE)
  expect_identical(names(f), c("a", "..x", "..y"))
  expect_equal(unname(f[["..x"]]), unname(sf::st_coordinates(pts)[, 1]))
  expect_equal(unname(f[["..y"]]), unname(sf::st_coordinates(pts)[, 2]))
})

test_that(".rf_frame produces identical columns for training and new data", {
  # ranger matches columns by name; a mismatch here would silently mispredict.
  tr <- mk_rf_pts(10)
  nd <- mk_rf_pts(4, seed = 99)
  nd$z <- NULL                                   # true out-of-sample
  # Pass the predictors in a DIFFERENT order for newdata: .rf_frame() must
  # impose the order it is given, so the caller controls alignment.  Handing
  # both calls the same vector could not have detected a reordering.
  for (ic in c(FALSE, TRUE)) {
    a <- .rf_frame(tr, c("a", "b"), include_coords = ic)
    b <- .rf_frame(nd, c("a", "b"), include_coords = ic)
    expect_identical(names(a), names(b))
    rev_b <- .rf_frame(nd, c("b", "a"), include_coords = ic)
    expect_identical(names(rev_b)[1:2], c("b", "a"))
    expect_false(identical(names(a), names(rev_b)))
  }
})

test_that(".rf_frame names the missing predictor", {
  pts <- mk_rf_pts(5)
  expect_error(.rf_frame(pts, c("a", "nope")), "absent from the data.*nope")
})


# --- fitting -----------------------------------------------------------------

test_that("fit_rf_model returns a spatial_fit with populated info", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts()
  fit <- fit_rf_model(dat, "z", c("a", "b", "noise"), num_trees = 200L)
  expect_s3_class(fit, "rf_fit")
  expect_s3_class(fit, "spatial_fit")
  expect_identical(fit$n, 200L)
  expect_identical(fit$info$num_trees, 200L)
  expect_false(fit$info$include_coords)
  expect_true(fit$info$fitted_are_oob)
  expect_true(is.finite(fit$info$oob_rmse))
  expect_true(is.finite(fit$info$oob_r_squared))
  expect_setequal(names(fit$info$importance), c("a", "b", "noise"))
})

test_that("permutation importance ranks the real predictors above noise", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts()
  fit <- fit_rf_model(dat, "z", c("a", "b", "noise"), num_trees = 300L)
  imp <- fit$info$importance
  expect_gt(imp[["a"]], imp[["noise"]])
  expect_gt(imp[["b"]], imp[["noise"]])
  expect_gt(imp[["a"]], imp[["b"]])              # coefficient 2 versus -1
})

test_that("fitted() returns out-of-bag predictions, not in-sample ones", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts()
  fit <- fit_rf_model(dat, "z", c("a", "b"), num_trees = 300L)
  y <- dat$z

  oob_rmse <- sqrt(mean((y - fitted(fit))^2))
  ins_rmse <- sqrt(mean((y - predict(fit, newdata = dat))^2))

  # In-sample predictions from a forest are close to memorisation.  If fitted()
  # returned those, summary() would report a fictitious R-squared.
  expect_lt(ins_rmse, oob_rmse)
  expect_equal(oob_rmse, fit$info$oob_rmse, tolerance = 1e-8)
  expect_equal(residuals(fit), y - fitted(fit))
  expect_length(fitted(fit), fit$n)
})

test_that("summary() on an rf_fit reports the out-of-bag numbers", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts()
  fit <- fit_rf_model(dat, "z", c("a", "b"), num_trees = 200L)
  s <- summary(fit)
  expect_equal(s$in_sample$RMSE, fit$info$oob_rmse, tolerance = 1e-8)
})

test_that("fit_rf_model is reproducible from its seed", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts(120)
  nd  <- mk_rf_pts(30, seed = 77)
  f1 <- fit_rf_model(dat, "z", c("a", "b"), num_trees = 200L, seed = 42L)
  f2 <- fit_rf_model(dat, "z", c("a", "b"), num_trees = 200L, seed = 42L)
  expect_identical(predict(f1, newdata = nd), predict(f2, newdata = nd))
  expect_identical(fitted(f1), fitted(f2))
  f3 <- fit_rf_model(dat, "z", c("a", "b"), num_trees = 200L, seed = 7L)
  expect_false(identical(fitted(f1), fitted(f3)))
})

test_that("include_coords logs a warning and shows the coordinates in the formula", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts(120)
  # .log_warn() writes through the logger; it raises no R condition, so
  # expect_warning() would pass vacuously here.  Capture the log instead.
  logged <- capture_spatialkit_log(
    fit <- fit_rf_model(dat, "z", "a", num_trees = 100L, include_coords = TRUE))
  expect_match(paste(logged, collapse = "\n"), "include_coords = TRUE")
  expect_match(paste(logged, collapse = "\n"), "memorising location")
  expect_true(fit$info$include_coords)
  expect_true(all(c("..x", "..y") %in% all.vars(fit$formula)))
  expect_setequal(names(fit$info$importance), c("..x", "..y", "a"))
})

test_that("fit_rf_model rejects a non-numeric response", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts(40)
  dat$grp <- letters[1:2]
  expect_error(fit_rf_model(dat, "grp", "a"), "not numeric")
})


# --- prediction --------------------------------------------------------------

test_that("predict.rf_fit aligns NA rows to their original positions", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts(120)
  fit <- fit_rf_model(dat, "z", c("a", "b"), num_trees = 100L)
  nd  <- mk_rf_pts(10, seed = 55)
  nd$a[3] <- NA
  p <- predict(fit, newdata = nd)
  expect_length(p, 10L)
  expect_true(is.na(p[3]))
  expect_false(anyNA(p[-3]))
})

test_that("predict.rf_fit reprojects newdata to the fitting CRS", {
  skip_if_not_installed("ranger")
  # Only discriminating when the coordinates are predictors: otherwise the CRS
  # cannot change the answer.  Degrees fed to a forest trained on metres would
  # land every point in the same corner of predictor space.
  dat <- mk_rf_pts(150)
  fit <- fit_rf_model(dat, "z", "a", num_trees = 200L, include_coords = TRUE)

  nd_proj <- mk_rf_pts(20, seed = 66)
  nd_ll   <- sf::st_transform(nd_proj, 4326)
  expect_true(sf::st_is_longlat(nd_ll))

  # Loose tolerance on purpose: a forest is a step function, so a nanometre of
  # round-trip reprojection error can flip a single tree across a split.  What
  # this catches is degrees being fed to a forest trained on metres, which
  # would collapse every point into one corner of predictor space.
  expect_equal(predict(fit, newdata = nd_ll), predict(fit, newdata = nd_proj),
               tolerance = 1e-3)

  # The tolerance above is wide enough that it depends on the installed PROJ
  # version, which varies across the CI matrix.  This is the discriminating
  # half: mislabel the lon/lat coordinates as metres -- exactly what skipping
  # the reprojection does -- and the predictions must NOT come back equal.
  # Without it, a `tolerance` that quietly grew too large would go unnoticed.
  nd_mislabelled <- sf::st_set_crs(
    sf::st_set_crs(nd_ll, NA_character_), sf::st_crs(nd_proj))
  p_wrong <- suppressWarnings(predict(fit, newdata = nd_mislabelled))
  p_right <- predict(fit, newdata = nd_proj)
  expect_length(p_wrong, length(p_right))
  expect_false(isTRUE(all.equal(unname(p_wrong), unname(p_right),
                                tolerance = 1e-3)))
  # ... and wrong by a lot, not by rounding: degrees put every point within
  # ~50 units of the origin, far outside the training extent.
  expect_gt(max(abs(p_wrong - p_right)), 1e-2 * stats::sd(p_right))
})

test_that("coef() on an rf_fit refuses rather than inventing something", {
  skip_if_not_installed("ranger")
  fit <- fit_rf_model(mk_rf_pts(60), "z", c("a", "b"), num_trees = 100L)
  expect_error(coef(fit), "no coefficients")
  expect_error(coef(fit), "importance")
})

test_that("print.rf_fit states the out-of-bag caveat", {
  skip_if_not_installed("ranger")
  fit <- fit_rf_model(mk_rf_pts(80), "z", c("a", "b"), num_trees = 100L)
  txt <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(txt, "Random Forest \\(ranger\\)")
  expect_match(txt, "OOB RMSE")
  expect_match(txt, "optimistic under spatial")
  expect_match(txt, "cv_rf")
  expect_match(txt, "Coords as predictors: no")
  invisible(utils::capture.output(vis <- withVisible(print(fit))))
  expect_false(vis$visible)
})


# --- integration with the rest of the package --------------------------------

test_that("cv_rf produces a spatially blocked estimate", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts(200)
  cv  <- suppressMessages(cv_rf(dat, "z", c("a", "b"), k = 4,
                                num_trees = 100L, seed = 1))
  expect_true(is.finite(cv$overall$RMSE))
  expect_true(is.finite(cv$overall$R2))
  expect_gte(nrow(cv$fold_metrics), 2L)
  expect_lte(nrow(cv$fold_metrics), 4L)
})

test_that("predict_surface works from an rf_fit", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts(150)
  fit <- fit_rf_model(dat, "z", "a", num_trees = 100L, include_coords = TRUE)

  # A bare grid carries geometry but no covariates, so the model's predictors
  # have to come from somewhere; predict_surface() says so rather than
  # silently predicting from nothing.
  expect_error(predict_surface(fit, n_cells = 200L), "absent from the grid")

  grd <- predict_surface(fit, n_cells = 200L, covariates = dat)
  expect_true(inherits(grd, "sf"))
  expect_true(".pred" %in% names(grd))
  expect_gt(nrow(grd), 0L)
  expect_false(all(is.na(grd$.pred)))
})

test_that("area_of_applicability accepts an rf_fit and its importance weights", {
  skip_if_not_installed("ranger")
  dat <- mk_rf_pts(150)
  fit <- fit_rf_model(dat, "z", c("a", "b", "noise"), num_trees = 200L)
  nd  <- sf::st_as_sf(
    data.frame(x = c(100, 200), y = c(100, 200),
               a = c(0, 20), b = c(0, 20), noise = c(0, 0)),
    coords = c("x", "y"), crs = 32632
  )
  res <- area_of_applicability(nd, model = fit,
                               weights = pmax(fit$info$importance, 0))
  expect_s3_class(res, "aoa")
  expect_identical(res$predictor_vars, c("a", "b", "noise"))
  expect_true(res$aoa$AOA[1])
  expect_false(res$aoa$AOA[2])
})

test_that("negative permutation importance produces an actionable message", {
  # Raw ranger permutation importance goes slightly negative for predictors
  # that do not help, so the naive `weights = fit$info$importance` will hit
  # this; the message has to say what to do about it.
  expect_error(.aoa_weight_vector(c(a = 1, b = -0.01), c("a", "b")),
               "pmax\\(importance, 0\\)")
})
