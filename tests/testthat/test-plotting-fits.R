# tests/testthat/test-plotting-fits.R
# ---------------------------------------------------------------------------
# plot.spatial_fit() and plot_folds().
#
# The package shipped print() and summary() for fitted models but no plot(),
# so the diagnostics most likely to reveal a problem -- residual structure, and
# whether blocks actually separate the folds -- had to be hand-rolled.
#
# ggplot objects are built lazily, so constructing one does not prove it will
# render.  ggplot2::ggplot_build() forces the computation, which is what these
# tests assert on.
# ---------------------------------------------------------------------------

test_that("residual and observed-predicted plots build", {
  skip_if_not_installed("ggplot2")
  fit <- lm_spatial_fit(surf_test_points(), predictor_vars = "w")

  for (ty in c("residuals", "observed_predicted")) {
    p <- plot(fit, type = ty)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))   # forces the layer computation
  }
})

test_that("the residual plot carries one point per observation", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 90)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")
  b <- ggplot2::ggplot_build(plot(fit, type = "residuals"))
  expect_equal(nrow(b$data[[1]]), 90)
})

test_that("observed = fitted + residual is what gets plotted", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 60)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")
  b <- ggplot2::ggplot_build(plot(fit, type = "observed_predicted"))
  pt_layer <- b$data[[2]]                        # layer 1 is the abline
  expect_equal(sort(pt_layer$x), sort(pts$z), tolerance = 1e-8)
})

test_that("an unknown plot type is rejected", {
  skip_if_not_installed("ggplot2")
  fit <- lm_spatial_fit(surf_test_points())
  # A bare expect_error() passes on ANY error, including a typo in the test
  # itself.  match.arg() gives a stable message that names the valid values,
  # so pin it -- and pin that the valid values are still the documented three.
  expect_error(plot(fit, type = "nonsense"), "'arg' should be one of")
  expect_error(plot(fit, type = "nonsense"), "residuals")
  expect_error(plot(fit, type = "nonsense"), "observed_predicted")
  expect_error(plot(fit, type = "nonsense"), "variogram")
})

test_that("a fit without geometry is rejected", {
  skip_if_not_installed("ggplot2")
  fit <- lm_spatial_fit(surf_test_points())
  fit$data_sf <- sf::st_drop_geometry(fit$data_sf)
  expect_error(plot(fit, type = "residuals"), "no training geometry")
})

test_that("the residual variogram plot builds", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  # A field with real structure at the lags the variogram resolves (the
  # default binning is cutoff/15 ~ 40 units here), so the range is
  # identifiable.  An earlier a = 25 field was flat by the first lag bin and
  # fitted only by luck.
  set.seed(5); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 50) + diag(1e-4, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)),
                      coords = c("x", "y"), crs = 3857)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  p <- plot(fit, type = "variogram")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  # The fitted model line and the range marker are on it.
  layers <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_true("GeomLine" %in% layers)
  expect_true("GeomVline" %in% layers)

  # Residuals with NO spatial structure give a flat, nugget-only variogram
  # that no model fits: every start is singular, or the optimiser never
  # converges.  That is the picture most worth seeing, and it used to be
  # refused with "could not be fitted; there may be too few finite
  # residuals": estimate_sac_range() now returns NA WITH the empirical
  # variogram attached whichever way the fit failed, and the plot draws the
  # points, with a model line only when a (non-converged) model exists.
  flat <- pts
  set.seed(6); flat$z <- rnorm(n)
  fit_flat <- lm_spatial_fit(flat, predictor_vars = "w")
  sac_flat <- suppressWarnings(estimate_sac_range(
    sf::st_sf(.resid = residuals(fit_flat), geometry = sf::st_geometry(flat)),
    ".resid"))
  expect_true(is.na(sac_flat))
  expect_s3_class(attr(sac_flat, "variogram"), "data.frame")
  expect_match(attr(sac_flat, "rejected_reason"),
               "no variogram model|did not converge")
  p2 <- plot(fit_flat, type = "variogram")
  expect_s3_class(p2, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p2))
  layers2 <- vapply(p2$layers, function(l) class(l$geom)[1L], character(1))
  expect_true("GeomPoint" %in% layers2)
  expect_false("GeomVline" %in% layers2)               # no range marker
  expect_equal("GeomLine" %in% layers2, !is.null(attr(sac_flat, "variogram_model")))
  expect_match(p2$labels$subtitle, "^No effective range")

  # The all-singular outcome, made deterministic: every fit.variogram() call
  # reports a singular model.  The NA carries the empirical variogram and no
  # model, and the plot says why there is no range.
  local_mocked_bindings(
    fit.variogram = function(object, model, ...) {
      model$psill <- c(0, 0); model$range <- c(0, 1)
      attr(model, "singular") <- TRUE
      model
    },
    .package = "gstat")
  sac_sing <- suppressWarnings(estimate_sac_range(
    sf::st_sf(.resid = residuals(fit_flat), geometry = sf::st_geometry(flat)),
    ".resid"))
  expect_true(is.na(sac_sing))
  expect_s3_class(attr(sac_sing, "variogram"), "data.frame")
  expect_null(attr(sac_sing, "variogram_model"))
  expect_match(attr(sac_sing, "rejected_reason"), "^no variogram model could be fitted")
  p3 <- plot(fit_flat, type = "variogram")
  expect_no_error(ggplot2::ggplot_build(p3))
  layers3 <- vapply(p3$layers, function(l) class(l$geom)[1L], character(1))
  expect_true("GeomPoint" %in% layers3)
  expect_false("GeomLine" %in% layers3)
  expect_match(p3$labels$subtitle, "no variogram model could be fitted")
})


test_that("plot_folds maps every point to a fold", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 100)
  f <- make_folds(pts, k = 4, method = "block_kfold", seed = 1)

  p <- plot_folds(f, pts)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$data[[1]]), 100)
  # colour is discrete by fold, so no more distinct colours than folds
  expect_lte(length(unique(b$data[[1]]$colour)), length(f$folds))
})

test_that("plot_folds rejects mismatched inputs", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 50)
  f <- make_folds(pts, k = 3, method = "block_kfold", seed = 1)

  expect_error(plot_folds(list(), pts), "returned by make_folds")
  expect_error(plot_folds(f, sf::st_drop_geometry(pts)), "must be an sf object")

  # Folds built from a different layer share no row ids.
  other <- surf_test_points(n = 50, seed = 99)
  other$..row_id <- seq_len(50) + 10000L
  expect_error(plot_folds(f, other), "no points matched")
})

test_that("plot_folds accepts a boundary", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 80)
  f <- make_folds(pts, k = 3, method = "block_kfold", seed = 1)
  bnd <- sf::st_as_sfc(sf::st_bbox(pts))
  p <- plot_folds(f, pts, boundary = bnd)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})


# --- regressions found in the full-package audit -----------------------------

test_that("a sill-less variogram is drawn rather than refused", {
  # estimate_sac_range() returns NA when the fitted range exceeds the longest
  # lag, and used to discard the variogram with it -- so plot() refused in
  # exactly the case worth looking at.
  skip_if_not_installed("gstat")
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(200)
  fit <- lm_spatial_fit(pts, predictor_vars = character(0))
  # z carries a strong linear x/y trend, so an intercept-only fit leaves it
  # all in the residuals: the variogram rises monotonically and never reaches a
  # sill.  The deterministic attribute check lives in test-sac-range.R; this
  # asserts the user-visible outcome, that a plot comes back.
  # gstat cannot converge on a variogram that never reaches a sill -- that is
  # the whole point of this case, so its warning is expected here.
  p <- suppressWarnings(plot(fit, type = "variogram"))
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_gt(nrow(b$data[[1]]), 0L)          # the empirical points are there
})

test_that("all-NA residuals error instead of painting a uniformly grey map", {
  # limits = c(Inf, -Inf) builds without complaint and maps every value to
  # na.value, which looks like a result.
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(40)
  fit <- lm_spatial_fit(pts)
  broken <- fit
  registerS3method("residuals", "broken_fit",
                   function(object, ...) rep(NA_real_, nrow(object$data_sf)))
  class(broken) <- c("broken_fit", class(fit))
  for (ty in c("residuals", "observed_predicted", "variogram"))
    expect_error(plot(broken, type = ty), "no finite residuals")
})
