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
  # A field with real short-range structure, so the variogram is identifiable.
  set.seed(5); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 25) + diag(1e-4, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)),
                      coords = c("x", "y"), crs = 3857)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  p <- plot(fit, type = "variogram")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
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
