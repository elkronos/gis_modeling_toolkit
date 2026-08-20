# tests/testthat/test-predict-surface.R
# ---------------------------------------------------------------------------
# predict_surface(): grid construction, clipping, covariate joins, chunking.
#
# Exercised through a deterministic lm-backed spatial_fit (helper-lmfit.R) so
# the whole path runs without GWmodel or Stan.
# ---------------------------------------------------------------------------

test_that("a grid covers the training extent at roughly the requested count", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)

  s <- predict_surface(fit, n_cells = 400)
  expect_s3_class(s, "sf")
  expect_true(".pred" %in% names(s))
  expect_true(all(is.finite(s$.pred)))
  # Cell count is approximate -- the grid is built from a derived cell size.
  expect_gt(nrow(s), 250)
  expect_lt(nrow(s), 600)

  bb_train <- sf::st_bbox(pts); bb_grid <- sf::st_bbox(s)
  expect_gte(bb_grid[["xmin"]], bb_train[["xmin"]])
  expect_lte(bb_grid[["xmax"]], bb_train[["xmax"]])
  expect_equal(sf::st_crs(s), sf::st_crs(pts))
})

test_that("cell_size takes precedence over n_cells", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)
  s <- predict_surface(fit, cell_size = 100, n_cells = 999999)
  expect_equal(attr(s, "cell_size"), 100)
  expect_lt(nrow(s), 200)          # ~10x10 over a 1000-unit extent
})

test_that("chunking does not change the result", {
  # The property that matters: rows do not interact, so chunk size is purely
  # a memory/throughput knob.
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)

  whole   <- predict_surface(fit, n_cells = 400, chunk_size = 10000)
  chunked <- predict_surface(fit, n_cells = 400, chunk_size = 37)

  expect_equal(nrow(whole), nrow(chunked))
  expect_equal(whole$.pred, chunked$.pred)
})

test_that("boundary clipping drops outside points", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)

  square <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(200, 200), c(600, 200), c(600, 600), c(200, 600), c(200, 200)))),
    crs = 3857)

  full <- predict_surface(fit, n_cells = 900)
  clip <- predict_surface(fit, n_cells = 900, boundary = square)

  expect_lt(nrow(clip), nrow(full))
  bb <- sf::st_bbox(clip)
  expect_gte(bb[["xmin"]], 200); expect_lte(bb[["xmax"]], 600)
  expect_gte(bb[["ymin"]], 200); expect_lte(bb[["ymax"]], 600)
})

test_that("a model with predictors demands covariates", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  expect_error(predict_surface(fit, n_cells = 100), "absent from the grid")
  expect_error(predict_surface(fit, n_cells = 100, covariates = pts["z"]),
               "lacks column")
})

test_that("covariates are joined from the nearest feature", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  s <- predict_surface(fit, n_cells = 200, covariates = pts)
  expect_true("w" %in% names(s))
  expect_true(all(is.finite(s$.pred)))
  # every joined value must come from the covariate layer
  expect_true(all(s$w %in% pts$w))
})

test_that("a supplied grid is used verbatim", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)
  g <- sf::st_as_sf(expand.grid(x = seq(100, 900, by = 200),
                                y = seq(100, 900, by = 200)),
                    coords = c("x", "y"), crs = 3857)
  s <- predict_surface(fit, grid = g)
  expect_equal(nrow(s), nrow(g))
})

test_that("se returns a posterior SD when the backend exposes draws", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)
  s <- predict_surface(fit, n_cells = 100, se = TRUE)
  expect_true(".pred_se" %in% names(s))
  expect_true(all(is.finite(s$.pred_se)))
  expect_true(all(s$.pred_se > 0))
})

test_that("non-spatial_fit input is rejected", {
  expect_error(predict_surface(lm(mpg ~ wt, mtcars)), "must be a spatial_fit")
})

test_that("an empty clip is an error, not an empty surface", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)
  far <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(1e6, 1e6), c(1e6 + 10, 1e6), c(1e6 + 10, 1e6 + 10),
    c(1e6, 1e6 + 10), c(1e6, 1e6)))), crs = 3857)
  expect_error(predict_surface(fit, n_cells = 100, boundary = far),
               "no grid points")
})


test_that("a cell_size wider than the extent yields one cell, not an error", {
  # seq(from, to, by) errors with "wrong sign in 'by' argument" when from > to
  # rather than returning length 0, so the previous length-0 guard was
  # unreachable and an oversized cell_size aborted with an internal message.
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)
  bb  <- sf::st_bbox(pts)
  huge <- 5 * max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
  g <- predict_surface(fit, cell_size = huge, covariates = pts)
  expect_true(inherits(g, "sf"))
  expect_identical(nrow(g), 1L)
})

test_that("grid spacing is exactly cell_size and stays inside the bbox", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)
  g   <- predict_surface(fit, cell_size = 100, covariates = pts)
  xy  <- sf::st_coordinates(g)
  ux  <- sort(unique(xy[, 1]))
  bb  <- sf::st_bbox(pts)
  expect_true(all(abs(diff(ux) - 100) < 1e-8))
  expect_gte(min(ux), as.numeric(bb[["xmin"]]))
  expect_lte(max(ux), as.numeric(bb[["xmax"]]))
})
