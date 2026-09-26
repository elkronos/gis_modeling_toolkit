# tests/testthat/test-review2-predict-surface.R
# ---------------------------------------------------------------------------
# predict_surface(): a polygon grid, `draws` passed through `...`, the grid's
# floating-point column count, and a stale .pred_se on a reused surface.
# Uses the lm-backed spatial_fit from helper-lmfit.R.
# ---------------------------------------------------------------------------

# Covariate points on a 20-unit lattice offset by 3, so no grid-cell centre is
# equidistant from two of them and the nearest one is unambiguous.
r2_lattice <- function() {
  xy <- expand.grid(x = seq(3, 1000, by = 20), y = seq(3, 1000, by = 20))
  pts <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 3857, remove = FALSE)
  pts$t <- pts$y / 200
  set.seed(2)
  pts$z <- 3 * pts$t + stats::rnorm(nrow(pts), 0, 0.1)
  pts
}


test_that("a polygon grid takes covariates at each cell's point and returns points", {
  cov <- r2_lattice()
  fit <- lm_spatial_fit(cov, predictor_vars = "t")
  polys <- sf::st_sf(geometry = sf::st_make_grid(
    sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000),
                              crs = sf::st_crs(3857))), n = c(5, 5)))

  s_poly <- predict_surface(fit, grid = polys, covariates = cov)
  s_pts  <- predict_surface(fit, grid = coerce_to_points(polys, "auto"),
                            covariates = cov)
  expect_true(all(sf::st_geometry_type(s_poly) == "POINT"))
  expect_identical(nrow(s_poly), 25L)
  # The cell's covariate is the one nearest its representative point, not an
  # arbitrary point inside it: the bottom-left cell's centre is at y = 100, so
  # t is about 0.5, where the old path could hand it anything from 0 to 1.
  expect_equal(s_poly$t, s_pts$t)
  expect_equal(s_poly$.pred, s_pts$.pred)
  expect_true(all(abs(s_poly$t - sf::st_coordinates(s_poly)[, 2] / 200) < 0.1))

  # And the answer no longer depends on the row order of `covariates`.
  set.seed(4)
  shuffled <- cov[sample(nrow(cov)), ]
  expect_equal(predict_surface(fit, grid = polys, covariates = shuffled)$.pred,
               s_poly$.pred)
})


test_that("`draws` passed through `...` is refused rather than flattened into .pred", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)
  expect_error(predict_surface(fit, n_cells = 50, draws = TRUE),
               "`draws` cannot be passed through `...`")
  expect_error(predict_surface(fit, n_cells = 50, dr = TRUE),
               "`draws` cannot be passed through `...`")
  expect_error(predict_surface(fit, n_cells = 50, se = TRUE, draws = TRUE),
               "`draws` cannot be passed through `...`")
})


test_that("a backend returning the wrong number of values is an error", {
  pts <- surf_test_points()
  fit <- lm_spatial_fit(pts)
  class(fit) <- c("r2_badlen_fit", class(fit))
  registerS3method("predict", "r2_badlen_fit",
                   function(object, newdata = NULL, ...) rep(1, 2L * nrow(newdata)))
  expect_error(predict_surface(fit, n_cells = 50),
               "returned \\d+ value\\(s\\) for the \\d+ rows")
})


test_that("the automatic grid does not lose a column to floating-point rounding", {
  grid_fn <- spatialkit:::.make_prediction_grid
  crs <- sf::st_crs(32632)
  # 0.3 / 0.1 is 2.9999999999999996, and floor() dropped the third column.
  bb <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 0.3, ymax = 0.3), crs = crs)
  g <- grid_fn(bb, crs, cell_size = 0.1)
  expect_equal(sort(unique(g$..grid_x)), c(0.05, 0.15, 0.25))
  expect_identical(nrow(g), 9L)

  # At the default n_cells a square's cell size is its width / 100, and the
  # same rounding gave 99 columns for some widths (2.11 among them).
  for (w in c(2.11, 2.48, 3.59, 11.36)) {
    bb <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = w, ymax = w), crs = crs)
    g <- grid_fn(bb, crs, n_cells = 10000)
    expect_identical(length(unique(g$..grid_x)), 100L)
    expect_identical(length(unique(g$..grid_y)), 100L)
  }

  # A width that is NOT a multiple keeps the floor: 100 / 30 is 3 cells.
  bb <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 100, ymax = 100), crs = crs)
  expect_equal(sort(unique(grid_fn(bb, crs, cell_size = 30)$..grid_x)),
               c(15, 45, 75))
})


test_that("a reused surface does not keep the previous model's .pred_se", {
  pts <- surf_test_points()
  fit_a <- lm_spatial_fit(pts)
  fit_b <- lm_spatial_fit(pts, predictor_vars = "w")
  s1 <- predict_surface(fit_a, n_cells = 100, se = TRUE, covariates = pts)
  expect_true(".pred_se" %in% names(s1))

  s2 <- predict_surface(fit_b, grid = s1, covariates = pts)
  expect_false(".pred_se" %in% names(s2))
  expect_false(isTRUE(all.equal(s2$.pred, s1$.pred)))

  # Asked for again, it is this model's.
  s3 <- predict_surface(fit_b, grid = s1, covariates = pts, se = TRUE)
  expect_true(".pred_se" %in% names(s3))
  expect_false(isTRUE(all.equal(s3$.pred_se, s1$.pred_se)))
})
