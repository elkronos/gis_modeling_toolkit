# ===========================================================================
# GWR regressions from the adversarial review.
# ===========================================================================

skip_if_not_installed("GWmodel")
skip_if_not_installed("sp")


# ---------------------------------------------------------------------------
# POINT Z / M geometry.  GWmodel is strictly 2-D: gw.dist() refuses a third
# data coordinate and reshapes prediction points with matrix(, ncol = 2), so
# an elevation that survived into the sp object scrambled predict()'s
# locations, failed or silently 3-D-ified the fit, and put 3-D distances
# into gwr_model_selection().
# ---------------------------------------------------------------------------

# The same 120 stations as a 2-D layer and with an elevation (0-3000 m).
.rg_xyz <- function(n = 120, seed = 3) {
  set.seed(seed)
  x <- runif(n, 0, 2000); y <- runif(n, 0, 2000)
  df <- data.frame(x = 5e5 + x, y = 5e6 + y, zel = runif(n, 0, 3000),
                   a = rnorm(n), b = rnorm(n))
  df$v <- 1 + (x / 1000) * df$a + rnorm(n, 0, 0.3)
  list(d2 = sf::st_as_sf(df, coords = c("x", "y"), crs = 32632),
       d3 = sf::st_as_sf(df, coords = c("x", "y", "zel"), crs = 32632))
}

test_that("prep_model_data() returns XY points for POINT Z and POINT M input", {
  xyz <- .rg_xyz()
  xy  <- sf::st_coordinates(prep_model_data(xyz$d2, "v", "a"))
  p3  <- prep_model_data(xyz$d3, "v", "a")
  expect_identical(colnames(sf::st_coordinates(p3)), c("X", "Y"))
  expect_identical(sf::st_coordinates(p3), xy)

  pm <- sf::st_sf(sf::st_drop_geometry(xyz$d2),
                  geometry = sf::st_sfc(lapply(seq_len(nrow(xy)), function(i)
                    sf::st_point(c(xy[i, ], 7), dim = "XYM")), crs = 32632))
  expect_identical(sf::st_coordinates(prep_model_data(pm, "v", "a")), xy)
})

test_that("predict() on POINT Z newdata predicts at the XY locations", {
  xyz <- .rg_xyz()
  fit <- suppressWarnings(fit_gwr_model(xyz$d2, "v", "a", bandwidth = 30))
  # A constant altitude (KML's 0) scrambled the locations just the same.
  xy <- sf::st_coordinates(xyz$d2)
  d0 <- sf::st_as_sf(data.frame(sf::st_drop_geometry(xyz$d2),
                                X = xy[, 1], Y = xy[, 2], Z = 0),
                     coords = c("X", "Y", "Z"), crs = 32632)
  expect_identical(predict(fit, newdata = xyz$d3), predict(fit, newdata = xyz$d2))
  expect_identical(predict(fit, newdata = d0), predict(fit, newdata = xyz$d2))
})

test_that("fit_gwr_model() fits POINT Z data on 2-D distances", {
  xyz <- .rg_xyz()
  f2 <- suppressWarnings(fit_gwr_model(xyz$d2, "v", "a", bandwidth = 30))
  f3 <- suppressWarnings(fit_gwr_model(xyz$d3, "v", "a", bandwidth = 30))
  expect_identical(coef(f3), coef(f2))
  expect_identical(f3$info$AICc, f2$info$AICc)

  # .already_prepped = TRUE skips prep_model_data(), so .to_sp() has to drop
  # the elevation itself -- for the fit and again for predict()'s training set.
  f3p <- suppressWarnings(fit_gwr_model(xyz$d3, "v", "a", bandwidth = 30,
                                        .already_prepped = TRUE))
  expect_identical(coef(f3p), coef(f2))
  expect_identical(predict(f3p, newdata = xyz$d2[1:10, ]),
                   predict(f2, newdata = xyz$d2[1:10, ]))
})

test_that("gwr_model_selection() and cv_gwr() use 2-D distances on POINT Z data", {
  xyz <- .rg_xyz()
  s3 <- suppressWarnings(gwr_model_selection(xyz$d3, "v", c("a", "b"), bandwidth = 30))
  s2 <- suppressWarnings(gwr_model_selection(xyz$d2, "v", c("a", "b"), bandwidth = 30))
  expect_identical(s3$table, s2$table)

  folds <- make_folds(xyz$d2, k = 3, method = "block_kfold", seed = 7)
  cv3 <- suppressWarnings(suppressMessages(
    cv_gwr(xyz$d3, "v", "a", folds = folds, bandwidth = 30)))
  cv2 <- suppressWarnings(suppressMessages(
    cv_gwr(xyz$d2, "v", "a", folds = folds, bandwidth = 30)))
  expect_equal(cv3$n_folds_succeeded, 3L)
  expect_identical(cv3$predictions, cv2$predictions)
})
