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


# ---------------------------------------------------------------------------
# predict() at new locations.  GWmodel::gwr.predict() returned EVERY value as
# NA when one location's window was empty or singular (inv() threw for the
# whole call), when nrow(train) + nrow(newdata) > 10000 ('DM3.given' not
# found) and when nrow(train) > 5000 ("No regression point is fixed").
# ---------------------------------------------------------------------------

.rg_pts <- function(n = 150, seed = 9, extent = 1000) {
  set.seed(seed)
  x <- runif(n, 0, extent); y <- runif(n, 0, extent)
  df <- data.frame(x = 5e5 + x, y = 5e6 + y, a = rnorm(n))
  df$v <- 1 + (x / (extent / 2)) * df$a + rnorm(n, 0, 0.3)
  sf::st_as_sf(df, coords = c("x", "y"), crs = 32632)
}

test_that("one empty fixed-bandwidth window makes only its own prediction NA", {
  d   <- .rg_pts()
  fit <- suppressWarnings(fit_gwr_model(d, "v", "a", adaptive = FALSE,
                                        bandwidth = 300))
  # 20 locations inside the data and, 11th, one 400 m beyond its edge, where
  # no training point lies within the bandwidth.
  set.seed(1)
  nd <- sf::st_as_sf(
    data.frame(x = 5e5 + c(runif(10, 100, 900), 1400, runif(10, 100, 900)),
               y = 5e6 + c(runif(10, 100, 900), 500, runif(10, 100, 900)),
               a = rnorm(21)),
    coords = c("x", "y"), crs = 32632)
  p_in <- predict(fit, newdata = nd[-11, ])
  expect_false(anyNA(p_in))
  expect_warning(p_all <- predict(fit, newdata = nd),
                 "1 of 21 location\\(s\\) have no estimable local regression")
  expect_true(is.na(p_all[11]))
  # The others are the values they get without the bad location, to the bit.
  expect_identical(p_all[-11], p_in)
})

test_that("fixed-bandwidth block CV scores the held-out points it can reach", {
  d <- .rg_pts(n = 200, seed = 10, extent = 2000)
  folds <- make_folds(d, k = 4, method = "block_kfold", block_size = 700,
                      seed = 1)
  cv <- suppressWarnings(suppressMessages(
    cv_gwr(d, "v", "a", folds = folds, adaptive = FALSE, bandwidth = 400)))
  expect_equal(cv$n_folds_succeeded, cv$n_folds_attempted)
  expect_gt(cv$overall$n_pred, 0)
})

test_that("predict() fills a 100 x 100 grid (train + newdata > 10000 rows)", {
  d   <- .rg_pts(n = 150, seed = 4)
  fit <- suppressWarnings(fit_gwr_model(d, "v", "a", bandwidth = 40))
  gx  <- seq(5, 995, length.out = 100)
  set.seed(2)
  g <- sf::st_as_sf(data.frame(x = 5e5 + rep(gx, 100),
                               y = 5e6 + rep(gx, each = 100),
                               a = rnorm(1e4)),
                    coords = c("x", "y"), crs = 32632)
  p <- predict(fit, newdata = g)
  expect_length(p, 1e4)
  expect_false(anyNA(p))
  some <- c(1:50, 9951:1e4)
  expect_identical(p[some], predict(fit, newdata = g[some, ]))
})

test_that("predict() works with more than 5000 training rows", {
  d <- .rg_pts(n = 5200, seed = 12, extent = 5000)
  # A GWR fit of 5200 rows takes minutes; predict() needs only the data and
  # the settings, so hand it those.
  fit <- new_spatial_fit("gwr_fit", engine = list(), formula = v ~ a,
                         response_var = "v", predictor_vars = "a", data_sf = d,
                         info = list(bandwidth = 60, adaptive = TRUE,
                                     kernel = "bisquare"))
  nd <- d[c(7, 1000, 4000), ]
  nd$a <- c(-1, 0.5, 2)
  p <- predict(fit, newdata = nd)
  # Local least squares by hand: bisquare weights on the distance to the
  # 60th-nearest TRAINING point.
  xy <- sf::st_coordinates(d)
  X  <- cbind(1, d$a)
  hand <- vapply(1:3, function(i) {
    xy0 <- sf::st_coordinates(nd)[i, ]
    dd  <- sqrt((xy[, 1] - xy0[1])^2 + (xy[, 2] - xy0[2])^2)
    w   <- spatialkit:::.gw_kernel_weights(dd, 60, "bisquare", TRUE)
    b   <- solve(crossprod(X, w * X), crossprod(X, w * d$v))
    sum(c(1, nd$a[i]) * b)
  }, numeric(1))
  expect_equal(p, hand, tolerance = 1e-10)
})

test_that("predict() keeps fitted() at the training points and NA rows in place", {
  # Not regressions -- the gwr.predict() path got these right too -- but what
  # its replacement has to preserve.
  d <- .rg_pts(n = 120, seed = 3)
  for (ad in c(TRUE, FALSE)) {
    fit <- suppressWarnings(fit_gwr_model(d, "v", "a", adaptive = ad,
                                          bandwidth = if (ad) 30 else 400))
    expect_identical(predict(fit, newdata = d), unname(fitted(fit)))
  }
  nd <- d[1:12, ]
  nd$a[c(2, 9)] <- c(NA, Inf)
  p <- predict(fit, newdata = nd)
  expect_identical(which(is.na(p)), c(2L, 9L))
  expect_identical(p[-c(2, 9)], predict(fit, newdata = nd[-c(2, 9), ]))
})
