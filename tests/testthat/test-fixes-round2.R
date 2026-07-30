# ---------------------------------------------------------------------------
# Regression tests for the second review pass:
#  - predict.gwr_fit() must align newdata to the training CRS
#  - prep_model_data() must coerce multi-vertex MULTIPOINT to POINT
#  - residual_morans_i() must accept Matrix-package (sparse) weights
#  - core-count handling must survive parallel::detectCores() returning NA
# ---------------------------------------------------------------------------

test_that("prep_model_data coerces multi-vertex MULTIPOINT to POINT", {
  skip_if_not_installed("sf")

  g <- sf::st_sfc(
    sf::st_multipoint(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10))),
    sf::st_multipoint(rbind(c(100, 100), c(120, 100))),
    sf::st_point(c(50, 50)),
    crs = 32632
  )
  dat <- sf::st_sf(y = c(1, 2, 3), x1 = c(10, 20, 30), geometry = g)

  out <- prep_model_data(dat, "y", "x1")

  # One row per feature, all plain POINT (no per-vertex explosion)
  expect_equal(nrow(out), 3L)
  gtypes <- as.character(sf::st_geometry_type(out, by_geometry = TRUE))
  expect_true(all(gtypes == "POINT"))

  # MULTIPOINT features collapse to their centroid
  xy <- sf::st_coordinates(out)
  expect_equal(unname(xy[1, ]), c(5, 5))
  expect_equal(unname(xy[2, ]), c(110, 100))
  expect_equal(unname(xy[3, ]), c(50, 50))

  # st_coordinates() now aligns 1:1 with attribute rows
  expect_equal(nrow(xy), nrow(out))
})


test_that("residual_morans_i accepts sparse Matrix weights (no kNN fallback)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("Matrix")

  set.seed(21)
  n <- 30
  pts <- sf::st_as_sf(
    data.frame(x = runif(n), y = runif(n), resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )
  resid_vec <- rnorm(n)
  fake_fit <- structure(
    list(data_sf = pts, residuals = resid_vec, engine = list()),
    class = "spatial_fit"
  )
  registerS3method("residuals", "spatial_fit",
                   function(object, ...) object$residuals)

  # Row-standardised weight matrix (3 neighbours each, weight 1/3)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nbrs <- sample(setdiff(seq_len(n), i), 3)
    W[i, nbrs] <- 1 / 3
  }
  W_sparse <- Matrix::Matrix(W, sparse = TRUE)

  res_dense  <- residual_morans_i(fake_fit, weights = W)
  res_sparse <- residual_morans_i(fake_fit, weights = W_sparse)

  # Sparse weights must be USED (identical results), not silently replaced
  # by the kNN fallback (which would give different I / sd / z).
  expect_equal(res_sparse$observed, res_dense$observed, tolerance = 1e-12)
  expect_equal(res_sparse$sd,       res_dense$sd,       tolerance = 1e-12)
  expect_equal(res_sparse$z,        res_dense$z,        tolerance = 1e-12)
  expect_equal(res_sparse$p_value,  res_dense$p_value,  tolerance = 1e-12)
})


test_that(".sanitize_core_count collapses NA/invalid input to the fallback", {
  f <- spatialkit:::.sanitize_core_count

  expect_identical(f(NA), 1L)
  expect_identical(f(NA_integer_ - 1L), 1L)   # detectCores() NA arithmetic
  expect_identical(f(0), 1L)
  expect_identical(f(-3), 1L)
  expect_identical(f(NULL), 1L)
  expect_identical(f("not a number"), 1L)
  expect_identical(f(7), 7L)
  expect_identical(f(3.9), 3L)
  expect_identical(f(NA, fallback = 2L), 2L)
})


test_that(".resolve_n_cores returns sane values for the sequential paths", {
  f <- spatialkit:::.resolve_n_cores

  expect_identical(f(parallel = FALSE), 1L)
  expect_identical(f(parallel = 1), 1L)
  expect_identical(f(parallel = NA), 1L)

  skip_on_os("windows")
  expect_identical(f(parallel = 3), 3L)
  expect_identical(f(parallel = FALSE, n_cores = 2), 2L)
})


test_that("fit_bayesian_spatial_model rejects MULTIPOINT when .already_prepped", {
  skip_if_not_installed("sf")
  skip_if_not_installed("brms")

  g <- sf::st_sfc(
    sf::st_multipoint(rbind(c(0, 0), c(10, 0))),
    sf::st_point(c(5, 5)),
    sf::st_point(c(2, 8)),
    crs = 32632
  )
  dat <- sf::st_sf(y = c(1, 2, 3), x1 = c(10, 20, 30), geometry = g)

  expect_error(
    fit_bayesian_spatial_model(dat, "y", "x1", .already_prepped = TRUE),
    "must be POINT"
  )
})


test_that("predict.gwr_fit aligns newdata CRS to the training CRS", {
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")
  skip_if_not_installed("GWmodel")

  set.seed(31)
  n <- 40
  px <- runif(n, 0, 5000)
  py <- runif(n, 0, 5000)
  df <- data.frame(
    y  = 2 + 0.001 * px + 0.5 * rnorm(n),
    x1 = rnorm(n),
    px = px, py = py
  )
  train_sf <- sf::st_as_sf(df, coords = c("px", "py"), crs = 32632)

  fit <- fit_gwr_model(train_sf, "y", "x1",
                       adaptive = TRUE, bandwidth = 20)

  newdata <- train_sf[1:10, ]
  pred_same <- predict(fit, newdata = newdata)

  # Same locations expressed in lon/lat: before the fix, ensure_projected()
  # auto-picked a UTM zone from the newdata centroid (zone 31 here, not the
  # training zone 32), silently mixing coordinate systems in gwr.predict().
  newdata_ll <- sf::st_transform(newdata, 4326)
  pred_ll <- predict(fit, newdata = newdata_ll)
  expect_equal(pred_ll, pred_same, tolerance = 1e-4)

  # Same locations in a different *projected* CRS: before the fix this was
  # left untouched (already projected), producing garbage distances.
  newdata_merc <- sf::st_transform(newdata, 3857)
  pred_merc <- predict(fit, newdata = newdata_merc)
  expect_equal(pred_merc, pred_same, tolerance = 1e-4)
})
