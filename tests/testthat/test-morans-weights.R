# ===========================================================================
# Spatial weight matrices for residual_morans_i().
#
# Two failure modes, both silent: a sparse matrix being dropped in favour of
# the kNN fallback (giving a different statistic than the one asked for),
# and non-row-standardised weights being treated as though they were.
# See test-morans-variance.R for the variance formula itself.
#
# The fit stub and its residuals method live in helper-moranstub.R, on a
# test-only subclass -- see the note there on why registering against
# "spatial_fit" itself was wrong.
# ===========================================================================

test_that("residual_morans_i accepts sparse Matrix weights (no kNN fallback)", {
  skip_if_not_installed("Matrix")

  set.seed(21)
  n <- 30
  pts <- sf::st_as_sf(
    data.frame(x = runif(n), y = runif(n), resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )
  fake_fit <- moran_stub_fit(pts, rnorm(n))

  # Row-standardised weight matrix (3 neighbours each, weight 1/3)
  W <- moran_stub_weights(n, k = 3, seed = 21)
  W_sparse <- Matrix::Matrix(W, sparse = TRUE)

  res_dense  <- residual_morans_i(fake_fit, weights = W)
  res_sparse <- residual_morans_i(fake_fit, weights = W_sparse)

  # Sparse weights must be USED (identical results), not silently replaced
  # by the kNN fallback (which would give different I / sd / z).
  expect_equal(res_sparse$observed, res_dense$observed, tolerance = 1e-12)
  expect_equal(res_sparse$sd,       res_dense$sd,       tolerance = 1e-12)
  expect_equal(res_sparse$z,        res_dense$z,        tolerance = 1e-12)
  expect_equal(res_sparse$p_value,  res_dense$p_value,  tolerance = 1e-12)

  # ... and the kNN fallback really would give something else, so the four
  # equalities above are not satisfiable by "both paths fell back".
  res_knn <- residual_morans_i(fake_fit)
  expect_false(isTRUE(all.equal(res_knn$observed, res_dense$observed)))
})


test_that("residual_morans_i warns on non-row-standardised user weights", {
  # No skip_if_not_installed("FNN"): this test always passes explicit
  # `weights =`, so the kNN builder is never reached at all.

  n <- 20
  set.seed(42)
  pts <- sf::st_as_sf(
    data.frame(x = runif(n), y = runif(n), resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )
  fake_fit <- moran_stub_fit(pts, rnorm(n))

  # A row-standardised matrix should NOT trigger the row-standardisation
  # warning.  .log_warn() raises no R condition, so this has to be asserted on
  # captured log lines -- expect_no_warning() would pass vacuously either way
  # (see helper-logging.R).
  W_good <- moran_stub_weights(n, k = 4, seed = 1, value = 1 / 4)
  good_lines <- capture_spatialkit_log(
    result <- residual_morans_i(fake_fit, weights = W_good)
  )
  expect_false(log_has(good_lines, "not row-standardised"))
  expect_true(is.list(result))
  expect_true(is.finite(result$observed))

  # A non-row-standardised (binary) matrix SHOULD emit the logger warning.
  W_bad <- moran_stub_weights(n, k = 4, seed = 1, value = 1)
  bad_lines <- capture_spatialkit_log(
    res_bad <- residual_morans_i(fake_fit, weights = W_bad)
  )
  expect_true(log_has(bad_lines, "not row-standardised"))
  # The message reports the observed row-sum range, which is exactly 4 here.
  expect_true(log_has(bad_lines, "row sums range from 4 to 4"))

  # The warning is advisory: the general Cliff & Ord formula is still applied,
  # so a usable statistic comes back.  W_bad is 4x W_good, and both I and the
  # variance are invariant to a positive rescaling of W, so the two agree.
  expect_true(is.finite(res_bad$observed))
  expect_equal(res_bad$observed, result$observed, tolerance = 1e-10)
  expect_equal(res_bad$z,        result$z,        tolerance = 1e-10)
})
