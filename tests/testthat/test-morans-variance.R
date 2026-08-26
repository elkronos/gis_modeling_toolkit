# ---------------------------------------------------------------------------
# Regression test for residual_morans_i() variance (Cliff & Ord formula)
#
# Guards against the denominator bug where an extra factor of n made Var(I)
# ~n x too small, inflating z-scores by ~sqrt(n).
#
# The fit stub and its residuals method live in helper-moranstub.R, on a
# test-only subclass -- see the note there on why registering against
# "spatial_fit" itself was wrong.
# ---------------------------------------------------------------------------

test_that("residual_morans_i variance matches the Cliff & Ord randomisation formula", {
  set.seed(7)
  n <- 25
  coords_mat <- cbind(runif(n), runif(n))
  pts <- sf::st_as_sf(
    data.frame(x = coords_mat[, 1], y = coords_mat[, 2], resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )
  resid_vec <- rnorm(n)

  fake_fit <- moran_stub_fit(pts, resid_vec)

  # Row-standardised weight matrix (4 neighbours each, weight 1/4)
  W <- moran_stub_weights(n, k = 4, seed = 7)

  res <- residual_morans_i(fake_fit, weights = W)
  expect_true(is.list(res))

  # --- Independent computation, written directly from Cliff & Ord (1981) ---
  e  <- resid_vec - mean(resid_vec)
  S0 <- sum(W)
  S1 <- 0.5 * sum((W + t(W))^2)
  S2 <- sum((rowSums(W) + colSums(W))^2)
  b2 <- (sum(e^4) / n) / ((sum(e^2) / n)^2)

  I_exp <- (n / S0) * sum(e * (W %*% e)) / sum(e^2)
  EI    <- -1 / (n - 1)
  VI    <- (n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2) -
              b2 * ((n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2)) /
           ((n - 1) * (n - 2) * (n - 3) * S0^2) - EI^2

  expect_equal(res$observed, I_exp, tolerance = 1e-10)
  expect_equal(res$expected, EI,    tolerance = 1e-12)
  expect_equal(res$sd, sqrt(VI),    tolerance = 1e-10)
  expect_equal(res$z, (I_exp - EI) / sqrt(VI), tolerance = 1e-10)

  # p-value must be the two-sided normal tail of the correct z
  expect_equal(res$p_value,
               2 * stats::pnorm(abs(res$z), lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("residual_morans_i does not flag white-noise residuals as significant", {
  # No skip_if_not_installed("FNN"): .build_knn_weights() has a dense
  # order()-based fallback that is refused only above n = 5000, and n = 100
  # here.  Guarding this on FNN skipped the very test the file exists for
  # wherever FNN was absent -- which is every CI job but one.

  # With the extra-n bug, z was inflated ~sqrt(n) and pure noise came out
  # "significant" almost always. With the correct variance, white noise
  # should rarely be significant.
  set.seed(123)
  n <- 100
  pts <- sf::st_as_sf(
    data.frame(x = runif(n), y = runif(n), resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )
  fake_fit <- moran_stub_fit(pts, rnorm(n))

  res <- residual_morans_i(fake_fit)
  expect_true(is.finite(res$z))
  # |z| for iid noise should be modest; the buggy version produced |z| ~ 10+
  expect_lt(abs(res$z), 4)
})

test_that("the default kNN weights agree whether or not FNN is available", {
  # The dense fallback must build the same row-standardised k-NN matrix as
  # the kd-tree path; otherwise the statistic silently depends on which
  # optional packages happen to be installed.
  set.seed(5)
  n <- 40
  coords <- cbind(runif(n), runif(n))

  W_dense <- spatialkit:::.build_knn_weights(coords, k = 6L,
                                             use_fnn = FALSE, use_matrix = FALSE)
  expect_equal(dim(W_dense), c(n, n))
  expect_equal(unname(rowSums(W_dense)), rep(1, n), tolerance = 1e-12)
  expect_equal(sum(diag(W_dense)), 0)               # no self-neighbours

  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  W_fast <- spatialkit:::.build_knn_weights(coords, k = 6L,
                                            use_fnn = TRUE, use_matrix = TRUE)
  expect_equal(as.matrix(W_fast), W_dense, tolerance = 1e-12,
               ignore_attr = TRUE)
})
