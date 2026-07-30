# ---------------------------------------------------------------------------
# Regression test for residual_morans_i() variance (Cliff & Ord formula)
#
# Guards against the denominator bug where an extra factor of n made Var(I)
# ~n x too small, inflating z-scores by ~sqrt(n).
# ---------------------------------------------------------------------------

test_that("residual_morans_i variance matches the Cliff & Ord randomisation formula", {
  skip_if_not_installed("sf")

  set.seed(7)
  n <- 25
  coords_mat <- cbind(runif(n), runif(n))
  pts <- sf::st_as_sf(
    data.frame(x = coords_mat[, 1], y = coords_mat[, 2], resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )
  resid_vec <- rnorm(n)

  fake_fit <- structure(
    list(data_sf = pts, residuals = resid_vec, engine = list()),
    class = "spatial_fit"
  )
  registerS3method("residuals", "spatial_fit",
                   function(object, ...) object$residuals)

  # Row-standardised weight matrix (4 neighbours each, weight 1/4)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nbrs <- sample(setdiff(seq_len(n), i), 4)
    W[i, nbrs] <- 1 / 4
  }

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
  skip_if_not_installed("sf")
  skip_if_not_installed("FNN")

  # With the extra-n bug, z was inflated ~sqrt(n) and pure noise came out
  # "significant" almost always. With the correct variance, white noise
  # should rarely be significant.
  set.seed(123)
  n <- 100
  pts <- sf::st_as_sf(
    data.frame(x = runif(n), y = runif(n), resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )
  fake_fit <- structure(
    list(data_sf = pts, residuals = rnorm(n), engine = list()),
    class = "spatial_fit"
  )
  registerS3method("residuals", "spatial_fit",
                   function(object, ...) object$residuals)

  res <- residual_morans_i(fake_fit)
  expect_true(is.finite(res$z))
  # |z| for iid noise should be modest; the buggy version produced |z| ~ 10+
  expect_lt(abs(res$z), 4)
})
