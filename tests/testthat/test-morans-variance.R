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

  # The one-sided alternatives: the same z, one tail each.  No test passed
  # `alternative =` at all, so "greater" returning the TWO-sided p went
  # unnoticed (mutation testing, pass 6).
  p_g <- residual_morans_i(fake_fit, weights = W, alternative = "greater")$p_value
  p_l <- residual_morans_i(fake_fit, weights = W, alternative = "less")$p_value
  expect_equal(p_g, stats::pnorm(res$z, lower.tail = FALSE), tolerance = 1e-12)
  expect_equal(p_l, stats::pnorm(res$z, lower.tail = TRUE),  tolerance = 1e-12)
  expect_equal(p_g + p_l, 1, tolerance = 1e-12)
  expect_equal(2 * min(p_g, p_l), res$p_value, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(p_g, res$p_value)))
})


# Cliff & Ord (1981, sec. 8.3) moments of I for OLS residuals, written with an
# explicit dense M = I - X (X'X)^-1 X' -- the textbook form, independent of the
# package's sparse trace reductions.
.cliff_ord_residual_moments <- function(W, X) {
  n  <- nrow(X); p <- ncol(X); S0 <- sum(W)
  M  <- diag(n) - X %*% solve(crossprod(X), t(X))
  MW <- M %*% W
  trMW <- sum(diag(MW))
  EI <- (n / S0) * trMW / (n - p)
  VI <- (n / S0)^2 * (sum(diag(MW %*% M %*% t(W))) + sum(diag(MW %*% MW)) + trMW^2) /
    ((n - p) * (n - p + 2)) - EI^2
  list(EI = EI, VI = VI)
}

test_that("the regression-residual null is the Cliff & Ord moments, to machine precision", {
  # Every residual_morans_i() test used a stub with no response_var /
  # predictor_vars, so .morans_ols_design() returned NULL and the
  # randomisation null was taken on every path: the Cliff & Ord residual
  # moments -- E[I] = (n/S0) tr(MW)/(n-p) and the matching variance -- were
  # never evaluated, and null = "auto" could never have selected them.  Two
  # independent references: the dense-M textbook formula above, always, and
  # spdep::lm.morantest() -- the agreement the documentation claims -- when
  # spdep is installed.
  set.seed(3)
  n <- 80
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  w <- 0.002 * x + rnorm(n)                  # a spatially smooth covariate
  z <- 0.5 * w + rnorm(n)                    # independent errors
  d <- sf::st_as_sf(data.frame(x = x, y = y, w = w, z = z),
                    coords = c("x", "y"), crs = 3857)
  fit <- lm_spatial_fit(d, response_var = "z", predictor_vars = "w")
  W   <- moran_stub_weights(n, k = 6, seed = 3)

  auto <- residual_morans_i(fit, weights = W)
  expect_identical(auto$null, "residual")   # "auto" recognises OLS residuals
  expect_equal(auto$df, n - 2)

  res <- residual_morans_i(fit, weights = W, null = "residual")
  expect_identical(res$null, "residual")
  expect_equal(res[c("observed", "expected", "sd", "z", "p_value")],
               auto[c("observed", "expected", "sd", "z", "p_value")])

  # The dense-M reference.
  X   <- cbind(1, w)
  co  <- .cliff_ord_residual_moments(W, X)
  e   <- stats::residuals(fit$engine); e <- e - mean(e)
  expect_equal(res$observed, (n / sum(W)) * sum(e * (W %*% e)) / sum(e^2), tolerance = 1e-12)
  expect_equal(res$expected, co$EI, tolerance = 1e-12)
  expect_equal(res$sd, sqrt(co$VI), tolerance = 1e-12)
  expect_equal(res$z, (res$observed - co$EI) / sqrt(co$VI), tolerance = 1e-12)
  expect_equal(res$p_value, 2 * stats::pnorm(abs(res$z), lower.tail = FALSE), tolerance = 1e-12)

  if (requireNamespace("spdep", quietly = TRUE)) {
    # W is row-standardised already, so style = "W" hands spdep the same
    # matrix (and avoids its warning about style "M").
    expect_equal(unname(rowSums(W)), rep(1, n))
    lw  <- spdep::mat2listw(W, style = "W")
    ref <- spdep::lm.morantest(fit$engine, listw = lw, alternative = "two.sided")
    expect_equal(res$observed, unname(ref$estimate[1]), tolerance = 1e-10)
    expect_equal(res$expected, unname(ref$estimate[2]), tolerance = 1e-10)
    expect_equal(res$sd,  sqrt(unname(ref$estimate[3])), tolerance = 1e-10)
    expect_equal(res$z,   as.numeric(ref$statistic),     tolerance = 1e-10)
    expect_equal(res$p_value, as.numeric(ref$p.value),   tolerance = 1e-10)
  }

  # The two nulls are distinguishable on this fit: a smooth covariate drags
  # E[I] below -1/(n - 1), which is why "auto" matters.
  rnd <- residual_morans_i(fit, weights = W, null = "randomisation")
  expect_identical(rnd$null, "randomisation")
  expect_equal(rnd$expected, -1 / (n - 1))
  expect_false(isTRUE(all.equal(rnd$expected, res$expected)))
  expect_false(isTRUE(all.equal(rnd$p_value, res$p_value)))

  # Intercept-only: tr(MW) = -S0/n for any zero-diagonal W, so the general
  # formula must reduce to the classical -1/(n - 1) exactly.
  fit0 <- lm_spatial_fit(d, response_var = "z", predictor_vars = character(0))
  res0 <- residual_morans_i(fit0, weights = W, null = "residual")
  expect_identical(res0$null, "residual")
  expect_equal(res0$expected, -1 / (n - 1), tolerance = 1e-12)
  co0 <- .cliff_ord_residual_moments(W, matrix(1, n, 1))
  expect_equal(res0$sd, sqrt(co0$VI), tolerance = 1e-12)
  if (requireNamespace("spdep", quietly = TRUE)) {
    ref0 <- spdep::lm.morantest(fit0$engine, listw = lw, alternative = "two.sided")
    expect_equal(res0$sd, sqrt(unname(ref0$estimate[3])), tolerance = 1e-10)
    expect_equal(res0$p_value, as.numeric(ref0$p.value), tolerance = 1e-10)
  }
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
