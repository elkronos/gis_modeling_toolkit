# The GP pre-fit check of fit_bayesian_spatial_model() (.gp_resid_moran(),
# behind .gp_short_lag_moran() and .gp_floor_moran()) and the
# standardisation of predictors and coordinates, at any origin and unit.
#
# The check regressed the raw response on the raw predictors: a predictor at
# 1e8 was dropped by the QR as collinear with the intercept, and a response
# at 1e15 left the rounding of the offset in the residuals, so I, z and the
# check's status moved with the origin of a column.  Its only test for
# residuals without spread was sum(e^2) <= 0, so it scored the rounding of
# an exact fit that estimate_sac_range() refuses, and could log that the
# response "is autocorrelated at short range".  standardize_predictors = TRUE
# floored a predictor's scale at .Machine$double.eps, so a predictor in small
# units (a standard deviation under 2.2e-16) was not standardised, and the
# coordinates' scales had the same floor.
#
# No Stan: brms::brm() is replaced by a recorder.

# 600 points; a predictor wq and a response zq with a spatial field and a
# trend on wq.  Both are the values as a double holds them at 1e15, so
# off + wq and off + zq are exact at every offset used here and hold the
# same spread to the last bit.
.fgpo_layer <- function(n = 600, seed = 3, wave = c(2000, 1500)) {
  withr::with_seed(seed, {
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4); w <- rnorm(n)
    z <- 2 * w + sin(x / wave[1]) + cos(y / wave[2]) + rnorm(n, sd = 0.5)
  })
  list(xy = cbind(x, y), wq = (1e15 + w) - 1e15, zq = (1e15 + z) - 1e15)
}

.fgpo_stat <- function(r) if (is.null(r)) NULL else c(I = r$I, z = r$z, p = r$p)

# .gp_resid_moran() as it was: the raw fit and the raw design.
.fgpo_old <- function(y, X, W) {
  n  <- length(y)
  e  <- stats::lm.fit(X, y)$residuals
  e  <- e - mean(e)
  S0 <- sum(W)
  if (!is.finite(S0) || S0 <= 0 || sum(e^2) <= 0) return(NULL)
  I   <- as.numeric((n / S0) * sum(e * (W %*% e)) / sum(e^2))
  mom <- spatialkit:::.morans_residual_moments(W = W, X = X, S0 = S0,
                                               is_sparse = inherits(W, "Matrix"))
  if (is.null(mom) || !is.finite(mom$VI) || mom$VI <= 0) return(NULL)
  z <- (I - mom$EI) / sqrt(mom$VI)
  list(I = I, z = z, p = stats::pnorm(z, lower.tail = FALSE))
}

.fgpo_capture <- function(...) {
  cap <- new.env()
  local_mocked_bindings(
    brm = function(...) { cap$args <- list(...); structure(list(), class = "fgpo_stub") },
    .package = "brms")
  fit <- suppressWarnings(suppressMessages(fit_bayesian_spatial_model(
    ..., compute_loo = FALSE, check_convergence = FALSE)))
  list(fit = fit, args = cap$args)
}


test_that("a predictor or a response at an offset gives the Moran's I of the same columns at 0 (BUG-94)", {
  L <- .fgpo_layer()
  sl <- function(y, X) .fgpo_stat(spatialkit:::.gp_short_lag_moran(L$xy, y, X))
  fl <- function(y, X) .fgpo_stat(spatialkit:::.gp_floor_moran(L$xy, y, X, floor = 2000))
  X0 <- cbind(p = L$wq)
  ref_s <- sl(L$zq, X0); ref_f <- fl(L$zq, X0)
  expect_gt(ref_s[["I"]], 0.7)
  for (o in c(1e8, 1e15)) {
    expect_equal(sl(L$zq, cbind(p = o + L$wq)), ref_s, tolerance = 1e-8)
    expect_equal(fl(L$zq, cbind(p = o + L$wq)), ref_f, tolerance = 1e-8)
  }
  for (o in c(1e12, 1e15)) {
    expect_equal(sl(o + L$zq, X0), ref_s, tolerance = 1e-8)
    expect_equal(fl(o + L$zq, X0), ref_f, tolerance = 1e-8)
  }
  # Both at once, and a response whose raw values show their spread through
  # a trend on the predictor (fitted again about its centre).
  expect_equal(sl(1e15 + L$zq, cbind(p = 1e8 + L$wq)), ref_s, tolerance = 1e-8)
  expect_equal(sl(1e15 + 1e5 * L$wq + L$zq, X0), sl(1e5 * L$wq + L$zq, X0),
               tolerance = 1e-8)
  # A trend too large for the centring rule (the centre under 1,000 times
  # the spread): the residuals of the raw fit lose their spread, and the fit
  # is taken again about the response's centre.  The fit's own rounding of
  # terms of 1e12 bounds the agreement (the raw fit at 1e15 gave I 0.722
  # for 0.770).
  expect_equal(sl(1e15 + 1e12 * L$wq + L$zq, X0), sl(1e12 * L$wq + L$zq, X0),
               tolerance = 1e-4)
  # No predictors.
  expect_equal(sl(1e15 + L$zq, NULL), sl(L$zq, NULL), tolerance = 1e-8)
  expect_equal(fl(1e15 + L$zq, NULL), fl(L$zq, NULL), tolerance = 1e-8)
})

test_that("the pre-fit check's status and tests do not move with a column's origin (BUG-94)", {
  # A layer whose variogram range (about 1,435) is below the floor (3,000).
  L <- .fgpo_layer(seed = 5, wave = c(300, 300))
  chk <- function(p, z) {
    d <- sf::st_as_sf(data.frame(x = L$xy[, 1], y = L$xy[, 2], z = z, p = p),
                      coords = c("x", "y"), crs = 32615)
    suppressWarnings(suppressMessages(spatialkit:::.gp_prefit_range_check(
      d, "z", "p", resolution = 3000, gp_c = 1.5, S = 1, xy_scale = 1)))
  }
  ref <- chk(L$wq, L$zq)
  expect_identical(ref$status, "below_floor")
  for (r in list(chk(1e8 + L$wq, L$zq), chk(L$wq, 1e15 + L$zq))) {
    expect_identical(r$status, ref$status)
    expect_equal(r$range, ref$range, tolerance = 1e-8)
    expect_equal(.fgpo_stat(r$moran), .fgpo_stat(ref$moran), tolerance = 1e-8)
    expect_equal(.fgpo_stat(r$moran_floor), .fgpo_stat(ref$moran_floor), tolerance = 1e-8)
  }
})

test_that("a constant or collinear predictor at an offset is still dropped (BUG-94)", {
  L <- .fgpo_layer()
  ref <- .fgpo_stat(spatialkit:::.gp_short_lag_moran(L$xy, L$zq, cbind(p = L$wq)))
  for (X in list(cbind(p = L$wq, c = 1e8),
                 cbind(p = 1e8 + L$wq, q = 1e8 + 2 * L$wq),
                 cbind(p = L$wq, q = 5 + 3 * L$wq))) {
    expect_equal(.fgpo_stat(spatialkit:::.gp_short_lag_moran(L$xy, L$zq, X)), ref,
                 tolerance = 1e-8)
  }
})

test_that("columns below the centring gate are tested as before, to the last bit (BUG-94, BUG-92)", {
  for (s in 1:6) {
    withr::with_seed(s, {
      n <- c(60, 200, 600)[(s %% 3) + 1]
      xy <- cbind(runif(n, 0, 1e4) + 5e5 * (s %% 2), runif(n, 0, 1e4) + 4.5e6 * (s %% 2))
      w <- rnorm(n, 10 * (s %% 2)); w2 <- runif(n, 0, 100)
      z <- (50 + 2 * w + sin(xy[, 1] / 700) + cos(xy[, 2] / 900) + rnorm(n, sd = 0.3 * s)) *
        10^c(-6, 0, 6)[(s %% 3) + 1]
    })
    for (X in list(matrix(1, n, 1L), cbind(1, w), cbind(1, w, w2), cbind(1, w, c = 3))) {
      Wk <- spatialkit:::.build_knn_weights(xy, k = 4L)
      expect_identical(spatialkit:::.gp_resid_moran(z, X, Wk), .fgpo_old(z, X, Wk))
      D <- as.matrix(stats::dist(xy))
      Wb <- (D >= 1125 & D < 1875) * 1
      diag(Wb) <- 0
      expect_identical(spatialkit:::.gp_resid_moran(z, X, Wb), .fgpo_old(z, X, Wb))
    }
  }
  # A trend of 1e10 on a spread of 1: the residuals' spread is 3e-11 of the
  # response's size, above estimate_sac_range()'s floor of 1e-12, so the fit
  # is not taken again.
  withr::with_seed(11, {
    n <- 300
    xy <- cbind(runif(n, 0, 1e4), runif(n, 0, 1e4))
    w <- rnorm(n); z <- 1e10 * w + sin(xy[, 1] / 700) + rnorm(n, sd = 0.3)
  })
  Wk <- spatialkit:::.build_knn_weights(xy, k = 4L)
  expect_identical(spatialkit:::.gp_resid_moran(z, cbind(1, w), Wk),
                   .fgpo_old(z, cbind(1, w), Wk))
})

# 300 points and two predictors; exact fits v = a + b w + c w2.
.fgpo_exact <- function() {
  withr::with_seed(21, {
    n <- 300
    xy <- cbind(5e5 + runif(n, 0, 5000), 4.5e6 + runif(n, 0, 5000))
    w <- rnorm(n); w2 <- rnorm(n)
  })
  cfs <- withr::with_seed(99, matrix(rnorm(3000, 0, 5), ncol = 3, byrow = TRUE))
  list(xy = xy, w = w, w2 = w2, cfs = cfs,
       mk = function(v) sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], v = v,
                                                w = w, w2 = w2),
                                     coords = c("X", "Y"), crs = 32615))
}

test_that("the residuals of an exact fit give no Moran test, at any scale (BUG-92)", {
  E <- .fgpo_exact()
  X <- cbind(w = E$w, w2 = E$w2)
  for (i in c(1:20, 384)) {
    cf <- E$cfs[i, ]
    v <- cf[1] + cf[2] * E$w + cf[3] * E$w2
    for (k in c(-60, 0, 60)) {
      expect_null(spatialkit:::.gp_short_lag_moran(E$xy, 2^k * v, X))
      expect_null(spatialkit:::.gp_floor_moran(E$xy, 2^k * v, X, floor = 1000))
    }
  }
  # Fit 384's rounding gave p = 0.0059; the check now has no test to report.
  cf <- E$cfs[384, ]
  g <- suppressWarnings(suppressMessages(spatialkit:::.gp_prefit_range_check(
    E$mk(cf[1] + cf[2] * E$w + cf[3] * E$w2), "v", c("w", "w2"),
    resolution = c(x = 400, y = 400), gp_c = 1.5, S = 1, xy_scale = 1)))
  expect_identical(g$status, "not_estimated")
  expect_match(g$reason, "explain the response exactly")
  expect_null(g$moran)
  expect_null(g$moran_floor)
  # Residuals with spread are scored, the same at every scale.
  v <- withr::with_seed(7, 1 + 2 * E$w - E$w2 + rnorm(length(E$w)))
  ref <- spatialkit:::.gp_short_lag_moran(E$xy, v, X)
  expect_false(is.null(ref))
  for (k in c(-1000, -600, -60, 60, 900)) {
    expect_identical(spatialkit:::.gp_short_lag_moran(E$xy, 2^k * v, X), ref)
    expect_identical(spatialkit:::.gp_floor_moran(E$xy, 2^k * v, X, floor = 1000),
                     spatialkit:::.gp_floor_moran(E$xy, v, X, floor = 1000))
  }
})

test_that("an exact fit at an offset, or with terms larger than the response, gives no Moran test (BUG-92)", {
  E <- .fgpo_exact()
  none <- function(v, X) {
    expect_null(spatialkit:::.gp_short_lag_moran(E$xy, v, X))
    expect_null(spatialkit:::.gp_floor_moran(E$xy, v, X, floor = 1000))
  }
  for (i in 1:40) {
    cf <- E$cfs[i, ]
    # The response at an offset: the rounding is that of the response as
    # given, not of the centred one.
    none(1e6 + cf[1] + cf[2] * E$w + cf[3] * E$w2, cbind(w = E$w, w2 = E$w2))
    none(1e10 + cf[1] + cf[2] * E$w + cf[3] * E$w2, cbind(w = E$w, w2 = E$w2))
    # A predictor at an offset.
    p8 <- 1e8 + E$w
    none(cf[1] + cf[2] * p8 + cf[3] * E$w2, cbind(p = p8, w2 = E$w2))
  }
  # Terms of 4.5e6 that cancel to a response of order 1: the rounding is
  # that of the predictors as given, with their centres put back.
  p <- 4.5e6 + E$w; q <- 4.5e6 + E$w2
  for (i in 1:60) {
    b <- E$cfs[i, 2]
    none(3 + b * p - b * q, cbind(p = p, q = q))
  }
  b <- E$cfs[1, 2]
  d <- sf::st_as_sf(data.frame(X = E$xy[, 1], Y = E$xy[, 2], v = 3 + b * p - b * q,
                               p = p, q = q), coords = c("X", "Y"), crs = 32615)
  g <- suppressWarnings(suppressMessages(spatialkit:::.gp_prefit_range_check(
    d, "v", c("p", "q"), resolution = c(x = 400, y = 400), gp_c = 1.5, S = 1,
    xy_scale = 1)))
  expect_identical(g$status, "not_estimated")
  expect_null(g$moran)
  expect_null(g$moran_floor)
})

test_that("the log says nothing about autocorrelation of an exact fit's rounding (BUG-92)", {
  skip_if_not_installed("brms")
  E <- .fgpo_exact()
  cf <- E$cfs[384, ]
  logged <- capture_spatialkit_log(
    r <- .fgpo_capture(E$mk(cf[1] + cf[2] * E$w + cf[3] * E$w2), "v", c("w", "w2"),
                       chains = 1, iter = 200, warmup = 100),
    level = logger::WARN)
  expect_identical(r$fit$info$gp_prefit$status, "not_estimated")
  expect_null(r$fit$info$gp_prefit$moran)
  expect_true(log_has(logged, "could not be estimated before fitting \\(the residuals on predictor_vars"))
  expect_false(log_has(logged, "autocorrelated at short range"))
  expect_false(log_has(logged, "Moran's I = "))
})

test_that("a predictor in small units is standardised to a standard deviation of 1 (BUG-91)", {
  skip_if_not_installed("brms")
  withr::with_seed(4, {
    n <- 120
    xy <- cbind(5e5 + runif(n, 0, 5000), 4.5e6 + runif(n, 0, 5000))
    w <- rnorm(n); y <- 1 + 2 * w + sin(xy[, 1] / 1500) + rnorm(n, 0, 0.5)
  })
  run <- function(p) .fgpo_capture(
    sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], y = y, p = p),
                 coords = c("X", "Y"), crs = 32615),
    "y", "p", standardize_predictors = TRUE, chains = 1, iter = 200, warmup = 100)
  ref <- run(w)
  ps0 <- ref$fit$info$predictor_scaling$p
  expect_identical(ps0$scale, sd(w, na.rm = TRUE))
  expect_identical(ps0$center, mean(w, na.rm = TRUE))
  b0 <- as.data.frame(ref$args$prior)
  for (k in c(-56, -60, -70, -600)) {
    r <- run(w * 2^k)
    expect_equal(sd(r$args$data$p), 1, tolerance = 1e-12)
    expect_equal(r$args$data$p, ref$args$data$p, tolerance = 1e-12)
    ps <- r$fit$info$predictor_scaling$p
    expect_equal(ps$scale / 2^k, ps0$scale, tolerance = 1e-12)
    expect_equal(ps$center / 2^k, ps0$center, tolerance = 1e-12)
    pr <- as.data.frame(r$args$prior)
    expect_identical(pr$prior[pr$class == "b"], b0$prior[b0$class == "b"])
  }
})

test_that("a predictor without spread keeps a scale of 1, and the coordinates read the same way (BUG-91)", {
  skip_if_not_installed("brms")
  withr::with_seed(4, {
    n <- 120
    xy <- cbind(runif(n, 0, 5000), runif(n, 0, 5000))
    w <- rnorm(n); y <- 1 + 2 * w + sin(xy[, 1] / 1500) + rnorm(n, 0, 0.5)
  })
  c1 <- rep(c(0.1, 0.3 / 3), length.out = n)   # one unit in the last place apart
  r <- .fgpo_capture(
    sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], y = y, p = w, c0 = 3, c1 = c1),
                 coords = c("X", "Y"), crs = 32615),
    "y", c("p", "c0", "c1"), standardize_predictors = TRUE,
    chains = 1, iter = 200, warmup = 100)
  ps <- r$fit$info$predictor_scaling
  expect_identical(ps$c0$scale, 1)
  expect_identical(ps$c1$scale, 1)
  expect_lt(max(abs(r$args$data$c1)), 1e-15)
  expect_identical(ps$p$scale, sd(w, na.rm = TRUE))

  # Coordinates whose standard deviation is under 2.2e-16 map units.
  run <- function(k) .fgpo_capture(
    sf::st_as_sf(data.frame(X = xy[, 1] * 2^k, Y = xy[, 2] * 2^k, y = y, p = w),
                 coords = c("X", "Y"), crs = 32615),
    "y", "p", chains = 1, iter = 200, warmup = 100)
  r0 <- run(0)
  cs0 <- r0$fit$info$coord_scaling
  expect_identical(cs0$x_scale, sd(xy[, 1]))
  expect_identical(cs0$y_center, mean(xy[, 2]))
  r70 <- run(-70)
  cs <- r70$fit$info$coord_scaling
  expect_equal(c(cs$x_scale, cs$y_scale) / 2^-70, c(cs0$x_scale, cs0$y_scale),
               tolerance = 1e-12)
  expect_equal(sd(r70$args$data$..x), 1, tolerance = 1e-12)
  expect_equal(sd(r70$args$data$..y), 1, tolerance = 1e-12)
  expect_equal(r70$fit$info$gp_resolution / 2^-70, r0$fit$info$gp_resolution,
               tolerance = 1e-8)
})

test_that("a coordinate axis without spread takes the other axis's scale (BUG-91)", {
  skip_if_not_installed("brms")
  withr::with_seed(8, {
    n <- 150
    yy <- sort(runif(n, 0, 3))    # a transect 3 map units long
    z <- sin(yy * 4) + rnorm(n, sd = 0.3)
  })
  run <- function(u) .fgpo_capture(
    sf::st_as_sf(data.frame(X = rep(1000, n) * u, Y = yy * u, z = z),
                 coords = c("X", "Y"), crs = 32615),
    "z", character(0), chains = 1, iter = 200, warmup = 100)
  r <- run(1)
  cs <- r$fit$info$coord_scaling
  expect_identical(cs$y_scale, sd(yy))
  expect_identical(cs$x_scale, cs$y_scale)
  expect_true(all(r$args$data$..x == 0))
  res <- r$fit$info$gp_resolution
  expect_identical(res[["x"]], res[["y"]])
  expect_identical(r$fit$info$gp_prefit$floor, res[["y"]])
  # The same transect in units 1,000 times smaller: the resolution and the
  # floor scale with the unit.
  r3 <- run(1000)
  expect_equal(r3$fit$info$gp_resolution / 1000, res, tolerance = 1e-12)
  expect_equal(r3$fit$info$gp_prefit$floor / 1000, r$fit$info$gp_prefit$floor,
               tolerance = 1e-12)
  expect_identical(r3$fit$info$gp_k, r$fit$info$gp_k)
})

test_that(".standardise_by() reads a column again only below eps", {
  sb <- spatialkit:::.standardise_by
  x <- c(1, 4, 9, 16)
  expect_identical(sb(x, mean(x), sd(x)), list(center = mean(x), scale = sd(x)))
  expect_identical(sb(rep(2, 5), 2, 0), list(center = 2, scale = 1))
  y <- x * 2^-80
  expect_equal(sb(y, mean(y), sd(y))$scale, sd(x) * 2^-80, tolerance = 1e-14)
  z <- x * 2^-1000
  expect_equal(sb(z, mean(z), sd(z))$scale / 2^-1000, sd(x), tolerance = 1e-14)
  big <- c(1, -1, 0.5) * 1e308
  out <- sb(big, mean(big), sd(big))
  expect_true(is.finite(out$scale) && out$scale > 1e307)
  expect_true(is.finite(out$center))
  # A mean and a standard deviation the caller could not take (a sum that
  # overflows): both are read again.
  huge <- c(1.7e308, 1.6e308, 1.5e308)
  out <- sb(huge, Inf, NaN)
  expect_equal(out$center, 1.6e308, tolerance = 1e-12)
  expect_equal(out$scale, 1e307, tolerance = 1e-12)
  expect_identical(sb(c(1.7e308, 1.7e308), Inf, NaN), list(center = 1.7e308, scale = 1))
  # A column without spread gets the scale `none`.
  expect_identical(sb(rep(2, 5), 2, 0, none = NA_real_), list(center = 2, scale = NA_real_))
})
