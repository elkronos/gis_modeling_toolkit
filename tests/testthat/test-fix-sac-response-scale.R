# tests/testthat/test-fix-sac-response-scale.R
# ---------------------------------------------------------------------------
# BUG-234: estimate_sac_range()'s structure test and semivariances, and
# kriging_adequacy()'s variogram, kriging and per-cell variances, squared a
# response near 1e160 or 1e-160 as given.  Near 1e160 the squares overflowed:
# no structure test, a range refused as "fewer than three lag bins hold point
# pairs", other folds under make_folds(auto_range = TRUE), and
# kriging_adequacy() refused for that reason.  Near 1e-160 they fell below
# the smallest normal double: Moran's I and the range moved, gstat kriged no
# cell and se moved.
#
# Scaling by a power of two (2^-530, 2^530) changes no bit of a mantissa, and
# every step of the fix scales by powers of two, so the unit-free results
# (range, Moran's I, p, folds) are compared with expect_identical() there.
# The products 1e-160 and 1e160 round each value once, to about 1e-16 of it;
# the results are compared at 1e-6 relative, far above that and far below the
# moves the bug made (1e-4 and more, or NA).
# ---------------------------------------------------------------------------

skip_if_not_installed("gstat")

scale_layer <- function() {
  set.seed(7)
  n <- 150
  x <- runif(n, 0, 3000); yy <- runif(n, 0, 3000)
  # On a grid of 2^-13, so the products below are exact or rounded once.
  z <- round((sin(2 * pi * x / 1000) + cos(2 * pi * yy / 1000) +
                0.3 * rnorm(n)) * 8192) / 8192
  pts <- sf::st_as_sf(data.frame(x = x, yy = yy, z = z), coords = c("x", "yy"),
                      crs = 32617)
  sq <- lapply(0:8, function(i) {
    a <- (i %% 3) * 1000; b <- (i %/% 3) * 1000
    sf::st_polygon(list(rbind(c(a, b), c(a + 1000, b), c(a + 1000, b + 1000),
                              c(a, b + 1000), c(a, b))))
  })
  cells <- sf::st_sf(poly_id = sprintf("C%d", 1:9), geometry = sf::st_sfc(sq, crs = 32617))
  pts$poly_id <- cells$poly_id[as.integer(sf::st_within(pts, cells))]
  list(pts = pts, cells = cells)
}

times <- function(d, s) { d$z <- d$z * s; d }

sac_of <- function(d) suppressMessages(estimate_sac_range(d, "z", seed = 1))

test_that("the structure test reads the response at any size", {
  L <- scale_layer()
  xy <- sf::st_coordinates(L$pts)
  base <- spatialkit:::.sac_structure_test(xy, L$pts$z)
  for (s in c(2^-530, 2^530)) {
    st <- spatialkit:::.sac_structure_test(xy, L$pts$z * s)
    expect_false(is.null(st))
    expect_identical(st$table, base$table)
    expect_identical(st$p, base$p)
  }
  for (s in c(1e-160, 1e160)) {
    st <- spatialkit:::.sac_structure_test(xy, L$pts$z * s)
    expect_false(is.null(st))
    expect_equal(st$table$morans_i, base$table$morans_i, tolerance = 1e-6)
  }
})

test_that("estimate_sac_range() gives the same range, test and reason at any size", {
  L <- scale_layer()
  base <- sac_of(L$pts)
  expect_true(is.finite(base))
  for (s in c(2^-530, 2^530)) {
    r <- sac_of(times(L$pts, s))
    expect_identical(as.numeric(r), as.numeric(base))
    expect_identical(attr(r, "morans_i"), attr(base, "morans_i"))
    expect_identical(attr(r, "structure_p"), attr(base, "structure_p"))
    expect_identical(attr(r, "rejected_reason"), attr(base, "rejected_reason"))
    expect_identical(attr(r, "short_range"), attr(base, "short_range"))
  }
  for (s in c(1e-160, 1e160)) {
    r <- sac_of(times(L$pts, s))
    expect_true(is.finite(r))
    expect_equal(as.numeric(r), as.numeric(base), tolerance = 1e-6)
    expect_equal(attr(r, "morans_i"), attr(base, "morans_i"), tolerance = 1e-6)
    expect_null(attr(r, "rejected_reason"))
  }
})

test_that("a refused range keeps its own reason at any size", {
  L <- scale_layer()
  # A cutoff of a fifth of the extent leaves the 1000 m structure past the
  # fitted lags: refused as such, not as "fewer than three lag bins".
  ref <- function(d) suppressMessages(estimate_sac_range(d, "z", cutoff = 0.2, seed = 1))
  base <- ref(L$pts)
  why <- attr(base, "rejected_reason")
  expect_true(is.character(why))
  for (s in c(2^-530, 2^530, 1e-160, 1e160)) {
    r <- ref(times(L$pts, s))
    expect_identical(attr(r, "rejected_reason"), why)
    expect_equal(attr(r, "rejected_range"), attr(base, "rejected_range"), tolerance = 1e-6)
  }
})

test_that("the variances come back in the response's units squared", {
  L <- scale_layer()
  base <- sac_of(L$pts)
  # 2^-500: the standard deviation (about 2e-151) is below 1e-150, so the
  # variogram is fitted in units of a power of two, and the variances
  # (about 1e-301) are still normal doubles.  Put back by powers of two,
  # they are the base's times 2^-1000.
  s <- 2^-500
  r <- sac_of(times(L$pts, s))
  for (a in c("nugget", "short_lag_semivariance", "nugget_floor"))
    expect_equal(attr(r, a) / s^2, attr(base, a), tolerance = 1e-12, info = a)
  expect_equal(attr(r, "variogram")$gamma / s^2, attr(base, "variogram")$gamma,
               tolerance = 1e-12)
  expect_equal(attr(r, "variogram_model")$psill / s^2,
               attr(base, "variogram_model")$psill, tolerance = 1e-12)
  expect_equal(sac_nugget(r) / s^2, sac_nugget(base), tolerance = 1e-12)
  # Near 1e160 the squares are past the largest double: Inf, as documented.
  big <- sac_of(times(L$pts, 2^530))
  expect_identical(attr(big, "short_lag_semivariance"), Inf)
  expect_identical(sum(attr(big, "variogram_model")$psill), Inf)
  expect_true(all(is.infinite(attr(big, "variogram")$gamma)))
})

test_that("make_folds(auto_range = TRUE) builds the same folds at any size", {
  L <- scale_layer()
  mf <- function(d) suppressMessages(make_folds(d, k = 4, method = "block_kfold", seed = 1,
                                                auto_range = TRUE, response_var = "z"))
  base <- mf(L$pts)$assignment
  for (s in c(2^-530, 2^530, 1e-160, 1e160)) {
    a <- mf(times(L$pts, s))$assignment
    expect_identical(a$fold[order(a$row_id)], base$fold[order(base$row_id)])
  }
})

test_that("kriging_adequacy() kriges every cell at any size, with se and kr_pred in the response's units", {
  L <- scale_layer()
  ka <- function(d) suppressMessages(kriging_adequacy(d, "z", L$cells, k = 3L, nmax = 20L))
  base <- ka(L$pts)
  expect_true(all(is.finite(base$kr_pred)))
  for (s in c(2^-530, 2^530, 1e-160, 1e160)) {
    k <- ka(times(L$pts, s))
    expect_true(all(is.finite(k$kr_pred)))
    expect_equal(k$se / s, base$se, tolerance = 1e-6)
    expect_equal(k$mean / s, base$mean, tolerance = 1e-6)
    expect_equal(k$kr_pred / s, base$kr_pred, tolerance = 1e-6)
    expect_equal(k$kr_ratio, base$kr_ratio, tolerance = 1e-6)
    expect_equal(k$kr_gain, base$kr_gain, tolerance = 1e-6)
    expect_equal(k$kr_shift, base$kr_shift, tolerance = 1e-6)
    expect_equal(attr(k, "range"), attr(base, "range"), tolerance = 1e-6)
    expect_equal(attr(k, "cv")$zscore_var, attr(base, "cv")$zscore_var, tolerance = 1e-6)
    expect_equal(attr(k, "cv")$rmse / s, attr(base, "cv")$rmse, tolerance = 1e-6)
  }
  # At 2^-500 the variances are normal doubles: the base's times 2^-1000.
  s <- 2^-500
  k <- ka(times(L$pts, s))
  expect_equal(k$kr_var / s^2, base$kr_var, tolerance = 1e-6)
  expect_equal(k$plain_var / s^2, base$plain_var, tolerance = 1e-6)
  expect_equal(attr(k, "sill") / s^2, attr(base, "sill"), tolerance = 1e-6)
  expect_equal(attr(k, "nugget") / s^2, attr(base, "nugget"), tolerance = 1e-6)
})

test_that("kriging_adequacy() reads a supplied variogram at any size", {
  L <- scale_layer()
  ka <- function(d) suppressMessages(kriging_adequacy(d, "z", L$cells, sac = sac_of(d),
                                                      k = 3L, nmax = 20L))
  base <- ka(L$pts)
  # 2^-500: the supplied model's variances are normal doubles, so dividing
  # them by the square of the unit loses nothing.
  s <- 2^-500
  k <- ka(times(L$pts, s))
  expect_true(all(is.finite(k$kr_pred)))
  expect_equal(k$kr_pred / s, base$kr_pred, tolerance = 1e-6)
  expect_equal(k$kr_var / s^2, base$kr_var, tolerance = 1e-6)
  # The model is returned as it was passed.
  expect_identical(attr(k, "variogram")$psill,
                   attr(sac_of(times(L$pts, s)), "variogram_model")$psill)
  # At 2^-530 the supplied model's variances are subnormal and carry fewer
  # bits; every cell is still kriged, within the precision they hold.
  s <- 2^-530
  k <- ka(times(L$pts, s))
  expect_true(all(is.finite(k$kr_pred)))
  expect_equal(k$kr_pred / s, base$kr_pred, tolerance = 1e-2)
})

test_that("a cell whose own spread is near 1e-160 gets its se from .cell_sd()", {
  L <- scale_layer()
  d <- L$pts
  # One cell's values are its deviations from their mean times 2^-530
  # (near 1e-160); the layer's spread stays ordinary, so only that cell's
  # variance underflows, and var() of it is off in its fifth digit.
  in1 <- d$poly_id == "C1"
  d$z[in1] <- (d$z[in1] - mean(d$z[in1])) * 2^-530
  k <- suppressWarnings(suppressMessages(
    kriging_adequacy(d, "z", L$cells, k = 3L, nmax = 20L)))
  zc <- d$z[in1] - spatialkit:::.regression_centre(d$z, bounded = FALSE)
  # As a ratio: expect_equal() compares numbers this small absolutely.
  expect_equal(k$se[k$poly_id == "C1"] / (spatialkit:::.cell_sd(zc) / sqrt(sum(in1))), 1,
               tolerance = 1e-12)
  expect_gt(k$se[k$poly_id == "C1"], 0)
})

# A layer with repeat visits: the first 60 locations twice, the second visit
# off by up to about 0.6 (on the same 2^-13 grid), so that the pairs at
# distance zero hold 60 degrees of freedom and the nugget has a floor.
visit_layer <- function() {
  L <- scale_layer()
  d <- rbind(L$pts, L$pts[1:60, ])
  set.seed(11)
  d$z[151:210] <- d$z[151:210] + round(0.2 * rnorm(60) * 8192) / 8192
  list(pts = d, cells = L$cells)
}

test_that("the nugget floor and the zero-lag semivariance come back times the scale squared", {
  V <- visit_layer()
  base <- sac_of(V$pts)
  expect_gt(attr(base, "nugget_floor"), 0)
  expect_gt(attr(base, "zero_lag_semivariance"), 0)
  s <- 2^-500
  r <- sac_of(times(V$pts, s))
  expect_identical(as.numeric(r), as.numeric(base))
  for (a in c("nugget_floor", "zero_lag_semivariance", "nugget"))
    expect_equal(attr(r, a) / s^2, attr(base, a), tolerance = 1e-12, info = a)
})

test_that("the directional fits come back in the response's units squared", {
  L <- scale_layer()
  dir_of <- function(d) suppressMessages(
    estimate_sac_range(d, "z", seed = 1, keep_directional_fits = TRUE))
  base <- attr(dir_of(L$pts), "directional_fits")
  expect_true(is.list(base) && length(base) > 0L)
  s <- 2^-500
  got <- attr(dir_of(times(L$pts, s)), "directional_fits")
  expect_identical(names(got), names(base))
  for (a in names(base)) {
    expect_equal(got[[a]]$variogram$gamma / s^2, base[[a]]$variogram$gamma,
                 tolerance = 1e-12, info = a)
    expect_equal(got[[a]]$model$psill / s^2, base[[a]]$model$psill,
                 tolerance = 1e-12, info = a)
  }
})

test_that("kriging_adequacy()'s record of the replicates is in the response's units squared", {
  V <- visit_layer()
  ka <- function(d) suppressWarnings(suppressMessages(
    kriging_adequacy(d, "z", V$cells, k = 3L, nmax = 20L)))
  base <- attr(ka(V$pts), "replicates")
  expect_gt(base$within_var, 0)
  s <- 2^-500
  got <- attr(ka(times(V$pts, s)), "replicates")
  for (a in c("within_var", "nugget", "nugget_used"))
    expect_equal(got[[a]] / s^2, base[[a]], tolerance = 1e-6, info = a)
  expect_identical(got$n_points, base$n_points)
  expect_identical(got$df, base$df)
  expect_identical(got$above_nugget, base$above_nugget)
})

# A layer with a predictor, for the REML detrending.
reml_layer <- function() {
  set.seed(7)
  n <- 150
  x <- runif(n, 0, 3000); yy <- runif(n, 0, 3000)
  a <- rnorm(n)
  z <- round((sin(2 * pi * x / 1000) + cos(2 * pi * yy / 1000) + 0.5 * a +
                0.3 * rnorm(n)) * 8192) / 8192
  sf::st_as_sf(data.frame(x = x, yy = yy, z = z, a = a), coords = c("x", "yy"),
               crs = 32617)
}

test_that("detrend = \"reml\" fits the response at any size", {
  skip_if_not_installed("nlme")
  d <- reml_layer()
  # The warnings each call raises, counted; near 1e160 the REML fit failed
  # and the range came from the OLS fallback, with a warning that REML "did
  # not converge".
  est <- function(e) {
    w <- character(0)
    r <- withCallingHandlers(
      suppressMessages(estimate_sac_range(e, "z", predictor_vars = "a",
                                          detrend = "reml", seed = 1)),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
      })
    list(r = r, w = w)
  }
  base <- est(d)
  expect_length(base$w, 0L)
  expect_false(is.null(attr(base$r, "reml")))
  expect_true(is.finite(base$r))
  for (s in c(2^-530, 2^530)) {
    o <- est(times(d, s))
    expect_length(o$w, 0L)
    expect_false(is.null(attr(o$r, "reml")))
    expect_identical(as.numeric(o$r), as.numeric(base$r))
    expect_identical(attr(o$r, "reml")$nugget_prop, attr(base$r, "reml")$nugget_prop)
  }
  for (s in c(1e-160, 1e160)) {
    o <- est(times(d, s))
    expect_length(o$w, 0L)
    expect_false(is.null(attr(o$r, "reml")))
    expect_equal(as.numeric(o$r), as.numeric(base$r), tolerance = 1e-6)
  }
  # At 2^-500 the model's variances are normal doubles: the base's times
  # 2^-1000, and so is the REML variance.
  s <- 2^-500
  o <- est(times(d, s))
  expect_equal(attr(o$r, "variogram_model")$psill / s^2,
               attr(base$r, "variogram_model")$psill, tolerance = 1e-12)
  expect_equal(attr(o$r, "reml")$sigma2 / s^2, attr(base$r, "reml")$sigma2,
               tolerance = 1e-12)
  expect_equal(attr(o$r, "reml")$logLik_gain, attr(base$r, "reml")$logLik_gain,
               tolerance = 1e-9)
  # Near 1e160 the variances are past the largest double, as documented.
  o <- est(times(d, 2^530))
  expect_identical(sum(attr(o$r, "variogram_model")$psill), Inf)
})

test_that("kriging_adequacy() refuses a supplied model whose sill overflowed, and says so", {
  L <- scale_layer()
  d <- times(L$pts, 2^530)
  sac <- sac_of(d)
  expect_identical(sum(attr(sac, "variogram_model")$psill), Inf)
  expect_error(suppressMessages(kriging_adequacy(d, "z", L$cells, sac = sac,
                                                 k = 3L, nmax = 20L)),
               "sill past the largest double")
  # sac = NULL estimates the variogram on the scaled response and kriges.
  k <- suppressMessages(kriging_adequacy(d, "z", L$cells, k = 3L, nmax = 20L))
  expect_true(all(is.finite(k$kr_pred)))
})
