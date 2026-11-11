# residual_morans_i() reads its response, design and residuals at any origin
# and in any unit.
#
# The test that the residuals are OLS residuals and the residual moments were
# built on the raw predictor and response columns: a predictor at 1e8 was
# dropped by the QR (the moments then had rank 1) and a response at 1e8 left
# the rounding of the offset in the reference residuals, so the null switched
# in silence and null = "residual" logged that the residuals were not OLS
# residuals.  The tolerance of that test floored sd(y) at eps, so the test
# passed residuals that are not OLS residuals in small units.  The test for
# residuals with no spread was the absolute floor sum(resid_c^2) < eps, which
# refused residuals of real spread in small units as "zero residual variance"
# and scored the rounding of an exact fit in large ones.

# A custom fit whose residuals and fitted values are those of an lm().
residuals.rms_lm_fit <- function(object, ...) as.numeric(stats::residuals(object$engine))
fitted.rms_lm_fit <- function(object, ...) as.numeric(stats::fitted(object$engine))
registerS3method("residuals", "rms_lm_fit", residuals.rms_lm_fit)
registerS3method("fitted", "rms_lm_fit", fitted.rms_lm_fit)

# A custom fit whose residuals are given.
residuals.rms_fixed_fit <- function(object, ...) object$info$res
fitted.rms_fixed_fit <- function(object, ...)
  sf::st_drop_geometry(object$data_sf)[[object$response_var]] - object$info$res
registerS3method("residuals", "rms_fixed_fit", residuals.rms_fixed_fit)
registerS3method("fitted", "rms_fixed_fit", fitted.rms_fixed_fit)

# 200 points with an exponential field z, and two predictors w, w2.
.rms_layer <- function(n = 200, seed = 11) {
  withr::with_seed(seed, {
    xy <- cbind(5e5 + runif(n, 0, 5000), 4.5e6 + runif(n, 0, 5000))
    D <- as.matrix(stats::dist(xy))
    z <- as.numeric(t(chol(exp(-D / 1500) + diag(1e-8, n))) %*% rnorm(n))
    w <- rnorm(n); w2 <- rnorm(n)
  })
  list(n = n, xy = xy, z = z, w = w, w2 = w2)
}

.rms_sf <- function(L, v, w = L$w, w2 = L$w2) {
  sf::st_as_sf(data.frame(X = L$xy[, 1], Y = L$xy[, 2], v = v, w = w, w2 = w2),
               coords = c("X", "Y"), crs = 32615)
}

# An OLS fit of v on the predictors.  centred = TRUE fits the columns about
# their means, whose residuals do not depend on the origin of any column;
# FALSE fits the raw columns, as a backend that does not centre them would.
.rms_lm <- function(L, v, pv = "w", w = L$w, centred = FALSE) {
  d  <- .rms_sf(L, v, w = w)
  df <- sf::st_drop_geometry(d)
  fml <- stats::as.formula(paste("v ~", paste(pv, collapse = " + ")))
  eng <- if (centred) {
    dc <- as.data.frame(lapply(df[c("v", pv)], function(c) c - mean(c)))
    stats::lm(fml, data = dc)
  } else stats::lm(fml, data = df)
  new_spatial_fit("rms_lm_fit", eng, fml, "v", pv, d)
}

.rms_fixed <- function(L, v, res, w = L$w) {
  new_spatial_fit("rms_fixed_fit", list(), v ~ w, "v", "w", .rms_sf(L, v, w = w),
                  info = list(res = res))
}

.rms_keys <- c("observed", "expected", "sd", "z", "p_value", "n", "null", "df",
               "kurtosis", "p", "exact")

.rms_quiet <- function(expr) {
  lines <- capture_spatialkit_log(out <- suppressWarnings(expr))
  attr(out, "log") <- lines
  out
}


test_that("a predictor at a large offset gives the null, moments and rank of the predictor at 0", {
  L <- .rms_layer()
  wq <- (1e15 + L$w) - 1e15
  v  <- 2 * wq + L$z
  for (nl in c("auto", "residual")) {
    ref <- residual_morans_i(.rms_lm(L, v, w = wq, centred = TRUE), null = nl)
    expect_identical(ref$null, "residual")
    expect_true(ref$exact)
    expect_identical(ref$p, 2L)
    for (o in c(1e8, 1e12, 1e15)) {
      got <- .rms_quiet(residual_morans_i(.rms_lm(L, v, w = o + wq, centred = TRUE),
                                          null = nl))
      expect_equal(unclass(got)[.rms_keys], unclass(ref)[.rms_keys], tolerance = 1e-8,
                   info = sprintf("predictor at %g, null = %s", o, nl))
      expect_false(log_has(attr(got, "log"), "not the OLS residuals"))
    }
  }
})


test_that("a response at a large offset gives the null and moments of the response at 0", {
  L <- .rms_layer()
  ys <- (1e15 + L$z + 0.5 * L$w) - 1e15
  # The same residuals at every offset: the OLS residuals of the stored spread.
  e0 <- as.numeric(stats::residuals(stats::lm(ys ~ L$w)))
  ref <- residual_morans_i(.rms_fixed(L, ys, e0))
  expect_identical(ref$null, "residual")
  expect_true(ref$exact)
  for (o in c(1e6, 1e8, 1e9, 1e10, 1e12, 1e15)) {
    for (nl in c("auto", "residual")) {
      got <- .rms_quiet(residual_morans_i(.rms_fixed(L, o + ys, e0), null = nl))
      expect_equal(unclass(got)[.rms_keys], unclass(ref)[.rms_keys], tolerance = 1e-8,
                   info = sprintf("response at %g, null = %s", o, nl))
      expect_false(log_has(attr(got, "log"), "not the OLS residuals"))
    }
  }
})


test_that("a backend that regressed the raw response at an offset still gives OLS residuals", {
  L <- .rms_layer()
  ys <- (1e15 + L$z + 0.5 * L$w) - 1e15
  for (po in c(0, 1e5)) for (o in c(0, 1e6, 1e8, 1e10, 1e12)) {
    r <- residual_morans_i(.rms_lm(L, o + ys, w = po + L$w))
    expect_identical(r$null, "residual", info = sprintf("predictor %g, response %g", po, o))
    expect_true(r$exact)
    expect_identical(r$p, 2L)
  }
})


test_that("residuals that are not OLS residuals fail the test at every origin and in every unit", {
  L <- .rms_layer()
  y  <- L$z + 0.5 * L$w
  e0 <- as.numeric(stats::residuals(stats::lm(y ~ L$w)))
  e_not <- e0 + 0.05 * stats::sd(e0) * as.numeric(scale(L$z))
  keep <- rep(TRUE, L$n)
  ref <- residual_morans_i(.rms_fixed(L, y, e_not))
  expect_identical(ref$null, "randomisation")
  expect_false(ref$exact)
  for (p in c(0, -30, -60, -80, -90, 40)) {
    s <- 2^p
    f <- .rms_fixed(L, y * s, e_not * s)
    des <- .morans_ols_design(f, e_not * s, keep, .morans_response(f, keep))
    expect_false(des$is_ols, info = sprintf("scale 2^%d", p))
    r <- residual_morans_i(f)
    expect_identical(r$null, "randomisation")
    expect_false(r$exact)
    expect_identical(r$observed, ref$observed)
    expect_identical(r$z, ref$z)
  }
  # And at a predictor or response offset.
  for (o in c(1e8, 1e10)) {
    expect_identical(residual_morans_i(.rms_fixed(L, o + y, e_not))$null, "randomisation")
    expect_identical(residual_morans_i(.rms_fixed(L, y, e_not, w = o + L$w))$null,
                     "randomisation")
  }
  # null = "residual" still says the moments are an approximation there.
  lines <- capture_spatialkit_log(r <- residual_morans_i(.rms_fixed(L, 1e8 + y, e_not),
                                                         null = "residual"))
  expect_false(r$exact)
  expect_true(log_has(lines, "not the OLS residuals"))
})


test_that("a constant or collinear predictor is still left out of the rank", {
  L <- .rms_layer()
  y <- L$z + 0.5 * L$w
  d <- .rms_sf(L, y)
  d$k  <- 1e8
  d$w3 <- 2 * d$w + 1e8
  df <- sf::st_drop_geometry(d)
  for (pv in list(c("w", "k"), c("w", "w3"))) {
    fml <- stats::as.formula(paste("v ~", paste(pv, collapse = " + ")))
    f <- new_spatial_fit("rms_lm_fit", stats::lm(fml, data = df), fml, "v", pv, d)
    r <- residual_morans_i(f)
    expect_identical(r$p, 2L, info = paste(pv, collapse = " + "))
    expect_identical(r$null, "residual")
  }
})


test_that("the design and moments of a fit near the origin are those of the raw columns", {
  L <- .rms_layer()
  y <- L$z + 0.5 * L$w
  f <- .rms_lm(L, y, pv = c("w", "w2"))
  keep <- rep(TRUE, L$n)
  e <- residuals(f)
  des <- .morans_ols_design(f, e, keep, .morans_response(f, keep))
  X <- stats::model.matrix(~ ., data = sf::st_drop_geometry(f$data_sf)[c("w", "w2")])
  expect_identical(des$X, X)
  expect_null(des$centres)
  expect_true(des$is_ols)
})


test_that("residuals in small or large units give the I and z they give in unit ones", {
  L <- .rms_layer()
  y <- L$z + 0.5 * L$w
  ref <- residual_morans_i(.rms_lm(L, y))
  expect_identical(ref$null, "residual")
  for (p in c(-20, -30, -33, -40, -56, -200, -600, 200, 600)) {
    lines <- capture_spatialkit_log(r <- residual_morans_i(.rms_lm(L, y * 2^p)))
    expect_false(is.null(r), info = sprintf("scale 2^%d", p))
    expect_identical(unclass(r)[.rms_keys], unclass(ref)[.rms_keys],
                     info = sprintf("scale 2^%d", p))
    expect_false(log_has(lines, "zero residual variance"))
  }
  cm0 <- suppressWarnings(compare_models(list(m = .rms_lm(L, y))))
  cm1 <- suppressWarnings(compare_models(list(m = .rms_lm(L, y * 2^-30))))
  for (col in c("resid_morans_I", "resid_morans_z", "resid_morans_p", "resid_morans_null"))
    expect_identical(cm1[[col]], cm0[[col]], info = col)
})


test_that("the residuals of an exact fit are refused at every scale, with the reason", {
  L <- .rms_layer()
  withr::with_seed(3, cfs <- replicate(4, rnorm(3, 0, 5), simplify = FALSE))
  for (cf in cfs) {
    v <- cf[1] + cf[2] * L$w + cf[3] * L$w2
    for (p in c(0, 13, 20, 27, 30, -40)) {
      expect_warning(
        lines <- capture_spatialkit_log(r <- residual_morans_i(
          .rms_lm(L, v * 2^p, pv = c("w", "w2")))),
        "the rounding of an exact fit of the response on the predictors")
      expect_null(r)
      expect_false(log_has(lines, "zero residual variance"))
    }
  }
  # On predictors at an offset, as given and about their centres.
  wq <- (1e15 + L$w) - 1e15
  for (o in c(1e8, 4.5e6)) {
    v <- 3 + 4.317 * (o + wq)
    expect_warning(r <- residual_morans_i(.rms_lm(L, v, w = o + wq, centred = TRUE)),
                   "rounding of an exact fit")
    expect_null(r)
  }
  # A field on top is scored.
  r <- residual_morans_i(.rms_lm(L, 3 + 4.317 * (1e8 + wq) + 1e-3 * L$z, w = 1e8 + wq,
                                centred = TRUE))
  expect_false(is.null(r))
  expect_identical(r$null, "residual")
})


test_that("residuals that are not OLS residuals and vary by rounding alone are refused", {
  L <- .rms_layer()
  y <- 1e15 + 64 * round(L$z * 4)
  # Residuals of a backend whose fitted values differ from the response by
  # one unit in its last place: rounding of the response, not OLS residuals.
  withr::with_seed(5, res <- sample(c(-1, 0, 1), L$n, replace = TRUE) * 0.125)
  for (p in c(0, -60, 60)) {
    expect_warning(r <- residual_morans_i(.rms_fixed(L, y * 2^p, res * 2^p)),
                   "do not vary beyond the rounding of the response")
    expect_null(r)
  }
})


test_that("degenerate weights are refused with their own reason", {
  L <- .rms_layer()
  f <- .rms_lm(L, L$z + 0.5 * L$w)
  expect_warning(r <- residual_morans_i(f, weights = matrix(0, L$n, L$n)),
                 "degenerate weights")
  expect_null(r)
})


test_that("a GWR at a global window keeps the residual null with a predictor at an offset", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  L <- .rms_layer()
  wq <- (1e15 + L$w) - 1e15
  v  <- 2 * wq + L$z
  out <- lapply(c(0, 1e8), function(o) {
    d <- .rms_sf(L, v, w = o + wq)
    g <- suppressWarnings(suppressMessages(fit_gwr_model(
      d, "v", "w", adaptive = TRUE, bandwidth = L$n, kernel = "boxcar")))
    residual_morans_i(g)
  })
  expect_identical(out[[1]]$null, "residual")
  expect_equal(unclass(out[[2]])[.rms_keys], unclass(out[[1]])[.rms_keys], tolerance = 1e-8)
})


test_that("a backend that regressed the raw columns passes on 10,000 rows", {
  skip_if_not_installed("Matrix")
  # lm() leaves rounding on the raw columns that grows with the rows: 18 to
  # 30 eps times the size of the fit's terms on 10,000 rows at a response
  # offset of 1e6 to 1e10, above the bound of a regression's rounding.
  n <- 10000L
  withr::with_seed(n, {
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4); w <- rnorm(n)
    zz <- sin(x / 900) + cos(y / 700) + rnorm(n, sd = 0.5)
  })
  for (pv in list("w", character(0))) {
    fml <- if (length(pv)) v ~ w else v ~ 1
    fit_at <- function(o) {
      d <- sf::st_as_sf(data.frame(x = x, y = y, v = o + zz + 0.5 * w, w = w),
                        coords = c("x", "y"), crs = 32615)
      new_spatial_fit("rms_lm_fit", stats::lm(fml, data = sf::st_drop_geometry(d)),
                      fml, "v", pv, d)
    }
    ref <- residual_morans_i(fit_at(0))
    expect_identical(ref$null, "residual")
    expect_true(ref$exact)
    for (o in c(1e6, 1e8)) for (nl in c("auto", "residual")) {
      lab <- sprintf("%s, response at %g, null = %s", format(fml), o, nl)
      got <- .rms_quiet(residual_morans_i(fit_at(o), null = nl))
      expect_identical(got$null, "residual", info = lab)
      expect_true(got$exact, info = lab)
      expect_equal(got$z, ref$z, tolerance = 1e-8, info = lab)
      expect_false(log_has(attr(got, "log"), "not the OLS residuals"), info = lab)
    }
  }
})


test_that("an intercept-only lm at 1e15 keeps the residual null", {
  L <- .rms_layer()
  f <- function(v) {
    d <- .rms_sf(L, v)
    new_spatial_fit("rms_lm_fit", stats::lm(v ~ 1, data = sf::st_drop_geometry(d)),
                    v ~ 1, "v", character(0), d)
  }
  ref <- residual_morans_i(f(round(L$z * 2)))
  r <- .rms_quiet(residual_morans_i(f(1e15 + round(L$z * 2)), null = "residual"))
  expect_identical(r$null, "residual")
  expect_true(r$exact)
  expect_identical(r$p, 1L)
  expect_equal(r$expected, ref$expected, tolerance = 1e-12)
  expect_false(log_has(attr(r, "log"), "not the OLS residuals"))
})


test_that("the raw columns are read only where they keep the design's rank", {
  # lm() drops a predictor at 1e8 with a spread of 1, so its residuals are
  # those of the intercept alone, not the OLS residuals on the design.
  L <- .rms_layer()
  wq <- (1e15 + L$w) - 1e15
  r <- residual_morans_i(.rms_lm(L, 2 * wq + L$z, w = 1e8 + wq))
  expect_identical(r$null, "randomisation")
  expect_false(r$exact)
})


test_that("the residuals of an exact fit are refused under the randomisation null too", {
  L <- .rms_layer()
  for (p in c(0, 20)) {
    expect_warning(
      r <- residual_morans_i(.rms_lm(L, (3 + 2 * L$w - L$w2) * 2^p, pv = c("w", "w2")),
                             null = "randomisation"),
      "rounding of an exact fit")
    expect_null(r)
  }
})


test_that("an intercept-only fit of a response that varies by its rounding gives that reason", {
  L <- .rms_layer()
  # Steps of 0.125, one unit in the last place at 1e15.
  v <- 1e15 + 0.125 * round(L$z * 2)
  d <- .rms_sf(L, v)
  f <- new_spatial_fit("rms_lm_fit", stats::lm(v ~ 1, data = sf::st_drop_geometry(d)),
                       v ~ 1, "v", character(0), d)
  for (nl in c("auto", "residual", "randomisation")) {
    expect_warning(r <- residual_morans_i(f, null = nl),
                   "the response does not vary beyond its rounding")
    expect_null(r)
  }
})


test_that("residuals beyond the double range of the response's scale are scored", {
  L <- .rms_layer()
  ref <- residual_morans_i(.rms_fixed(L, L$z + 3, L$z))
  for (k in c(150, 200)) {
    r <- residual_morans_i(.rms_fixed(L, 10^-k * (L$z + 3), 10^k * L$z))
    expect_false(is.null(r), info = k)
    expect_equal(r$observed, ref$observed, tolerance = 1e-12, info = k)
    expect_identical(r$null, "randomisation", info = k)
  }
})


test_that(".pow2_below() gives the power of two at or below a number", {
  for (k in c(-1000, -60, -1, 0, 1, 52, 1000)) {
    expect_identical(.pow2_below(2^k), 2^k, info = k)
    expect_identical(.pow2_below(2^k * (1 - .Machine$double.eps / 2)), 2^(k - 1), info = k)
    expect_identical(.pow2_below(2^k * 1.5), 2^k, info = k)
  }
  expect_identical(.pow2_below(.Machine$double.xmax), 2^1023)
  for (a in list(0, -1, NA_real_, Inf, NaN)) expect_identical(.pow2_below(a), 1)
})
