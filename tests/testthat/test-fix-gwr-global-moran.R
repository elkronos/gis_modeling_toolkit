# fit_gwr_model() on a response whose centre is more than 1000 times its
# spread.  A GWR at a global window is OLS, and residual_morans_i() reads its
# residuals as OLS residuals.  Fitted on the raw response at 1e8 or 1e10, the
# local solve's residuals carried the rounding of the offset, which grows
# with the rows: on 6000 rows they missed the OLS check, and null = "auto"
# switched to randomisation in silence.  The response is now taken about its
# centre past the predictors' gate (.PRED_CENTRE_RATIO).

skip_if_not_installed("sp")

# 6000 points, a predictor w and a response with a spatial field, stored as
# at 1e10 so that every offset used holds the same values.
.gm_layer <- function(n = 6000L) {
  withr::with_seed(n, {
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4); w <- rnorm(n)
    zz <- sin(x / 900) + cos(y / 700) + rnorm(n, sd = 0.5)
  })
  v0 <- (1e10 + (zz + 0.5 * w)) - 1e10
  list(n = n, v0 = v0, x = x, y = y,
       mk = function(o) sf::st_as_sf(data.frame(x = x, y = y, v = o + v0, w = w),
                                     coords = c("x", "y"), crs = 32615))
}

.gm_rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300), na.rm = TRUE)

test_that("a GWR at a global window gives the residual null and the z of offset 0 at 1e8 and 1e10", {
  L <- .gm_layer()
  fit <- function(o) suppressWarnings(suppressMessages(
    fit_gwr_model(L$mk(o), "v", "w", adaptive = TRUE, bandwidth = L$n,
                  kernel = "boxcar")))
  ref <- fit(0)
  m0  <- suppressWarnings(residual_morans_i(ref))
  expect_identical(m0$null, "residual")
  expect_true(m0$exact)
  nd <- withr::with_seed(11, sf::st_as_sf(
    data.frame(x = runif(20, 1e3, 9e3), y = runif(20, 1e3, 9e3),
               w = rnorm(20)), coords = c("x", "y"), crs = 32615))
  p0 <- suppressWarnings(predict(ref, nd))
  ms <- list()
  for (o in c(1e8, 1e10)) {
    f <- fit(o)
    expect_identical(f$engine$centres$response, .spread_scaled(o + L$v0)$centre)
    m <- suppressWarnings(residual_morans_i(f))
    expect_identical(m$null, "residual")
    expect_true(m$exact)
    expect_lt(abs(m$z - m0$z), 1e-6)
    ms[[length(ms) + 1L]] <- m
    expect_lt(.gm_rel(f$engine$SDF$w, ref$engine$SDF$w), 1e-8)
    expect_lt(abs(f$info$AICc - ref$info$AICc) / abs(ref$info$AICc), 1e-8)
    expect_lt(abs(f$engine$GW.diagnostic$RSS.gw - ref$engine$GW.diagnostic$RSS.gw) /
                ref$engine$GW.diagnostic$RSS.gw, 1e-8)
    expect_lt(max(abs(residuals(f) - residuals(ref))), 1e-8 * max(abs(residuals(ref))))
    # Fitted values and predictions are in the response's units, to the
    # rounding a double holds at the offset.
    ulp <- 2 * .Machine$double.eps * o
    expect_lte(max(abs(fitted(f) - (o + fitted(ref)))), 2 * ulp)
    expect_lte(max(abs(suppressWarnings(predict(f, nd)) - (o + p0))), 2 * ulp)
  }
  for (k in c("observed", "expected", "sd", "z", "p_value"))
    expect_identical(ms[[2L]][[k]], ms[[1L]][[k]])
})

test_that("a global GWR of the intercept alone at -1e8 gives the residual null and the z of offset 0", {
  L <- .gm_layer()
  fit <- function(o) suppressWarnings(suppressMessages(
    fit_gwr_model(L$mk(o), "v", character(0), adaptive = TRUE,
                  bandwidth = L$n, kernel = "boxcar")))
  m0 <- suppressWarnings(residual_morans_i(fit(0)))
  expect_identical(m0$null, "residual")
  expect_true(m0$exact)
  f <- fit(-1e8)
  expect_identical(f$engine$centres$response, .spread_scaled(-1e8 + L$v0)$centre)
  expect_lt(f$engine$centres$response, 0)
  m <- suppressWarnings(residual_morans_i(f))
  expect_identical(m$null, "residual")
  expect_true(m$exact)
  expect_lt(abs(m$z - m0$z), 1e-6)
})

test_that("a response within 1000 times its spread of 0 is fitted as it is given", {
  L <- .gm_layer(300L)
  sc <- .spread_scaled(L$v0)
  for (r in c(-999, 0, 500, 999)) {
    # A shift that puts the centre at r times the spread.
    o <- r * sc$spread - sc$centre
    d <- L$mk(o)
    df <- sf::st_drop_geometry(d)
    s <- .spread_scaled(df$v)
    expect_lte(abs(s$centre), 1e3 * s$spread)
    dsn <- .gwr_design(v ~ w, df)
    cen <- .gwr_centre(dsn$X, dsn$y)
    expect_identical(cen$y_centre, 0)
    expect_identical(cen$y, dsn$y)
    f <- suppressWarnings(fit_gwr_model(d, "v", "w", bandwidth = 60))
    sv <- .gwr_local_solve(dsn$X, dsn$y, sf::st_coordinates(d), 60L,
                           "bisquare", TRUE, what = "fit")
    expect_identical(fitted(f), sv$fitted)
    expect_identical(residuals(f), df$v - sv$fitted)
    expect_identical(f$engine$centres$response, 0)
  }
  # At exactly 1000 times its spread the response is left as it is, as a
  # predictor is (.centre_columns()); one unit further it is centred.
  X1 <- matrix(1, 300L, 1L)
  for (c0 in c(1000, -1000)) {
    y <- c0 + rep(c(-1, 0, 1), 100L)
    expect_identical(.spread_scaled(y)$spread, 1)
    cen <- .gwr_centre(X1, y)
    expect_identical(cen$y_centre, 0)
    expect_identical(cen$y, y)
    y <- y + sign(c0)
    expect_identical(.gwr_centre(X1, y)$y_centre, c0 + sign(c0))
  }
})

test_that("a response past 1000 times its spread is fitted about its centre", {
  L <- .gm_layer(300L)
  ref <- suppressWarnings(fit_gwr_model(L$mk(0), "v", "w", bandwidth = 60))
  sc <- .spread_scaled(L$v0)
  for (o in c(1001 * sc$spread - sc$centre, 1e5, 1e8, -1e8)) {
    d <- L$mk(o)
    df <- sf::st_drop_geometry(d)
    s <- .spread_scaled(df$v)
    expect_gt(abs(s$centre), 1e3 * s$spread)
    f <- suppressWarnings(fit_gwr_model(d, "v", "w", bandwidth = 60))
    expect_identical(f$engine$centres$response, s$centre)
    expect_lt(max(abs(residuals(f) - residuals(ref))), 1e-9)
    expect_lt(abs(f$info$AICc - ref$info$AICc), 1e-8)
    expect_lt(.gm_rel(f$engine$SDF$w, ref$engine$SDF$w), 1e-8)
    # A GWR at a working bandwidth is not OLS, and keeps the randomisation
    # null.
    expect_identical(suppressWarnings(residual_morans_i(f))$null,
                     "randomisation")
  }
})

test_that("a GWR of the intercept alone past 1000 times its spread is fitted about its centre", {
  L <- .gm_layer(300L)
  fit <- function(o) suppressWarnings(
    fit_gwr_model(L$mk(o), "v", character(0), bandwidth = 60))
  ref <- fit(0)
  expect_identical(ref$engine$centres$response, 0)
  for (o in c(1e8, -1e8, 1e10)) {
    f <- fit(o)
    s <- .spread_scaled(o + L$v0)
    expect_identical(f$engine$centres$response, s$centre)
    expect_identical(sign(s$centre), sign(o))
    expect_lt(max(abs(residuals(f) - residuals(ref))), 1e-9)
    expect_lt(abs(f$info$AICc - ref$info$AICc), 1e-8)
  }
})

test_that("gwr_model_selection() on a response at 1e8 gives the table at 0", {
  L <- .gm_layer(300L)
  sel <- function(o) {
    d <- L$mk(o)
    d$q <- withr::with_seed(5, rnorm(nrow(d)))
    suppressWarnings(gwr_model_selection(d, "v", c("w", "q")))
  }
  ref <- sel(0)
  s <- sel(1e8)
  expect_identical(s$table$variables, ref$table$variables)
  expect_lt(max(abs(s$table$criterion - ref$table$criterion)), 1e-8)
  expect_identical(s$bandwidth, ref$bandwidth)
})
