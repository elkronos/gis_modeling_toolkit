# residual_morans_i() on a GWR at a working bandwidth whose response sits at
# a large offset.  fit_gwr_model() regresses such a response about its
# centre and residuals.gwr_fit() returns the residuals of the centred
# response, which carry no rounding of the offset.  The second OLS test of
# .morans_ols_design() still allowed them the rounding of a regression at
# the offset (10 eps times the offset), so at 1e15 a GWR whose residuals
# differ from the OLS residuals by 0.9 response units passed as OLS and
# null = "auto" took the residual null, exact TRUE.  The bound now reads the
# offset the fit regressed at.

skip_if_not_installed("sp")

# 300 points, a predictor p and a response z = k (2 w + a field + noise),
# stored as at `at` so that every offset used holds the same values.
.mgc_layer <- function(k = 1, at = 1e15) {
  withr::with_seed(3, {
    n <- 300
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4); w <- rnorm(n)
    z <- 2 * w + sin(x / 2000) + cos(y / 1500) + rnorm(n, sd = 0.5)
  })
  wq <- (1e15 + w) - 1e15
  zq <- (at + k * z) - at
  list(n = n, zq = zq,
       mk = function(o) sf::st_as_sf(data.frame(x = x, y = y, z = o + zq, p = wq),
                                     coords = c("x", "y"), crs = 32615))
}

.mgc_fit <- function(L, o, ...) suppressWarnings(suppressMessages(
  fit_gwr_model(L$mk(o), "z", "p", ...)))

.mgc_rel <- function(a, b) abs(a - b) / abs(b)

# The GWR at offset `o` gives the null, moments and p-value of offset 0,
# logs no warning, and compare_models() gives the same resid_morans columns.
.mgc_check <- function(L, offsets, bw) {
  f0 <- .mgc_fit(L, 0, bandwidth = bw)
  m0 <- residual_morans_i(f0)
  expect_identical(m0$null, "randomisation")
  expect_false(m0$exact)
  # compare_models() logs its caution on the residuals; the lines are kept
  # out of the test output.
  cmp <- function(g) {
    capture_spatialkit_log(out <- suppressWarnings(compare_models(list(g = g))))
    out
  }
  c0 <- cmp(f0)
  for (o in offsets) {
    what <- sprintf("bandwidth %g, offset %g", bw, o)
    f <- .mgc_fit(L, o, bandwidth = bw)
    expect_false(f$engine$centres$response == 0, info = what)
    lines <- capture_spatialkit_log(m <- residual_morans_i(f), level = logger::WARN)
    expect_length(lines, 0L)
    expect_false(is.null(m), info = what)
    expect_identical(m$null, "randomisation", info = what)
    expect_false(m$exact, info = what)
    expect_identical(m$df, m0$df, info = what)
    for (k in c("observed", "expected", "sd", "z", "p_value"))
      expect_lt(.mgc_rel(m[[k]], m0[[k]]), 1e-8, label = paste(what, k))
    cm <- cmp(f)
    expect_identical(cm$resid_morans_null, c0$resid_morans_null, info = what)
    for (col in c("resid_morans_I", "resid_morans_z", "resid_morans_p"))
      expect_lt(.mgc_rel(cm[[col]], c0[[col]]), 1e-8, label = paste(what, col))
  }
}

test_that("a GWR at a working bandwidth with the response at 1e13 to 1e15 keeps the null of offset 0", {
  L <- .mgc_layer()
  for (bw in c(30, 60)) .mgc_check(L, c(1e13, 1e14, 1e15, -1e15), bw)
})

test_that("a GWR at a working bandwidth with the response at 4e15 keeps the null of offset 0", {
  # Residuals four times wider, so that they vary beyond three units in the
  # last place of a response at 4e15 (0.5), which the no-spread test reads.
  L <- .mgc_layer(k = 4, at = 4e15)
  for (bw in c(30, 60)) .mgc_check(L, c(1e15, 4e15), bw)
})

test_that("a GWR at a global window keeps the residual null at 0, 1e8 and 1e15", {
  L <- .mgc_layer()
  for (o in c(0, 1e8, 1e15)) {
    f <- .mgc_fit(L, o, adaptive = TRUE, bandwidth = L$n, kernel = "boxcar")
    m <- suppressWarnings(residual_morans_i(f))
    expect_identical(m$null, "residual", info = sprintf("offset %g", o))
    expect_true(m$exact, info = sprintf("offset %g", o))
  }
})

test_that("only a GWR that keeps its residuals about a recorded centre reads that centre", {
  L <- .mgc_layer()
  f <- .mgc_fit(L, 1e15, bandwidth = 60)
  expect_identical(.morans_fit_centre(f), f$engine$centres$response)
  # Centre 0: the response was regressed as it is.
  expect_identical(.morans_fit_centre(.mgc_fit(L, 0, bandwidth = 60)), 0)
  # No residual column: residuals.gwr_fit() returns the raw values less
  # the fitted values, which carry the rounding of the offset.
  g <- f
  g$engine$SDF@data$residual <- NULL
  expect_identical(.morans_fit_centre(g), 0)
  # A residual column not of the data's length: residuals.gwr_fit() again
  # returns the raw values less the fitted values, so the centre is not read.
  g2 <- f
  g2$engine$SDF <- g2$engine$SDF[1:10, ]
  expect_identical(.morans_fit_centre(g2), 0)
  expect_identical(residuals(g2),
                   sf::st_drop_geometry(g2$data_sf)[["z"]] - fitted(g2))
  # A fit of another class with the same field is read at the offset.
  h <- f
  class(h) <- setdiff(class(h), "gwr_fit")
  expect_identical(.morans_fit_centre(h), 0)
  g$engine$centres$response <- NA_real_
  expect_identical(.morans_fit_centre(g), 0)
})
