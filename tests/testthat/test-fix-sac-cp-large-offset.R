# tests/testthat/test-fix-sac-cp-large-offset.R
# ---------------------------------------------------------------------------
# estimate_sac_range() and the test for `cp` in resolution_profile() read a
# response again, centred on its median and scaled, where its raw values
# show no spread, and read residuals against the rounding of the regression
# and of the raw response.  A response at a large offset (1e15 plus a spread
# of 1, eight units in the last place) is fitted and scored as it is at 0,
# with its regressions taken about its centre; a response of one value,
# values that differ only by their rounding and the residuals of an exact
# fit stay constant (BUG-36).
# ---------------------------------------------------------------------------

.lo_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# 800 points: an exponential field (range 600 m) plus 0.5 w plus noise.
.lo_layer <- function() {
  .lo_seed(5)
  n <- 800
  xy <- data.frame(x = stats::runif(n, 0, 4000), y = stats::runif(n, 0, 4000))
  D <- as.matrix(stats::dist(xy))
  f <- as.numeric(t(chol(exp(-D / 600) + diag(1e-8, n))) %*% stats::rnorm(n))
  xy$w  <- stats::rnorm(n)
  xy$z  <- f + 0.5 * xy$w + stats::rnorm(n, sd = 0.3)
  xy$w2 <- stats::rnorm(n)
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32615)
}

.lo_catch <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w, classes = vapply(w, function(c) class(c)[1L], ""))
}

.lo_profile <- function(L, v, pv = NULL)
  .lo_catch(resolution_profile(L, v, predictor_vars = pv, n_levels = 4, nstart = 2, seed = 1))

test_that("estimate_sac_range() fits 1e15 plus a spread of 1 as it fits the spread", {
  skip_if_not_installed("gstat")
  L <- .lo_layer()
  L$off  <- 1e15 + L$z
  L$grid <- L$off - 1e15              # the values as stored at 1e15, exact
  for (pv in list(NULL, "w")) {
    info <- if (is.null(pv)) "no predictors" else "predictor w"
    a <- suppressWarnings(estimate_sac_range(L, "off", predictor_vars = pv))
    expect_null(attr(a, "rejected_reason"), info = info)
    expect_true(is.finite(a), info = info)
    # The same answer as the values stored at 1e15 give at 0.
    g <- suppressWarnings(estimate_sac_range(L, "grid", predictor_vars = pv))
    expect_equal(as.numeric(a), as.numeric(g), tolerance = 1e-6, info = info)
    expect_equal(sac_nugget(a), sac_nugget(g), tolerance = 1e-6, info = info)
    expect_equal(attr(a, "structure_p"), attr(g, "structure_p"), tolerance = 1e-6, info = info)
    # And, to the eighth they are stored to, as the spread itself.
    z <- suppressWarnings(estimate_sac_range(L, "z", predictor_vars = pv))
    expect_equal(as.numeric(a), as.numeric(z), tolerance = 1e-3, info = info)
    expect_equal(sac_nugget(a), sac_nugget(z), tolerance = if (is.null(pv)) 0.01 else 0.05,
                 info = info)
  }
})

test_that("REML detrending of 1e15 plus a spread of 1 is that of the spread", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  L <- .lo_layer()
  L$off  <- 1e15 + L$z
  L$grid <- L$off - 1e15
  a <- suppressWarnings(estimate_sac_range(L, "off", predictor_vars = "w", detrend = "reml"))
  g <- suppressWarnings(estimate_sac_range(L, "grid", predictor_vars = "w", detrend = "reml"))
  expect_null(attr(a, "rejected_reason"))
  expect_equal(as.numeric(a), as.numeric(g), tolerance = 1e-6)
})

test_that("the profile scores cp on 1e15 plus a spread of 1, with and without predictors", {
  skip_if_not_installed("gstat")
  L <- .lo_layer()
  L$off <- 1e15 + L$z
  for (pv in list(NULL, "w")) {
    info <- if (is.null(pv)) "no predictors" else "predictor w"
    r0 <- .lo_profile(L, "z", pv)
    r  <- .lo_profile(L, "off", pv)
    p  <- r$value
    expect_false(attr(p, "one_cell")$constant, info = info)
    expect_false("spatialkit_constant_response" %in% r$classes, info = info)
    expect_true(all(is.finite(p$cp)), info = info)
    expect_identical(p$levels, r0$value$levels, info = info)
    # The OLS fit is taken about the response's centre: on the raw values
    # its rounding at 1e15 gave a sum of squares about 150 times the
    # spread's.
    expect_equal(p$rss, r0$value$rss, tolerance = 0.01, info = info)
    expect_identical(select_resolution(p, "cp")$best, select_resolution(r0$value, "cp")$best,
                     info = info)
  }
})

test_that("a trend the predictors explain at a large offset leaves the spread to fit", {
  # The raw values show their spread through the trend, so the response is
  # not centred for its own sake; the residuals of the fit on the raw values
  # hold the rounding of the offset, and the fit is taken again about the
  # centre (BUG-36 review).
  skip_if_not_installed("gstat")
  L <- .lo_layer()
  cases <- list(
    "1e15 + 1e5 w2 + z"          = list(off = 1e15, y = 1e15 + 1e5 * L$w2 + L$z),
    "1e8 + 1e-3 w2 + 1e-7 z"     = list(off = 1e8, y = 1e8 + 1e-3 * L$w2 + 1e-7 * L$z))
  for (nm in names(cases)) {
    d <- L
    d$v <- cases[[nm]]$y
    d$g <- d$v - cases[[nm]]$off        # the values as stored, at 0, exact
    expect_identical(.regression_centre(d$v), 0, info = nm)
    a <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = c("w", "w2")))
    g <- suppressWarnings(estimate_sac_range(d, "g", predictor_vars = c("w", "w2")))
    expect_null(attr(a, "rejected_reason"), info = nm)
    expect_true(is.finite(a), info = nm)
    expect_equal(as.numeric(a), as.numeric(g), tolerance = 1e-6, info = nm)
    expect_equal(sac_nugget(a), sac_nugget(g), tolerance = 1e-6, info = nm)
    r  <- .lo_profile(d, "v", c("w", "w2"))
    r0 <- .lo_profile(d, "g", c("w", "w2"))
    p  <- r$value
    expect_false(attr(p, "one_cell")$constant, info = nm)
    expect_false("spatialkit_constant_response" %in% r$classes, info = nm)
    expect_true(all(is.finite(p$cp)), info = nm)
    expect_identical(p$levels, r0$value$levels, info = nm)
    expect_equal(p$rss, r0$value$rss, tolerance = 1e-6, info = nm)
    expect_equal(p$cp, r0$value$cp, tolerance = 1e-6, info = nm)
    expect_identical(select_resolution(p, "cp")$best, select_resolution(r0$value, "cp")$best,
                     info = nm)
  }
})

test_that("REML detrending of a trend at a large offset is that of the values at 0", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  L <- .lo_layer()
  L$v <- 1e15 + 1e5 * L$w2 + L$z
  L$g <- L$v - 1e15
  a <- suppressWarnings(estimate_sac_range(L, "v", predictor_vars = c("w", "w2"),
                                           detrend = "reml"))
  g <- suppressWarnings(estimate_sac_range(L, "g", predictor_vars = c("w", "w2"),
                                           detrend = "reml"))
  expect_null(attr(a, "rejected_reason"))
  # The fit is taken about the response's median, not about 1e15, and the
  # REML optimiser's tolerance follows the intercept: 6e-6 apart.
  expect_equal(as.numeric(a), as.numeric(g), tolerance = 1e-4)
})

test_that("a response of one value, rounding and the residuals of exact fits stay constant", {
  skip_if_not_installed("gstat")
  L <- .lo_layer()
  .lo_seed(7)
  cases <- list(
    "one value 7"          = list(y = rep(7, nrow(L)), pv = NULL),
    "one value 1e15"       = list(y = rep(1e15, nrow(L)), pv = NULL),
    "0.1 beside 0.3 / 3"   = list(y = sample(c(0.1, 0.3 / 3), nrow(L), TRUE), pv = NULL),
    "1e15 + z / 4"         = list(y = 1e15 + L$z / 4, pv = NULL),
    "exact 7 + 2 w"        = list(y = 7 + 2 * L$w, pv = "w"),
    "exact 1e15 + 2 w"     = list(y = 1e15 + 2 * L$w, pv = "w"),
    "exact 1e6 + w - 3 w2" = list(y = 1e6 + L$w - 3 * L$w2, pv = c("w", "w2")),
    # Exact fits whose residuals on the raw values come out above three
    # units in the last place, fitted again about the centre.
    "exact 5 + 2 w"        = list(y = 5 + 2 * L$w, pv = "w"),
    "exact 100 + 3 w"      = list(y = 100 + 3 * L$w, pv = "w"),
    # Exact fits about a centre near 0, whose residuals in the profile come
    # out above three units in the last place: held by the bound of the
    # regression's rounding.
    "exact w"              = list(y = L$w, pv = "w"),
    "exact 1e6 w"          = list(y = 1e6 * L$w, pv = "w"),
    # A trend the predictors explain at a large offset, fitted again about
    # the response's centre.
    "exact 1e15 + 1e5 w2 + 2 w" = list(y = 1e15 + 1e5 * L$w2 + 2 * L$w, pv = c("w", "w2")))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    d <- L; d$v <- cs$y
    s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = cs$pv))
    expect_identical(attr(s, "rejected_reason"),
                     if (is.null(cs$pv)) "the response is constant"
                     else paste0("the residuals on predictor_vars are constant: the predictors ",
                                 "explain the response exactly"), info = nm)
    r <- .lo_catch(resolution_profile(d, "v", predictor_vars = cs$pv, n_levels = 3, nstart = 1,
                                      seed = 1))
    expect_true(attr(r$value, "one_cell")$constant, info = nm)
    expect_true(all(is.na(r$value$cp)), info = nm)
    expect_true("spatialkit_constant_response" %in% r$classes, info = nm)
  }
})

test_that("a response whose raw values show their spread is computed as before", {
  # The centre is 0 there, and both readings are those of the raw values.
  y <- c(1, 2, 3, 4) + 1e9
  expect_identical(.regression_centre(y), 0)
  # A standard deviation just above 1e-12 of the size.
  .lo_seed(3)
  expect_identical(.regression_centre(1e8 + 1e-3 * stats::rnorm(50)), 0)
  # Constant columns are not centred, bounded or not.
  expect_identical(.regression_centre(rep(1e15, 5)), 0)
  expect_identical(.regression_centre(rep(1e15, 5), bounded = FALSE), 0)
  expect_identical(.regression_centre(c(0.1, 0.3 / 3, 0.1, 0.1)), 0)
  # The centre a fit is taken about again: that of .spread_scaled().
  expect_identical(.residual_centre(1e15 + c(0, 0.125, 0.25, 1, 2)), 1e15 + 0.25)
  expect_identical(.residual_centre(c(NA, 3, 5, Inf, 7)), 5)
  expect_identical(.residual_centre(rep(0, 4)), 0)
  expect_identical(.residual_centre(7), 0)
  expect_identical(.residual_centre(1e160 * (1 + 1e-15 * c(-2, -1, 0, 1, 2))), 0)
  expect_identical(.residual_centre(1e160 * (1 + 1e-15 * c(-2, -1, 0, 1, 2)), bounded = FALSE),
                   1e160)
  expect_identical(.regression_centre(1e15 + c(0, 0.125, 0.25, 1, 2)), 1e15 + 0.25)
  # Beyond the pooling's arithmetic the profile's centre stays 0, and
  # estimate_sac_range()'s does not.
  big <- 1e160 * (1 + 1e-15 * c(-2, -1, 0, 1, 2))
  expect_identical(.regression_centre(big), 0)
  expect_identical(.regression_centre(big, bounded = FALSE), 1e160)
  expect_identical(.ulp_of(1e15), 0.125)
  expect_identical(.ulp_of(1), .Machine$double.eps)
  expect_identical(.ulp_of(0), 0)
})
