# estimate_sac_range(detrend = "reml") on a predictor or a response at an
# offset, and estimate_sac_range() on a lattice in another unit or origin.
# nlme::gls() failed on a predictor or a response at 1e8 plus a spread of 1
# or 2, and the estimate fell back to OLS with the reason "did not
# converge", so detrend_method, the range and the nugget moved with the
# origin of the column.  And .reml_trend() and the n_max subsample sorted
# the rows by the exact first coordinate before their draws, so on a
# lattice whose first coordinates differ in their last bits after another
# unit or origin the rows drawn, the range and the nugget moved with them.

# A smooth field, a predictor w and noise on 300 points; w as stored at
# 1e15, so that 1e8 + wq and 1e15 + wq hold the same spread.
.ro_layer <- function(n = 300, seed = 5) {
  withr::with_seed(seed, {
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4)
    w <- rnorm(n)
    z <- 2 * w + sin(x / 300) + cos(y / 300) + rnorm(n, sd = 0.5)
  })
  wq <- (1e15 + w) - 1e15
  list(xy = cbind(x, y), wq = wq,
       sf = sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"),
                         crs = 32615))
}

# Warnings caught and returned with the value.
.ro_run <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

test_that("REML detrending on a predictor at a large offset gives the fit of the predictor at 0", {
  skip_if_not_installed("nlme")
  lay <- .ro_layer()
  d <- lay$sf; d$p <- lay$wq
  ref <- .ro_run(estimate_sac_range(d, "z", predictor_vars = "p", detrend = "reml",
                                    seed = 1))
  expect_identical(attr(ref$value, "detrend_method"), "reml")
  expect_false(any(grepl("REML detrending", ref$warnings)))
  for (off in c(1e5, 1e8, 1e15)) {
    d$p <- off + lay$wq
    got <- .ro_run(estimate_sac_range(d, "z", predictor_vars = "p", detrend = "reml",
                                      seed = 1))
    info <- format(off)
    expect_identical(attr(got$value, "detrend_method"), "reml", info = info)
    expect_false(any(grepl("did not converge|falling back", got$warnings)), info = info)
    expect_identical(is.na(as.numeric(got$value)), is.na(as.numeric(ref$value)), info = info)
    expect_equal(as.numeric(got$value), as.numeric(ref$value), tolerance = 1e-6,
                 info = info)
    expect_equal(attr(got$value, "nugget"), attr(ref$value, "nugget"), tolerance = 1e-6,
                 info = info)
    expect_equal(attr(got$value, "reml")$nugget_prop, attr(ref$value, "reml")$nugget_prop,
                 tolerance = 1e-6, info = info)
    expect_identical(attr(got$value, "rejected_reason"), attr(ref$value, "rejected_reason"),
                     info = info)
  }
})

test_that("the REML fit is given predictors at an offset about their centres, and the rest as they are", {
  lay <- .ro_layer()
  d <- lay$sf
  seen <- NULL
  local_mocked_bindings(.reml_trend = function(mf, ...) {
    seen <<- mf
    NULL
  })
  # Below the gate (a year-like column, and w itself) the model frame holds
  # the values as given, to the last bit.
  d$p <- 2020 + round(lay$wq * 5)
  d$q <- lay$wq
  suppressWarnings(estimate_sac_range(d, "z", predictor_vars = c("p", "q"),
                                      detrend = "reml"))
  ord <- order(seen$.sac_x)
  expect_identical(seen$p[ord], d$p[order(lay$xy[, 1])])
  expect_identical(seen$q[ord], d$q[order(lay$xy[, 1])])
  # At an offset the column goes less its centre, and keeps its spread.
  d$p <- 1e8 + lay$wq
  suppressWarnings(estimate_sac_range(d, "z", predictor_vars = c("p", "q"),
                                      detrend = "reml"))
  expect_lt(max(abs(seen$p)), 10)
  expect_equal(stats::sd(seen$p), stats::sd(lay$wq), tolerance = 1e-6)
  expect_identical(seen$q[order(seen$.sac_x)], d$q[order(lay$xy[, 1])])
})

test_that("an exact REML fit on a northing fitted about its centre is still constant", {
  # The northing 4.5e6 plus a span of 4,000 is centred before the REML fit,
  # and a response computed from it carries the rounding of 4.317 times
  # 4.5e6: the test for an exact fit reads the terms of the predictor as
  # given (.rounding_design()), as the OLS path does.  The REML fit is
  # replaced by the exact coefficients on the frame it is given.
  lay <- .ro_layer()
  N <- 4.5e6 + lay$xy[, 2] * 0.4
  v <- 1.7 - 4.317 * N + 4.317 * mean(N)
  d <- lay$sf; d$N <- N; d$v <- v
  local_mocked_bindings(.reml_trend = function(mf, fml, ...) {
    list(beta = stats::coef(stats::lm(fml, data = mf)), range = 1000,
         nugget_prop = 0.1, sigma2 = 1, n_used = nrow(mf), subsampled = FALSE,
         n_colocated = 0L, pair_floor = 10, logLik = 0, logLik_null = -1,
         n_fits = 1L)
  })
  s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = "N",
                                           detrend = "reml"))
  expect_identical(attr(s, "rejected_reason"),
                   "the residuals on predictor_vars are constant: the predictors explain the response exactly")
  # A factor ahead of the northing: the centres go to the columns of the
  # design through attr(X, "assign"), the factor's two columns at 0.
  d$g <- factor(withr::with_seed(11, sample(c("a", "b", "c"), nrow(d), TRUE)))
  d$v <- v + 2 * as.integer(d$g)
  s <- .ro_run(estimate_sac_range(d, "v", predictor_vars = c("g", "N"),
                                  detrend = "reml"))
  expect_false(any(grepl("longer object length", s$warnings)))
  expect_identical(attr(s$value, "rejected_reason"),
                   "the residuals on predictor_vars are constant: the predictors explain the response exactly")
  # The response at an offset as well.  An exact fit is still constant.
  d$g <- NULL
  d$v <- 1e8 + v
  s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = "N",
                                           detrend = "reml"))
  expect_identical(attr(s, "rejected_reason"),
                   "the residuals on predictor_vars are constant: the predictors explain the response exactly")
  # And a response at an offset with a spread of its own is taken about its
  # centre for the REML fit, with the centre added back to the intercept
  # the test for a constant reads: the coefficients it is given give the
  # response it is given, less the residuals.
  real <- .residuals_rounding
  args <- list()
  local_mocked_bindings(.residuals_rounding = function(e, y, yd, X = NULL, b = NULL) {
    if (!is.null(X)) args[[length(args) + 1L]] <<- list(e = e, yd = yd, X = X, b = b)
    real(e, y, yd, X, b)
  })
  d$v <- 1e8 + v + withr::with_seed(12, rnorm(nrow(d)))
  s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = "N",
                                           detrend = "reml"))
  expect_false(identical(attr(s, "rejected_reason"),
                         "the residuals on predictor_vars are constant: the predictors explain the response exactly"))
  expect_gte(length(args), 1L)
  for (a in args) {
    fitted <- drop(as.matrix(a$X) %*% as.numeric(a$b))
    expect_lt(max(abs(a$yd - fitted - a$e)), 1e-4)
  }
})

test_that("REML detrending of a response at a large offset gives the fit of the response at 0", {
  # gls() failed on a response at 1e8 plus a spread of about 2 as it did on
  # such a predictor, and moved in the fifth digit at 1e7.  The response is
  # stored at 1e12, so that every offset below holds the same values.
  skip_if_not_installed("nlme")
  lay <- .ro_layer()
  w <- lay$wq
  zq <- (1e12 + lay$sf$z) - 1e12
  d <- lay$sf; d$z <- zq; d$p <- w
  ref <- .ro_run(estimate_sac_range(d, "z", predictor_vars = "p", detrend = "reml",
                                    seed = 1))
  expect_identical(attr(ref$value, "detrend_method"), "reml")
  expect_true(is.finite(as.numeric(ref$value)))
  for (off in c(1e7, 1e8, 1e12)) {
    d$z <- off + zq
    got <- .ro_run(estimate_sac_range(d, "z", predictor_vars = "p", detrend = "reml",
                                      seed = 1))
    info <- format(off)
    expect_identical(attr(got$value, "detrend_method"), "reml", info = info)
    expect_false(any(grepl("did not converge|falling back", got$warnings)), info = info)
    expect_equal(as.numeric(got$value), as.numeric(ref$value), tolerance = 1e-6,
                 info = info)
    expect_equal(attr(got$value, "nugget"), attr(ref$value, "nugget"), tolerance = 1e-6,
                 info = info)
    expect_equal(attr(got$value, "reml")$nugget_prop, attr(ref$value, "reml")$nugget_prop,
                 tolerance = 1e-6, info = info)
  }
  # The response and the predictor both at 1e8.
  d$z <- 1e8 + zq; d$p <- 1e8 + ((1e15 + w) - 1e15)
  got <- .ro_run(estimate_sac_range(d, "z", predictor_vars = "p", detrend = "reml",
                                    seed = 1))
  expect_identical(attr(got$value, "detrend_method"), "reml")
  expect_false(any(grepl("did not converge|falling back", got$warnings)))
  # Below the bound the response reaches the REML fit as given.
  seen <- NULL
  local_mocked_bindings(.reml_trend = function(mf, ...) {
    seen <<- mf
    NULL
  })
  d$z <- 1e3 + zq; d$p <- w
  suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p", detrend = "reml"))
  expect_identical(sort(seen$z), sort(d$z))
  d$z <- 1e8 + zq
  suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p", detrend = "reml"))
  expect_lt(max(abs(seen$z)), 100)
  expect_equal(stats::sd(seen$z), stats::sd(zq), tolerance = 1e-9)
})

# A 20 x 15 lattice at 250 m in EPSG:26915, and the same with the false
# easting and northing moved by 0.1 m: there the first coordinates of a
# lattice column differ in their last bits.
.ro_lattice <- function(jitter = 0) {
  g <- expand.grid(i = 0:19, j = 0:14)
  jx <- if (jitter > 0) withr::with_seed(99, runif(nrow(g), -jitter, jitter)) else 0
  jy <- if (jitter > 0) withr::with_seed(98, runif(nrow(g), -jitter, jitter)) else 0
  x <- 4.5e5 + g$i * 250 + jx; y <- 4.55e6 + g$j * 250 + jy
  withr::with_seed(1, {
    w <- (x - min(x)) / 1000 + rnorm(length(x))
    z <- sin((x - min(x)) / 1500) + cos((y - min(y)) / 1250) + 0.5 * w +
      rnorm(length(x), 0, 1)
  })
  m <- sf::st_as_sf(data.frame(x = x, y = y, z = z, w = w), coords = c("x", "y"),
                    crs = 26915)
  m01 <- sf::st_transform(m, sf::st_crs(paste(
    "+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1 +y_0=0.1",
    "+datum=NAD83 +units=m +no_defs")))
  list(m = m, m01 = m01)
}

test_that("the REML subsample of a lattice is the same rows after a move of the origin", {
  skip_if_not_installed("nlme")
  L <- .ro_lattice()
  xy <- lapply(L, function(l) sf::st_coordinates(l))
  # The exact order of the moved lattice is not that of the lattice, so the
  # test reaches the case.
  expect_false(identical(order(xy$m[, 1], xy$m[, 2]), order(xy$m01[, 1], xy$m01[, 2])))
  est <- lapply(L, function(l) suppressWarnings(
    estimate_sac_range(l, "z", predictor_vars = "w", detrend = "reml", reml_max_n = 100L)))
  expect_identical(attr(est$m, "detrend_method"), "reml")
  expect_true(isTRUE(attr(est$m, "reml")$subsampled))
  expect_equal(attr(est$m01, "nugget"), attr(est$m, "nugget"), tolerance = 1e-6)
  # The REML nugget proportion ends near 0 here, where the optimiser's end
  # point moves with rounding; the nugget read off the residuals does not.
  expect_identical(attr(est$m01, "detrend_method"), "reml")
  expect_lt(abs(attr(est$m01, "reml")$nugget_prop - attr(est$m, "reml")$nugget_prop),
            1e-4)
  expect_identical(is.na(as.numeric(est$m01)), is.na(as.numeric(est$m)))
  expect_equal(as.numeric(est$m01), as.numeric(est$m), tolerance = 1e-6)
})

test_that(".reml_trend() orders a lattice's rows by columns, and other layers as before", {
  # The rows a mocked gls() is given are those of the subsample; on the
  # lattice and the moved lattice they are the same points, and on a
  # jittered lattice the order is the exact order of the coordinates.
  skip_if_not_installed("nlme")
  rows_of <- function(l, max_n) {
    xy <- sf::st_coordinates(l)
    mf <- data.frame(z = l$z, w = l$w, .sac_x = xy[, 1], .sac_y = xy[, 2])
    got <- NULL
    local_mocked_bindings(gls = function(model, data, ...) {
      got <<- data
      stop("not fitted")
    }, .package = "nlme")
    .reml_trend(mf, z ~ w, extent = 5000, max_n = max_n, seed = 7L)
    got
  }
  L <- .ro_lattice()
  a <- rows_of(L$m, 100L); b <- rows_of(L$m01, 100L)
  expect_identical(a$z, b$z)
  expect_identical(a$w, b$w)
  J <- .ro_lattice(jitter = 60)
  xy <- sf::st_coordinates(J$m)
  full <- rows_of(J$m, 1000L)
  o <- order(xy[, 1], xy[, 2])
  expect_identical(full$z, J$m$z[o])
})

# The lattice of the item's reproduction: 250 m spacing in EPSG:26915, a
# field z2 of short range and a predictor w.
.ro_grid <- function(nx, ny, seed = 1, jitter = 0) {
  g <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  if (jitter > 0) {
    x <- x + withr::with_seed(99, runif(length(x), -jitter, jitter))
    y <- y + withr::with_seed(98, runif(length(y), -jitter, jitter))
  }
  w <- withr::with_seed(seed + 1, (x - min(x)) / 1000 + rnorm(length(x)))
  z2 <- withr::with_seed(seed + 2, {
    u <- sin((x - min(x)) / 300 + runif(1)) * cos((y - min(y)) / 260)
    u + rnorm(length(x), 0, 0.4)
  })
  sf::st_as_sf(data.frame(x = x, y = y, z2 = z2, w = w), coords = c("x", "y"),
               crs = 26915)
}
.ro_moved <- function(m) sf::st_transform(m, sf::st_crs(paste(
  "+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1 +y_0=0.1",
  "+datum=NAD83 +units=m +no_defs")))

test_that("an identified REML range on a lattice is the same after a move of the origin", {
  # Test 4's range is refused as NA in both frames; here it is a number.
  skip_if_not_installed("nlme")
  m <- .ro_grid(40, 30); m$zz <- m$z2 + 0.5 * m$w
  L <- list(m = m, m01 = .ro_moved(m))
  est <- lapply(L, function(l) suppressWarnings(
    estimate_sac_range(l, "zz", predictor_vars = "w", detrend = "reml",
                       reml_max_n = 300L)))
  expect_identical(attr(est$m, "detrend_method"), "reml")
  expect_identical(attr(est$m01, "detrend_method"), "reml")
  expect_true(isTRUE(attr(est$m, "reml")$subsampled))
  expect_true(is.finite(as.numeric(est$m)))
  expect_equal(as.numeric(est$m01), as.numeric(est$m), tolerance = 1e-6)
  expect_equal(attr(est$m01, "reml")$sigma2, attr(est$m, "reml")$sigma2, tolerance = 1e-6)
  expect_equal(attr(est$m01, "nugget"), attr(est$m, "nugget"), tolerance = 1e-6)
})

test_that("the one-row-a-location draw of a lattice with repeat visits takes the same rows after a move of the origin", {
  # Every location visited twice and reml_max_n above the number of
  # locations: no subsample, only the draw of one row a location, which
  # reads the same order.
  skip_if_not_installed("nlme")
  l2 <- .ro_grid(20, 15)
  l2b <- l2; l2b$z2 <- l2b$z2 + withr::with_seed(4, rnorm(nrow(l2b), 0, 0.3))
  ld <- rbind(l2, l2b)
  L <- list(m = ld, m01 = .ro_moved(ld))
  rows_of <- function(l) {
    xy <- sf::st_coordinates(l)
    mf <- data.frame(z2 = l$z2, w = l$w, .sac_x = xy[, 1], .sac_y = xy[, 2])
    got <- NULL
    local_mocked_bindings(gls = function(model, data, ...) {
      got <<- data
      stop("not fitted")
    }, .package = "nlme")
    .reml_trend(mf, z2 ~ w, extent = 5000, max_n = 1000L, seed = 7L)
    got
  }
  a <- rows_of(L$m); b <- rows_of(L$m01)
  expect_identical(nrow(a), 300L)
  expect_identical(a$z2, b$z2)
  expect_identical(a$w, b$w)
})

test_that("the n_max subsample of a lattice is the same rows after a move of the origin", {
  # The variogram's own subsample, drawn by position from the rows in the
  # order of their coordinates.
  m <- .ro_grid(40, 30)
  L <- list(m = m, m01 = .ro_moved(m))
  xy <- lapply(L, sf::st_coordinates)
  expect_false(identical(order(xy$m[, 1], xy$m[, 2]), order(xy$m01[, 1], xy$m01[, 2])))
  est <- lapply(L, function(l) suppressWarnings(estimate_sac_range(l, "z2", n_max = 300L)))
  expect_true(is.finite(as.numeric(est$m)))
  expect_equal(as.numeric(est$m01), as.numeric(est$m), tolerance = 1e-6)
  expect_equal(attr(est$m01, "nugget"), attr(est$m, "nugget"), tolerance = 1e-6)
  # The rows go to the projection in the same order in both frames, and a
  # jittered lattice, with no two first coordinates that close, in the
  # order of its coordinates.
  rows_in <- function(l) {
    seen <- NULL
    local_mocked_bindings(ensure_projected = function(x, ...) {
      if (is.null(seen)) seen <<- x
      x
    })
    suppressWarnings(estimate_sac_range(l, "z2", n_max = 300L))
    seen
  }
  expect_identical(rows_in(L$m)$z2, rows_in(L$m01)$z2)
  j <- .ro_grid(40, 30, jitter = 60)
  xj <- sf::st_coordinates(j)
  expect_identical(rows_in(j[nrow(j):1, ])$z2, j$z2[order(xj[, 1], xj[, 2])])
})
