# estimate_sac_range(detrend = "reml") where the REML range is not
# identified.  On 100 to 150 points of a lattice the REML range runs to 1e8
# m and more with a nugget proportion near 1e-6, the likelihood is flat, and
# whether nlme::gls() returned a fit or failed on the same rows was set by
# the last bits of the coordinates.  So detrend_method, the nugget and the
# fallback warning moved with the unit and the origin, and with the unit and
# the last bits of the response and the predictors.  .reml_trend() now hands
# gls() the coordinates of .kmeans_frame() (centred on the middle of the
# bounding box, divided by its diagonal, rounded to .KMEANS_STEP), and the
# response and each numeric predictor in the same frame of their range; it
# returns the range, pair_floor, the coefficients, the variance and the
# log-likelihoods in the caller's units.  The extent of the starting ranges
# has no floor of 1 in the caller's unit.

# A lattice at 250 m in EPSG:26915, the same in UTM 15N in US feet, and with
# the false easting and northing moved by 0.1 m.
.rf_frames <- function(m) {
  list(m = m,
       usft = sf::st_transform(m, sf::st_crs(
         "+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")),
       m01 = sf::st_transform(m, sf::st_crs(paste(
         "+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1 +y_0=0.1",
         "+datum=NAD83 +units=m +no_defs"))))
}

# The item's lattice: a predictor w with a trend in x, and z of a smooth
# field plus 0.5 w plus noise of sd `noise`.  `s` scales the spacing and the
# field with it (0.1: 25 m, a layer under 1 km across).
.rf_lattice <- function(nx, ny, noise, seed = 1, s = 1) {
  g <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  x <- 4.5e5 + g$i * 250 * s; y <- 4.55e6 + g$j * 250 * s
  withr::with_seed(seed, {
    w <- (x - min(x)) / (1000 * s) + rnorm(length(x))
    z <- sin((x - min(x)) / (1500 * s)) + cos((y - min(y)) / (1250 * s)) +
      0.5 * w + rnorm(length(x), 0, noise)
  })
  sf::st_as_sf(data.frame(x = x, y = y, z = z, w = w), coords = c("x", "y"),
               crs = 26915)
}

# The lattice of test 7 of test-fix-reml-offset.R, whose REML range is
# identified, with a jitter of up to `jitter` metres for the control.
.rf_grid <- function(nx, ny, seed = 1, jitter = 0) {
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
  sf::st_as_sf(data.frame(x = x, y = y, zz = z2 + 0.5 * w, w = w),
               coords = c("x", "y"), crs = 26915)
}

# The estimate and the warnings it gave.
.rf_run <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

test_that("a REML range that runs off gives the same detrend_method, warning and nugget in every unit and origin", {
  skip_if_not_installed("nlme")
  # (nx, ny, noise sd, reml_max_n).  Before, the three frames gave
  # reml/ols/ols, ols/ols/reml and ols/reml/reml on the first three, and the
  # fourth's nugget moved in the seventh digit.
  cfgs <- list(c(20, 15, 0.5, 100), c(20, 15, 0.5, 150), c(24, 18, 1, 120),
               c(20, 15, 1, 150))
  for (cfg in cfgs) {
    L <- .rf_frames(.rf_lattice(cfg[1], cfg[2], cfg[3]))
    est <- lapply(L, function(l) .rf_run(estimate_sac_range(
      l, "z", predictor_vars = "w", detrend = "reml", reml_max_n = cfg[4])))
    meth <- vapply(est, function(e) attr(e$value, "detrend_method"), "")
    fell <- vapply(est, function(e) any(grepl("REML detrending on w", e$warnings)),
                   TRUE)
    nug <- vapply(est, function(e) attr(e$value, "nugget"), 0)
    info <- paste(cfg, collapse = ", ")
    expect_identical(unname(meth[c("usft", "m01")]), unname(meth[c("m", "m")]),
                     info = info)
    expect_identical(unname(fell[c("usft", "m01")]), unname(fell[c("m", "m")]),
                     info = info)
    expect_identical(unname(fell[["m"]]), unname(meth[["m"]] == "ols"), info = info)
    expect_equal(nug[["usft"]], nug[["m"]], tolerance = 1e-8, info = info)
    expect_equal(nug[["m01"]], nug[["m"]], tolerance = 1e-8, info = info)
    expect_true(all(is.na(vapply(est, function(e) as.numeric(e$value), 0))),
                info = info)
  }
})

test_that("gls() is given the same coordinates and starting ranges in every unit and origin", {
  skip_if_not_installed("nlme")
  L <- .rf_frames(.rf_lattice(20, 15, 0.5))
  seen <- function(l) {
    got <- list()
    local_mocked_bindings(gls = function(model, data, correlation = NULL, ...) {
      got[[length(got) + 1L]] <<- list(
        x = data$.sac_x, y = data$.sac_y, z = data$z,
        start = if (is.null(correlation)) NULL else as.numeric(correlation))
      stop("not fitted")
    }, .package = "nlme")
    suppressWarnings(estimate_sac_range(l, "z", predictor_vars = "w",
                                        detrend = "reml", reml_max_n = 100L))
    got
  }
  a <- lapply(L, seen)
  # The grid start and the three extent-based starts; with no fit there is
  # no null model to compare it with.
  expect_length(a$m, 4L)
  for (nm in c("usft", "m01")) {
    expect_identical(a[[nm]], a$m, info = nm)
  }
  # The coordinates are those of the frame: centred and in diagonals.
  expect_lte(max(abs(c(a$m[[1]]$x, a$m[[1]]$y))), 0.5)
  expect_gt(diff(range(a$m[[1]]$x)), 0.5)
})

test_that(".reml_trend() returns the range and pair_floor in the caller's unit", {
  skip_if_not_installed("nlme")
  L <- .rf_frames(.rf_grid(40, 30))
  ext <- c(m = 9750, usft = 9750 * 3937 / 1200, m01 = 9750)
  r <- lapply(names(L), function(nm) {
    xy <- sf::st_coordinates(L[[nm]])
    mf <- data.frame(zz = L[[nm]]$zz, w = L[[nm]]$w, .sac_x = xy[, 1], .sac_y = xy[, 2])
    .reml_trend(mf, zz ~ w, extent = ext[[nm]], max_n = 200L, seed = 123L)
  })
  names(r) <- names(L)
  usft <- 3937 / 1200
  # The range parameter is short of the layer: identified, in metres.
  expect_true(r$m$range > 50 && r$m$range < 2000)
  expect_equal(r$usft$range, r$m$range * usft, tolerance = 1e-8)
  expect_equal(r$m01$range, r$m$range, tolerance = 1e-8)
  # The 30th-shortest distance of the rows used is the lattice spacing.
  expect_equal(r$m$pair_floor, 250, tolerance = 1e-8)
  expect_equal(r$usft$pair_floor, 250 * usft, tolerance = 1e-8)
  for (nm in c("usft", "m01")) {
    expect_equal(r[[nm]]$nugget_prop, r$m$nugget_prop, tolerance = 1e-6, info = nm)
    expect_equal(r[[nm]]$sigma2, r$m$sigma2, tolerance = 1e-8, info = nm)
    expect_equal(r[[nm]]$beta, r$m$beta, tolerance = 1e-8, info = nm)
    expect_equal(r[[nm]]$logLik, r$m$logLik, tolerance = 1e-8, info = nm)
  }
  # And the effective range of estimate_sac_range() in US feet.
  e <- lapply(L[c("m", "usft")], function(l) suppressWarnings(estimate_sac_range(
    l, "zz", predictor_vars = "w", detrend = "reml", reml_max_n = 200L)))
  expect_identical(attr(e$usft, "detrend_method"), "reml")
  expect_true(is.finite(as.numeric(e$m)))
  expect_equal(as.numeric(e$usft), as.numeric(e$m) * usft, tolerance = 1e-8)
})

test_that("a jittered lattice keeps its REML range", {
  # The control: no ties, an identified range, and the same answer as
  # before the change to rounding (632.98363 m before it).
  skip_if_not_installed("nlme")
  j <- .rf_grid(40, 30, jitter = 60)
  e <- suppressWarnings(estimate_sac_range(j, "zz", predictor_vars = "w",
                                           detrend = "reml", reml_max_n = 150L))
  expect_identical(attr(e, "detrend_method"), "reml")
  expect_equal(as.numeric(e), 632.983630502, tolerance = 1e-6)
  expect_equal(attr(e, "reml")$sigma2, 0.369090647275, tolerance = 1e-6)
})

# The run-off lattices of test 1 and the one of test 4 of
# test-fix-reml-offset.R, (nx, ny, noise sd, reml_max_n).
.rf_cfgs <- list(c(20, 15, 0.5, 100), c(20, 15, 0.5, 150), c(24, 18, 1, 120),
                 c(20, 15, 1, 150), c(20, 15, 1, 100))

# detrend_method, the fallback warning and the nugget divided by `sc` of
# each layer, compared with those of the first.
.rf_same <- function(layers, sc, max_n, info) {
  est <- lapply(layers, function(l) .rf_run(estimate_sac_range(
    l, "z", predictor_vars = "w", detrend = "reml", reml_max_n = max_n)))
  meth <- vapply(est, function(e) attr(e$value, "detrend_method"), "")
  fell <- vapply(est, function(e) any(grepl("REML detrending on w", e$warnings)),
                 TRUE)
  nug <- vapply(est, function(e) attr(e$value, "nugget"), 0) / sc
  for (k in seq_along(layers)[-1L]) {
    ik <- paste(info, names(layers)[k])
    expect_identical(meth[[k]], meth[[1L]], info = ik)
    expect_identical(fell[[k]], fell[[1L]], info = ik)
    expect_equal(nug[[k]], nug[[1L]], tolerance = 1e-8, info = ik)
  }
  expect_identical(unname(fell), unname(meth == "ols"), info = info)
}

test_that("a layer under 1 km across gives the same REML outcome in metres and in km", {
  # The extent that bounds the starting ranges was floored at 1 in the
  # caller's unit: 475 in metres and 1 (for 0.475) in km.  The first
  # lattice gave ols / reml and the fourth reml / ols.
  skip_if_not_installed("nlme")
  km <- sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=km +no_defs")
  for (cfg in .rf_cfgs) {
    m <- .rf_lattice(cfg[1], cfg[2], cfg[3], s = 0.1)
    .rf_same(list(m = m, km = sf::st_transform(m, km)), c(1, 1), cfg[4],
             paste(cfg, collapse = ", "))
  }
})

test_that("the REML outcome does not move with the unit or the last bits of the response", {
  # Before, z, 1000 z, z (1 + 2 eps) and z / 1000 gave ols / reml / reml /
  # reml on the first lattice, reml / reml / reml / ols on the second and
  # reml / reml / ols / reml on the fifth.
  skip_if_not_installed("nlme")
  eps <- .Machine$double.eps
  for (cfg in .rf_cfgs) {
    m <- .rf_lattice(cfg[1], cfg[2], cfg[3])
    v <- list(as_is = m, k = within(m, z <- z * 1000),
              ulp = within(m, z <- z * (1 + 2 * eps)),
              milli = within(m, z <- z / 1000))
    .rf_same(v, c(1, 1e6, 1, 1e-6), cfg[4], paste(cfg, collapse = ", "))
  }
})

test_that("the REML outcome does not move with the unit or the last bits of a predictor", {
  # Before, w, 1000 w, w (1 + 2 eps) and w / 1000 gave reml / ols / reml /
  # ols on the second lattice, reml / reml / ols / reml on the fourth and
  # reml / ols / reml / reml on the fifth.
  skip_if_not_installed("nlme")
  eps <- .Machine$double.eps
  for (cfg in .rf_cfgs) {
    m <- .rf_lattice(cfg[1], cfg[2], cfg[3])
    v <- list(as_is = m, k = within(m, w <- w * 1000),
              ulp = within(m, w <- w * (1 + 2 * eps)),
              milli = within(m, w <- w / 1000))
    .rf_same(v, rep(1, 4), cfg[4], paste(cfg, collapse = ", "))
  }
})

test_that(".reml_trend() returns the coefficients, variance and log-likelihoods of the data as given", {
  # The fit reads the response and the predictor in the frame of their
  # range; what it returns is that of gls() on the data as given at the
  # fitted correlation.  A response at 50 in a unit 1000 times smaller, a
  # predictor in another unit at an offset, and a factor left as it is.
  skip_if_not_installed("nlme")
  l <- .rf_grid(15, 12)
  xy <- sf::st_coordinates(l)
  f <- factor(rep(c("a", "b", "c"), length.out = nrow(l)))
  mf <- data.frame(zz = 50 + l$zz * 1000, w = 7 + l$w / 1000, f = f,
                   .sac_x = xy[, 1], .sac_y = xy[, 2])
  r <- .reml_trend(mf, zz ~ w + f, extent = 3500, max_n = 400L, seed = 123L)
  expect_false(is.null(r))
  expect_false(r$subsampled)
  g <- nlme::gls(zz ~ w + f, data = mf, method = "REML",
                 correlation = nlme::corExp(value = c(r$range, r$nugget_prop),
                                            form = ~ .sac_x + .sac_y,
                                            nugget = TRUE, fixed = TRUE))
  expect_equal(r$beta, stats::coef(g), tolerance = 1e-7)
  expect_equal(r$sigma2, g$sigma^2, tolerance = 1e-7)
  expect_equal(r$logLik, as.numeric(stats::logLik(g)), tolerance = 1e-7)
  g0 <- nlme::gls(zz ~ w + f, data = mf, method = "REML")
  expect_equal(r$logLik_null, as.numeric(stats::logLik(g0)), tolerance = 1e-7)
})

test_that(".reml_trend() reads an extent that is not positive as the diagonal", {
  skip_if_not_installed("nlme")
  l <- .rf_grid(15, 12)
  xy <- sf::st_coordinates(l)
  mf <- data.frame(zz = l$zz, w = l$w, .sac_x = xy[, 1], .sac_y = xy[, 2])
  dg <- sqrt(sum((apply(xy[, 1:2], 2, max) - apply(xy[, 1:2], 2, min))^2))
  a <- .reml_trend(mf, zz ~ w, extent = dg, max_n = 400L)
  for (e in list(0, NA_real_, -1)) {
    b <- .reml_trend(mf, zz ~ w, extent = e, max_n = 400L)
    expect_identical(b, a, info = format(e))
  }
})

test_that(".kmeans_frame() returns the centre and the scale of its frame", {
  xy <- cbind(c(4.5e5, 4.51e5, 4.505e5), c(4.55e6, 4.552e6, 4.553e6))
  fr <- .kmeans_frame(xy)
  expect_identical(fr$mid, c(4.5e5, 4.55e6) / 2 + c(4.51e5, 4.553e6) / 2)
  expect_identical(fr$scale, sqrt(1e3^2 + 3e3^2))
  # Rounded to 2^-30 of the diagonal.
  expect_lte(max(abs(fr$z - sweep(xy, 2L, fr$mid) / fr$scale)), 2^-31)
  one <- .kmeans_frame(cbind(c(5, 5), c(7, 7)))
  expect_identical(one$scale, 1)
})
