# tests/testthat/test-fix-variogram-invariance.R
# ---------------------------------------------------------------------------
# estimate_sac_range() is a function of the data: its verdict, range and
# nugget do not change with the unit of the response or of the coordinates,
# with the order of the rows, or with the last bits of the values.
#
# gstat's iterative fit was neither.  On these layers (500 points on a 400 x
# 300 km extent; white noise, and a 25 km range under as much noise) the same
# layer asked eight ways gave more than one verdict in 6 of 10 white-noise
# layers and 2 of 10 short-range ones, the range moved by more than 5 percent
# in 9 and in 6, and gstat printed "a possible solution MIGHT be to scale
# semivariances and/or distances" 291 times to the console.
# ---------------------------------------------------------------------------

# The layers of the review's reproduction, by seed.
fvi_layer <- function(seed, type = c("short", "white"), n = 500L) {
  type <- match.arg(type)
  set.seed(seed)
  x <- stats::runif(n, 0, 400e3); y <- stats::runif(n, 0, 300e3)
  D <- as.matrix(stats::dist(cbind(x, y)))
  z <- switch(type,
    white = stats::rnorm(n),
    short = drop(crossprod(chol(exp(-3 * D / 25e3) + diag(1e-8, n)), stats::rnorm(n))) +
      stats::rnorm(n))
  sf::st_as_sf(data.frame(x = 5e5 + x, y = 4e6 + y, v = z), coords = c("x", "y"),
               crs = 32615)
}
# What a result says, on the scale of the layer as generated: `by` is what
# the response was multiplied by, `len` what the coordinates were.
fvi_says <- function(s, by = 1, len = 1) {
  vm <- attr(s, "variogram_model")
  list(verdict   = if (is.finite(s)) "identified" else attr(s, "rejected_reason"),
       range     = as.numeric(s) / len,
       rejected  = (attr(s, "rejected_range") %||% NA_real_) / len,
       nugget    = sac_nugget(s) / by^2,
       sill      = if (is.data.frame(vm)) sum(vm$psill) / by^2 else NA_real_,
       converged = if (is.data.frame(vm)) attr(vm, "converged") else NA,
       at_bound  = attr(s, "nugget_at_bound"),
       p         = attr(s, "structure_p"))
}
expect_same_answer <- function(a, b, info) {
  expect_identical(b$verdict, a$verdict, info = info)
  expect_identical(b$converged, a$converged, info = info)
  expect_identical(b$at_bound, a$at_bound, info = info)
  expect_equal(b$range, a$range, tolerance = 1e-6, info = info)
  expect_equal(b$rejected, a$rejected, tolerance = 1e-6, info = info)
  expect_equal(b$nugget, a$nugget, tolerance = 1e-6, info = info)
  expect_equal(b$sill, a$sill, tolerance = 1e-6, info = info)
  expect_equal(b$p, a$p, tolerance = 1e-6, info = info)
}
# The layers on which the old answer moved: short-range 2, 5 and 6, white 1.
FVI_LAYERS <- list(list(102, "short"), list(105, "short"), list(106, "short"),
                   list(101, "white"))

test_that("the answer does not depend on the unit of the response", {
  skip_if_not_installed("gstat")
  for (spec in FVI_LAYERS) {
    L <- fvi_layer(spec[[1]], spec[[2]])
    base <- fvi_says(suppressWarnings(estimate_sac_range(L, "v")))
    # 1e-9: a response whose standard deviation fell under 1.5e-8 was called
    # constant, whatever its values, so nanometres got no answer where
    # metres got one.
    for (by in c(1000, 0.001, 1e-9)) {
      d <- L; d$v <- L$v * by
      expect_same_answer(base, fvi_says(suppressWarnings(estimate_sac_range(d, "v")), by = by),
                         info = sprintf("%s %d, response times %g", spec[[2]], spec[[1]], by))
    }
    # A shift costs six digits of the values and nothing else.
    d <- L; d$v <- L$v + 1e6
    expect_same_answer(base, fvi_says(suppressWarnings(estimate_sac_range(d, "v"))),
                       info = sprintf("%s %d, response plus 1e6", spec[[2]], spec[[1]]))
  }
})

test_that("a response is constant when its spread is rounding, in any unit", {
  skip_if_not_installed("gstat")
  L <- fvi_layer(102, "short")
  r <- suppressWarnings(estimate_sac_range(L, "v"))
  expect_true(is.finite(r))
  # A small unit, a large one, and a level far above the spread: all data.
  for (by in c(1e-12, 1e-150, 1e150)) {
    d <- L; d$v <- L$v * by
    rb <- suppressWarnings(estimate_sac_range(d, "v"))
    expect_null(attr(rb, "rejected_reason"), info = format(by))
    expect_equal(as.numeric(rb), as.numeric(r), tolerance = 1e-6, info = format(by))
  }
  d <- L; d$v <- L$v + 1e8          # it was "the response is constant"
  rs <- suppressWarnings(estimate_sac_range(d, "v"))
  expect_null(attr(rs, "rejected_reason"))
  expect_equal(as.numeric(rs), as.numeric(r), tolerance = 1e-4)
  # Constant, and constant up to the last bits, in any unit.
  for (level in c(3, 3e-9, 3e9)) {
    d <- L; d$v <- level
    expect_identical(attr(estimate_sac_range(d, "v"), "rejected_reason"),
                     "the response is constant", info = format(level))
    set.seed(1)
    d$v <- level * (1 + 4 * .Machine$double.eps * stats::rnorm(nrow(L)))
    expect_identical(attr(estimate_sac_range(d, "v"), "rejected_reason"),
                     "the response is constant", info = format(level))
  }
  # Residuals of an exact fit are rounding against the response they are of.
  d <- L; set.seed(2); d$w <- stats::rnorm(nrow(L)); d$v <- 1e-9 * (2 * d$w + 1)
  expect_match(attr(estimate_sac_range(d, "v", "w"), "rejected_reason"),
               "^the residuals on predictor_vars are constant")
})

test_that("the nugget and the sill scale by the square of the multiplier", {
  skip_if_not_installed("gstat")
  L <- fvi_layer(102, "short")
  r1 <- estimate_sac_range(L, "v")
  d <- L; d$v <- 1000 * L$v
  r2 <- estimate_sac_range(d, "v")
  expect_true(is.finite(r1))
  expect_equal(as.numeric(r2), as.numeric(r1), tolerance = 1e-8)
  expect_equal(sac_nugget(r2), 1e6 * sac_nugget(r1), tolerance = 1e-8)
  expect_equal(attr(r2, "variogram_model")$psill, 1e6 * attr(r1, "variogram_model")$psill,
               tolerance = 1e-8)
  expect_equal(attr(r2, "variogram")$gamma, 1e6 * attr(r1, "variogram")$gamma,
               tolerance = 1e-10)
  expect_equal(attr(r2, "short_lag_semivariance"), 1e6 * attr(r1, "short_lag_semivariance"),
               tolerance = 1e-10)
  # The bins hold the same pairs at the same lags.
  expect_identical(attr(r2, "variogram")$np, attr(r1, "variogram")$np)
  expect_equal(attr(r2, "variogram")$dist, attr(r1, "variogram")$dist, tolerance = 1e-12)
})

test_that("the answer does not depend on the unit of the coordinates", {
  skip_if_not_installed("gstat")
  for (spec in FVI_LAYERS) {
    L <- fvi_layer(spec[[1]], spec[[2]])
    base <- fvi_says(suppressWarnings(estimate_sac_range(L, "v")))
    xy <- sf::st_coordinates(L)
    # Kilometres, and US survey feet: neither is an exact binary multiple.
    for (len in c(1 / 1000, 3937 / 1200)) {
      d <- sf::st_as_sf(data.frame(x = xy[, 1] * len, y = xy[, 2] * len, v = L$v),
                        coords = c("x", "y"), crs = NA)
      expect_same_answer(base, fvi_says(suppressWarnings(estimate_sac_range(d, "v")), len = len),
                         info = sprintf("%s %d, coordinates times %g", spec[[2]], spec[[1]], len))
    }
  }
})

test_that("the same rows in another order give the identical object", {
  skip_if_not_installed("gstat")
  for (spec in FVI_LAYERS) {
    L <- fvi_layer(spec[[1]], spec[[2]])
    a <- suppressWarnings(estimate_sac_range(L, "v"))
    set.seed(7000 + spec[[1]])
    b <- suppressWarnings(estimate_sac_range(L[sample(nrow(L)), ], "v"))
    xy <- sf::st_coordinates(L)
    s <- suppressWarnings(estimate_sac_range(L[order(xy[, 1], xy[, 2]), ], "v"))
    info <- sprintf("%s %d", spec[[2]], spec[[1]])
    # Identical, attributes and all: the value, the variogram, the model, the
    # directional sweep and the structure test.
    expect_identical(b, a, info = info)
    expect_identical(s, a, info = info)
  }
  # With predictors the detrending fit sees the rows in the same order too.
  L <- fvi_layer(102, "short")
  set.seed(1); L$w <- rnorm(nrow(L))
  set.seed(2); perm <- sample(nrow(L))
  expect_identical(suppressWarnings(estimate_sac_range(L[perm, ], "v", "w")),
                   suppressWarnings(estimate_sac_range(L, "v", "w")))
})

test_that("lon/lat rows in another order get the same projection and the identical object", {
  skip_if_not_installed("gstat")
  # 400 points over the conterminous United States, a field of effective
  # range 600 km.  The projection for a lon/lat layer is chosen by scoring
  # candidates on rows taken at even spacing down the layer, so the rows had
  # to be sorted before that choice and not after it: as generated this layer
  # got a Lambert azimuthal CRS and a range of 505,530, and in four of six
  # other row orders an Albers one and 506,838.
  set.seed(14); n <- 400
  lon <- runif(n, -120, -70); lat <- runif(n, 25, 50)
  ll <- sf::st_as_sf(data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326)
  xy <- sf::st_coordinates(sf::st_transform(ll, 5070))
  ll$z <- as.numeric(crossprod(chol(0.8 * exp(-3 * as.matrix(stats::dist(xy)) / 600e3) +
                                      diag(0.2, n)), rnorm(n)))
  a <- suppressWarnings(estimate_sac_range(ll, "z"))
  expect_true(is.finite(a))
  expect_false(sf::st_is_longlat(attr(a, "crs")))
  for (s in 1:6) {
    set.seed(s)
    b <- suppressWarnings(estimate_sac_range(ll[sample(n), ], "z"))
    expect_true(attr(b, "crs") == attr(a, "crs"), info = paste("order", s))
    expect_identical(b, a, info = paste("order", s))
  }
  # A layer with no range to return keeps its evidence in any order too.
  set.seed(13)
  wn <- sf::st_as_sf(data.frame(x = runif(300, -120, -70), y = runif(300, 25, 50)),
                     coords = c("x", "y"), crs = 4326)
  wn$z <- rnorm(300) + sf::st_coordinates(wn)[, 2] / 5
  a <- suppressWarnings(estimate_sac_range(wn, "z"))
  set.seed(3)
  b <- suppressWarnings(estimate_sac_range(wn[sample(300), ], "z"))
  expect_identical(b, a)
  expect_identical(attr(b, "max_dist"), attr(a, "max_dist"))
})

test_that("a change in the last bits of the values changes the answer by as little", {
  skip_if_not_installed("gstat")
  for (spec in FVI_LAYERS) {
    L <- fvi_layer(spec[[1]], spec[[2]])
    base <- fvi_says(suppressWarnings(estimate_sac_range(L, "v")))
    set.seed(7000 + spec[[1]])
    d <- L; d$v <- L$v * (1 + 1e-12 * stats::rnorm(nrow(L)))
    expect_same_answer(base, fvi_says(suppressWarnings(estimate_sac_range(d, "v"))),
                       info = sprintf("%s %d", spec[[2]], spec[[1]]))
  }
})

test_that("nothing is printed to the console, whatever the scale of the response", {
  skip_if_not_installed("gstat")
  # gstat::fit.variogram() print()s a hint when a fit is singular and the
  # semivariances are small against the distances, which no quiet setting,
  # suppressWarnings() or suppressMessages() holds back.
  old <- spatialkit_quiet(TRUE)
  on.exit(spatialkit_quiet(old), add = TRUE)
  for (spec in FVI_LAYERS) {
    L <- fvi_layer(spec[[1]], spec[[2]])
    for (by in c(1, 0.1, 0.001)) {
      d <- L; d$v <- L$v * by
      out <- utils::capture.output(r <- suppressWarnings(estimate_sac_range(d, "v")))
      expect_identical(out, character(0),
                       info = sprintf("%s %d, response times %g", spec[[2]], spec[[1]], by))
    }
  }
})

test_that("blocks sized from the range do not change with the unit of the response", {
  skip_if_not_installed("gstat")
  L <- fvi_layer(102, "short")
  d <- L; d$v <- 1000 * L$v
  f1 <- suppressWarnings(suppressMessages(
    make_folds(L, k = 4, method = "block_kfold", auto_range = TRUE, response_var = "v", seed = 1)))
  f2 <- suppressWarnings(suppressMessages(
    make_folds(d, k = 4, method = "block_kfold", auto_range = TRUE, response_var = "v", seed = 1)))
  expect_true(is.finite(f1$params$sac_range))
  expect_equal(as.numeric(f2$params$sac_range), as.numeric(f1$params$sac_range), tolerance = 1e-8)
  expect_identical(c(f2$params$grid_nx, f2$params$grid_ny), c(f1$params$grid_nx, f1$params$grid_ny))
  expect_identical(f2$folds, f1$folds)
})
