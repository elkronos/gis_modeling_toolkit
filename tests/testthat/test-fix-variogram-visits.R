# tests/testthat/test-fix-variogram-visits.R
# ---------------------------------------------------------------------------
# Rows that share a location (repeat visits) are read for what they are.
#
# estimate_sac_range() sorts its rows by coordinates and then by the
# response, so the rows of one location arrive in the order of their values.
# Three steps read meaning into that order or into the rows being neighbours:
#   - the REML fit kept "the first" row of each location, which is its
#     smallest response: with four visits of unit noise the intercept came
#     out 1.04 low and the nugget at 0.41 for a truth of 1;
#   - the short-lag semivariance paired each row with the next at its
#     location, which pairs neighbouring order statistics: 0.205 for a
#     variance between visits of 0.300 at 8 visits a location;
#   - the test for spatial structure took the other visits to a location as
#     its nearest neighbours, and found "structure" in 40 of 40 layers of
#     100 locations with independent effects, which then got a range (median
#     96 against 47 between locations) and a nugget of 0, because the pairs
#     at zero separation were pooled into the first lag bin.
# ---------------------------------------------------------------------------

# `sites` locations on a 1000 m square, `visits` rows at each.  The location
# effect has variance 1 and effective range `R` (independent when R is NULL);
# the variance between visits is `noise`; the mean is `mu`.
fvv_layer <- function(seed, sites = 100L, visits = 5L, R = NULL, noise = 0.3, mu = 0) {
  set.seed(seed)
  sx <- runif(sites, 0, 1000); sy <- runif(sites, 0, 1000)
  eff <- if (is.null(R)) rnorm(sites) else
    drop(crossprod(chol(exp(-3 * as.matrix(stats::dist(cbind(sx, sy))) / R) +
                          diag(1e-8, sites)), rnorm(sites)))
  m <- if (length(visits) == 1L) rep(visits, sites) else visits
  site <- rep(seq_len(sites), m)
  d <- data.frame(x = 5e5 + sx[site], y = 5e6 + sy[site], site = site)
  d$w <- rnorm(nrow(d))
  d$z <- mu + eff[site] + rnorm(nrow(d), sd = sqrt(noise))
  d
}
fvv_sf <- function(d) sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)

test_that("the REML fit takes one visit a location at random, not its smallest response", {
  skip_if_not_installed("nlme")
  rt <- spatialkit:::.reml_trend
  b0 <- nug <- numeric(8)
  for (i in 1:8) {
    d <- fvv_layer(8300 + i, visits = 5L, R = 300, noise = 1, mu = 10)
    mf <- data.frame(z = d$z, w = d$w, .sac_x = d$x, .sac_y = d$y)
    # The order estimate_sac_range() hands the rows over in: by coordinates,
    # then by the response, so each location's smallest response is first.
    mf <- mf[order(mf$.sac_x, mf$.sac_y, mf$z), ]
    fit <- rt(mf, z ~ w, extent = 1000)
    expect_identical(fit$n_used, 100L)
    expect_identical(fit$n_colocated, 400L)
    b0[i] <- fit$beta[["(Intercept)"]]
    nug[i] <- fit$sigma2 * fit$nugget_prop
    if (i <= 2L) {
      # The same rows in any order give the same fit: the rows are put in a
      # fixed order before one is drawn at each location.
      set.seed(i); perm <- sample(nrow(mf))
      again <- rt(mf[perm, ], z ~ w, extent = 1000)
      expect_identical(again$beta, fit$beta)
      expect_identical(again$range, fit$range)
      # Another seed draws other visits.
      other <- rt(mf, z ~ w, extent = 1000, seed = 7L)
      expect_false(isTRUE(all.equal(other$beta, fit$beta)))
    }
  }
  # The intercept is 10: the eight fits average 9.93.  The smallest of five
  # visits of unit noise lies 1.16 below its location's level, and fits to
  # those average 8.81.
  expect_gt(mean(b0), 9.6)
  expect_lt(mean(b0), 10.4)
  # The nugget is the variance between visits, 1.  A single fit of 100
  # locations places it anywhere from 0 to 1.4; the eight average 0.74,
  # against 0.39 for fits to the minima (a fit to the means of five visits
  # would aim at a fifth).
  expect_gt(mean(nug), 0.55)
})

test_that("a layer with no repeated location gets the REML subsample it always got", {
  skip_if_not_installed("nlme")
  rt <- spatialkit:::.reml_trend
  set.seed(77); n <- 90
  mf <- data.frame(.sac_x = runif(n, 0, 1000), .sac_y = runif(n, 0, 1000), w = rnorm(n))
  mf$z <- 2 + 0.5 * mf$w + rnorm(n)
  set.seed(5); before <- runif(1)
  set.seed(5); fit <- rt(mf, z ~ w, extent = 1000, max_n = 60L); after <- runif(1)
  expect_identical(after, before)                 # the caller's stream is untouched
  expect_identical(fit$n_used, 60L)
  expect_identical(fit$n_colocated, 0L)
  expect_true(fit$subsampled)
  # The rows the fit used are those of a seeded draw over the sorted rows.
  srt <- mf[order(mf$.sac_x, mf$.sac_y, mf$w, mf$z), ]
  set.seed(123); idx <- sample.int(n, 60L)
  D <- as.matrix(stats::dist(srt[idx, c(".sac_x", ".sac_y")]))
  expect_equal(fit$pair_floor, sort(D[upper.tri(D)])[30])
})

test_that("estimate_sac_range(detrend = 'reml') on repeat visits fits the locations, and says so", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  d <- fvv_layer(8101, sites = 120L, visits = 4L, R = 300, noise = 1, mu = 10)
  d$z <- d$z + 0.5 * d$w
  set.seed(3); d <- d[sample(nrow(d)), ]
  lines <- capture_spatialkit_log(
    r <- suppressWarnings(estimate_sac_range(fvv_sf(d), "z", "w", detrend = "reml",
                                             structure_alpha = 1)))
  expect_identical(attr(r, "detrend_method"), "reml")
  expect_identical(attr(r, "reml")$n_used, 120L)
  expect_identical(attr(r, "n_used"), 480L)       # every row is detrended and binned
  expect_true(log_has(lines, "120 location\\(s\\) have more than one row"))
  expect_true(log_has(lines, "one row a location, drawn at random under `seed` \\(120 of 480"))
  # The fit behind it is .reml_trend()'s on those rows: the draw is made on
  # the rows in their sorted order, whatever order they came in.
  mf <- data.frame(z = d$z, w = d$w, .sac_x = d$x, .sac_y = d$y)
  direct <- spatialkit:::.reml_trend(
    mf, z ~ w, extent = max(diff(range(d$x)), diff(range(d$y))))
  expect_equal(attr(r, "reml")$sigma2, direct$sigma2, tolerance = 1e-6)
  expect_equal(attr(r, "reml")$nugget_prop, direct$nugget_prop, tolerance = 1e-6)
  expect_equal(attr(r, "reml")$logLik, direct$logLik, tolerance = 1e-6)
})

test_that("the semivariance at zero separation is the variance between visits, from every pair", {
  skip_if_not_installed("gstat")
  for (v in c(2L, 5L, 8L)) {
    d <- fvv_layer(9900 + v, visits = v)
    r <- suppressWarnings(estimate_sac_range(fvv_sf(d), "z"))
    within <- mean(tapply(d$z, d$site, stats::var))
    # With as many visits at every location, the mean over all pairs of rows
    # at one location is the mean of the locations' variances, exactly.
    expect_equal(attr(r, "zero_lag_semivariance"), within, tolerance = 1e-10)
    expect_identical(attr(r, "zero_lag_pairs"), as.integer(100 * v * (v - 1) / 2))
    # It was 0.3 by construction: within sampling error of it, not short of
    # it as the order-statistic pairs were.
    expect_lt(abs(attr(r, "zero_lag_semivariance") - 0.3), 0.08)
    # The nearest-neighbour pairs are between locations: one pair for each
    # location and its nearest, counted once, at a positive separation.
    expect_lte(attr(r, "short_lag_pairs"), 100L)
    expect_gte(attr(r, "short_lag_pairs"), 50L)
    expect_gt(attr(r, "short_lag_dist"), 0)
    # They carry the variation between locations as well, so they bound the
    # zero-lag value from above.
    expect_gt(attr(r, "short_lag_semivariance"), attr(r, "zero_lag_semivariance"))
  }
  # A layer with no repeated location has no such pairs.
  one <- suppressWarnings(estimate_sac_range(fvv_sf(fvv_layer(9911, visits = 1L, R = 300)), "z"))
  expect_identical(attr(one, "zero_lag_pairs"), 0L)
  expect_true(is.na(attr(one, "zero_lag_semivariance")))
})

test_that(".sac_short_lag() matches its definition with unequal numbers of visits", {
  set.seed(21); sites <- 25L
  m <- sample(1:4, sites, replace = TRUE)
  d <- fvv_layer(22, sites = sites, visits = m)
  xy <- cbind(d$x, d$y); v <- d$z
  got <- spatialkit:::.sac_short_lag(xy, v)
  # Zero separation: every pair of rows at one location.
  D <- as.matrix(stats::dist(xy)); ut <- upper.tri(D)
  half <- outer(v, v, "-")^2 / 2
  expect_identical(got$zero_pairs, sum(ut & D == 0))
  expect_equal(got$zero_gamma, mean(half[ut & D == 0]))
  expect_equal(got$zero_ss, sum(half[ut & D == 0]))
  # Nearest neighbours: each location and its nearest other location, each
  # such pair once; every pairing of their rows, the pair of locations
  # weighing one.
  sxy <- unique(xy); S <- as.matrix(stats::dist(sxy)); diag(S) <- Inf
  nn <- apply(S, 1, which.min)
  pr <- unique(t(apply(cbind(seq_len(nrow(sxy)), nn), 1, sort)))
  key <- function(a) paste(a[, 1], a[, 2])
  loc <- match(key(xy), key(sxy))
  g <- vapply(seq_len(nrow(pr)), function(k)
    mean(half[loc == pr[k, 1], loc == pr[k, 2]]), numeric(1))
  expect_equal(got$gamma, mean(g))
  expect_identical(got$pairs, nrow(pr))
  expect_equal(got$dist, stats::median(S[pr]))
  # The order of the rows changes nothing: in particular, sorting the rows
  # of each location by their value does not.
  o <- order(xy[, 1], xy[, 2], v)
  srt <- spatialkit:::.sac_short_lag(xy[o, ], v[o])
  expect_equal(srt[c("gamma", "zero_gamma")], got[c("gamma", "zero_gamma")])
  expect_identical(srt[c("pairs", "zero_pairs")], got[c("pairs", "zero_pairs")])
  # One location: nothing at a positive separation, the visits still count.
  one <- spatialkit:::.sac_short_lag(cbind(rep(1, 4), rep(2, 4)), c(1, 2, 4, 7))
  expect_true(is.na(one$gamma)); expect_identical(one$pairs, 0L)
  expect_identical(one$zero_pairs, 6L)
  expect_equal(one$zero_gamma, stats::var(c(1, 2, 4, 7)))
})

test_that("pairs at zero separation are in no bin of the variogram", {
  skip_if_not_installed("gstat")
  d <- fvv_layer(9905, sites = 80L, visits = 4L, R = 300)
  r <- suppressWarnings(estimate_sac_range(fvv_sf(d), "z", lags = "equal", n_lags = 10,
                                           structure_alpha = 1))
  vg <- attr(r, "variogram")
  # The variogram by hand, over pairs at a positive separation only.
  D <- as.matrix(stats::dist(cbind(d$x, d$y))); ut <- upper.tri(D)
  half <- outer(d$z, d$z, "-")^2 / 2
  b <- attr(r, "lag_boundaries")
  lo <- c(0, b[-length(b)])
  hand <- t(vapply(seq_along(b), function(i) {
    k <- ut & D > lo[i] & D <= b[i]
    c(sum(k), mean(D[k]), mean(half[k]))
  }, numeric(3)))
  expect_equal(vg$np, hand[, 1])
  expect_equal(vg$dist, hand[, 2], tolerance = 1e-8)
  expect_equal(vg$gamma, hand[, 3], tolerance = 1e-8)
  # 480 pairs of visits were at zero separation; gstat puts them in the
  # first bin, where they pulled its mean distance and semivariance down.
  expect_identical(attr(r, "zero_lag_pairs"), 480L)
  expect_equal(attr(r, "first_lag"), hand[1, 2], tolerance = 1e-8)
  # The same holds for each direction, and for the adaptive bins.
  ra <- suppressWarnings(estimate_sac_range(fvv_sf(d), "z", structure_alpha = 1,
                                            keep_directional_fits = TRUE))
  va <- attr(ra, "variogram")
  expect_equal(sum(va$np), sum(ut & D > 0 & D <= attr(ra, "cutoff_dist")))
  expect_gt(min(va$dist), min(D[ut & D > 0]) - 1e-9)
  for (az in names(attr(ra, "directional_fits"))) {
    vd <- attr(ra, "directional_fits")[[az]]$variogram
    if (!is.null(vd)) expect_gt(min(vd$dist), min(D[ut & D > 0]) - 1e-9)
  }
})

test_that(".sac_drop_zero_lag() takes the pairs out of the first row and leaves the rest", {
  f <- spatialkit:::.sac_drop_zero_lag
  vg <- data.frame(np = c(100, 50), dist = c(6, 20), gamma = c(0.7, 1.2))
  # 40 pairs at zero separation with half squared differences summing to 12
  # (a semivariance of 0.3): the 60 others sit at 10 with semivariance 29/30.
  out <- f(vg, zero_pairs = 40, zero_ss = 12, first_bound = 15)
  expect_equal(out$np, c(60, 50))
  expect_equal(out$dist, c(10, 20))
  expect_equal(out$gamma, c((70 - 12) / 60, 1.2))
  # Nothing to take out.
  expect_identical(f(vg, 0L, 0, 15), vg)
  # A first row with nothing but those pairs goes.
  only <- data.frame(np = c(40, 50), dist = c(0, 20), gamma = c(0.3, 1.2))
  expect_equal(f(only, 40, 12, 15)$np, 50)
  # A first row that cannot hold them (fewer pairs, or a remainder outside
  # the bin or below zero) is left as it is.
  expect_identical(f(vg, 400, 12, 15), vg)
  expect_identical(f(vg, 40, 12, first_bound = 8), vg)
  expect_identical(f(vg, 40, 500, 15), vg)
})

test_that("location effects alone are not spatial structure", {
  skip_if_not_installed("gstat")
  # 100 locations visited 5 times, independent location effects (variance 1)
  # and a variance between visits of 0.3.  Rows at one location are alike;
  # locations near each other are not.
  res <- lapply(1:12, function(i)
    suppressWarnings(estimate_sac_range(fvv_sf(fvv_layer(9900 + i)), "z")))
  p <- vapply(res, function(r) attr(r, "structure_p"), numeric(1))
  expect_lte(sum(p < 0.05), 2L)                   # it was 12 of 12
  refused <- vapply(res, function(r) identical(
    attr(r, "rejected_reason"), "no spatial structure detected at the nearest-neighbour scale"),
    logical(1))
  expect_gte(sum(refused), 10L)
  expect_lte(sum(vapply(res, is.finite, logical(1))), 1L)
  # The test is the one made on the location means.
  d <- fvv_layer(9901)
  xy <- cbind(d$x, d$y)
  first <- !duplicated(d$site)
  on_means <- spatialkit:::.sac_structure_test(xy[first, ], as.numeric(tapply(d$z, d$site, mean)))
  on_rows  <- spatialkit:::.sac_structure_test(xy, d$z)
  expect_equal(on_rows$p, on_means$p, tolerance = 1e-10)
  expect_equal(on_rows$table$morans_i, on_means$table$morans_i, tolerance = 1e-10)
  expect_identical(on_rows$n_locations, 100L)
  expect_equal(attr(res[[1]], "structure_p"), on_means$p, tolerance = 1e-8)
  # With a trend's design, the design is averaged over each location too.
  X <- cbind(1, d$w)
  with_x <- spatialkit:::.sac_structure_test(xy, d$z, X)
  Xm <- cbind(1, as.numeric(tapply(d$w, d$site, mean)))
  expect_equal(with_x$p, spatialkit:::.sac_structure_test(
    xy[first, ], as.numeric(tapply(d$z, d$site, mean)), Xm)$p, tolerance = 1e-10)
  # The refusal says what the test was made on.
  lines <- capture_spatialkit_log(estimate_sac_range(fvv_sf(d), "z"))
  expect_true(log_has(lines, "The 500 rows lie at 100 locations, and the test is made on the location means"))
})

test_that("structure between visited locations is still found, with its range", {
  skip_if_not_installed("gstat")
  # The location effect has an effective range of 300.
  res <- lapply(1:6, function(i)
    suppressWarnings(estimate_sac_range(fvv_sf(fvv_layer(9950 + i, R = 300)), "z")))
  p <- vapply(res, function(r) attr(r, "structure_p"), numeric(1))
  expect_true(all(p < 0.05))
  rng <- vapply(res, as.numeric, numeric(1))
  # 501, 343, 226, 461 and 343, and one refused as past the lags.
  expect_gte(sum(is.finite(rng)), 4L)
  expect_gt(stats::median(rng, na.rm = TRUE), 120)
  expect_lt(stats::median(rng, na.rm = TRUE), 700)
  # The same layer, rows in another order, gives the identical object.
  d <- fvv_layer(9951, R = 300)
  set.seed(4); perm <- sample(nrow(d))
  expect_identical(suppressWarnings(estimate_sac_range(fvv_sf(d[perm, ]), "z")), res[[1]])
})

test_that("the nugget is not let under the variance between visits", {
  skip_if_not_installed("gstat")
  # 100 locations visited 5 times, a location effect of range 300 and a
  # variance between visits of 0.3.  Extrapolated from the positive lags
  # alone, the nugget came out under that variance in 22 of 27 fits and
  # exactly 0 in 7, while 1000 pairs of visits measure it at 0.28 to 0.32.
  res <- lapply(1:8, function(i)
    suppressWarnings(estimate_sac_range(fvv_sf(fvv_layer(9950 + i, R = 300)), "z")))
  for (r in res) {
    z0 <- attr(r, "zero_lag_semivariance")
    expect_identical(attr(r, "nugget_floor"), z0)
    expect_gte(sac_nugget(r), z0)
    # A nugget on its floor is a measurement, not the zero bound.
    expect_identical(attr(r, "nugget_at_bound"), FALSE)
  }
  nug <- vapply(res, sac_nugget, numeric(1))
  expect_lt(abs(stats::median(nug) - 0.3), 0.08)
  # Where the floor holds, the nugget is the floor to the last bit.
  on_floor <- vapply(res, function(r) sac_nugget(r) == attr(r, "nugget_floor"), logical(1))
  expect_true(any(on_floor))
  # Fewer than 30 rows beyond one a location: too little to bound with.
  d <- fvv_layer(9960, sites = 200L, visits = 1L, R = 300)
  few <- rbind(d, d[1:20, ])
  few$z[201:220] <- few$z[201:220] + rnorm(20, sd = 0.5)
  rf <- suppressWarnings(estimate_sac_range(fvv_sf(few), "z"))
  expect_identical(attr(rf, "zero_lag_pairs"), 20L)
  expect_true(is.finite(attr(rf, "zero_lag_semivariance")))
  expect_identical(attr(rf, "nugget_floor"), 0)
  # No repeated location: no floor.
  expect_identical(attr(suppressWarnings(estimate_sac_range(fvv_sf(d), "z")), "nugget_floor"), 0)
})

test_that(".sac_fit_model() and .sac_fit_nested() hold the nugget to a floor exactly", {
  skip_if_not_installed("gstat")
  fit <- spatialkit:::.sac_fit_model
  # A variogram that rises from near zero: nugget 0.02, partial sill 1,
  # effective range 200.
  h <- c(10, 20, 40, 80, 160, 320, 500)
  vg <- data.frame(np = c(60, 120, 300, 800, 2000, 5000, 8000), dist = h,
                   gamma = 0.02 + (1 - exp(-3 * h / 200)))
  free <- fit(vg, cutoff_dist = 500)
  expect_equal(free$psill[1], 0.02, tolerance = 1e-6)
  expect_equal(3 * free$range[2], 200, tolerance = 1e-5)
  # A floor under the fitted nugget changes nothing.
  expect_equal(fit(vg, 500, nugget_min = 0.01)$psill, free$psill, tolerance = 1e-10)
  # A floor above it holds, to the last bit, and the fit is the best one
  # with the nugget there: no better at any other range.
  held <- fit(vg, 500, nugget_min = 0.1)
  expect_identical(held$psill[1], 0.1)
  expect_gt(attr(held, "SSErr"), attr(free, "SSErr"))
  sse <- function(r, c1) {
    w <- vg$np / (h / 500)^2; w <- w / sum(w)
    gbar <- sum(vg$np * vg$gamma) / sum(vg$np)
    sum(w * (vg$gamma / gbar - 0.1 / gbar - c1 / gbar * (1 - exp(-3 * h / r)))^2)
  }
  at <- sse(3 * held$range[2], held$psill[2])
  expect_equal(at, attr(held, "SSErr"), tolerance = 1e-8)
  for (r in 3 * held$range[2] * c(0.8, 1.25)) for (c1 in held$psill[2] * c(0.9, 1, 1.1))
    expect_gte(sse(r, c1), at * (1 - 1e-9))
  # A floor above every semivariance leaves the constant at the floor.
  flat <- fit(vg, 500, nugget_min = 5)
  expect_identical(as.character(flat$model), "Nug")
  expect_identical(flat$psill, 5)
  # One above the weighted mean (0.37 here, most of the weight being on the
  # short lags) and under the long lags: the nugget sits on it, and the
  # structure takes what is left.
  mid <- fit(vg, 500, nugget_min = 0.5)
  expect_identical(mid$psill[1], 0.5)
  expect_gt(mid$psill[2], 0)
  # The two-structure fit takes the same floor.
  nest <- spatialkit:::.sac_fit_nested
  g2 <- 0.02 + 0.5 * (1 - exp(-3 * h / 30)) + 0.5 * (1 - exp(-3 * h / 400))
  vg2 <- transform(vg, gamma = g2)
  a <- nest(vg2, 500, 30, 400)
  expect_equal(a$psill, c(0.02, 0.5, 0.5), tolerance = 1e-6)
  b <- nest(vg2, 500, 30, 400, nugget_min = 0.1)
  expect_identical(b$psill[1], 0.1)
  expect_true(all(b$psill[2:3] > 0))
  expect_equal(nest(vg2, 500, 30, 400, nugget_min = 0.01)$psill, a$psill, tolerance = 1e-10)
})
