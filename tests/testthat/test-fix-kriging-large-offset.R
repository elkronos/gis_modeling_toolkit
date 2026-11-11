# tests/testthat/test-fix-kriging-large-offset.R
# ---------------------------------------------------------------------------
# kriging_adequacy() takes a response whose raw values lose their spread to
# rounding (1e15 plus a spread of 1, stored to an eighth) about the centre
# estimate_sac_range() regresses it about, so the replicate means, the
# kriging, the plain means and the cross-validation read the spread and not
# the rounding of the offset.  kr_shift, its printed counts and the
# cross-validation statistic are those of the same stored values at 0;
# `mean` and `kr_pred` are in the response's units (BUG-60).
# ---------------------------------------------------------------------------

.ko_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# 800 points: an exponential field (range 600 m) plus 0.5 w plus noise, as
# stored at 1e15 and brought back to 0, so both offsets hold the same values.
# With `replicates`, the first 200 locations are visited a second time.
.ko_layer <- function(replicates = FALSE) {
  .ko_seed(5)
  n <- 800
  xy <- data.frame(x = stats::runif(n, 0, 4000), y = stats::runif(n, 0, 4000))
  D <- as.matrix(stats::dist(xy))
  f <- as.numeric(t(chol(exp(-D / 600) + diag(1e-8, n))) %*% stats::rnorm(n))
  xy$w <- stats::rnorm(n)
  xy$z <- f + 0.5 * xy$w + stats::rnorm(n, sd = 0.3)
  L <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32615)
  if (replicates) {
    L2 <- L[1:200, ]
    L2$z <- L2$z + stats::rnorm(200, sd = 0.3)
    L <- rbind(L, L2)
  }
  L$zg <- (1e15 + L$z) - 1e15
  L
}

.ko_cells <- function(L) {
  g <- sf::st_make_grid(L, n = c(6, 6))
  sf::st_sf(poly_id = seq_along(g), geometry = g)
}

.ko_run <- function(L, cells, off, sac) {
  p <- L
  p$zz <- off + p$zg
  a <- suppressWarnings(assign_features_to_polygons(p, cells))
  suppressWarnings(kriging_adequacy(a, "zz", cells, sac = sac))
}

.ko_shift_line <- function(ka) {
  out <- utils::capture.output(print(ka))
  out[grepl("|shift| > 1", out, fixed = TRUE)]
}

.ko_compare <- function(L, info) {
  cells <- .ko_cells(L)
  sac0 <- suppressWarnings(estimate_sac_range(L, "zg"))
  # At 0 the raw values show their spread, so nothing is centred there.
  expect_identical(.regression_centre(L$zg, bounded = FALSE), 0)
  expect_false(.regression_centre(1e15 + L$zg, bounded = FALSE) == 0)
  for (mode in c("sac NULL", "sac given")) {
    inf <- paste(info, mode)
    s <- if (mode == "sac given") sac0 else NULL
    k0 <- .ko_run(L, cells, 0, s)
    k1 <- .ko_run(L, cells, 1e15, s)
    d0 <- sf::st_drop_geometry(k0); d1 <- sf::st_drop_geometry(k1)
    # The evidence: kr_shift, the counts print() gives, the cross-validation.
    expect_identical(is.finite(d1$kr_shift), is.finite(d0$kr_shift), info = inf)
    ok <- is.finite(d0$kr_shift)
    expect_lt(max(abs(d1$kr_shift[ok] - d0$kr_shift[ok])), 1e-6)
    expect_identical(.ko_shift_line(k1), .ko_shift_line(k0), info = inf)
    cv0 <- attr(k0, "cv"); cv1 <- attr(k1, "cv")
    for (nm in c("zscore_var", "zscore_mean", "rmse"))
      expect_lt(abs(cv1[[nm]] - cv0[[nm]]), 1e-8)
    expect_identical(cv1$n_pred, cv0$n_pred, info = inf)
    expect_equal(d1$kr_var, d0$kr_var, tolerance = 1e-6, info = inf)
    expect_equal(d1$plain_var, d0$plain_var, tolerance = 1e-6, info = inf)
    expect_equal(d1$se, d0$se, tolerance = 1e-9, info = inf)
    expect_equal(attr(k1, "replicates")$within_var, attr(k0, "replicates")$within_var,
                 tolerance = 1e-9, info = inf)
    # `mean` and `kr_pred` in the response's units, to the eighth a double
    # holds them to at 1e15.
    expect_lte(max(abs((d1$mean - 1e15) - d0$mean)), 0.125)
    expect_lte(max(abs((d1$kr_pred - 1e15) - d0$kr_pred)), 0.125)
    expect_true(all(abs(d1$mean - 1e15) < 10))
  }
}

test_that("kriging_adequacy() at 1e15 reads the spread, not the rounding of the offset", {
  skip_if_not_installed("gstat")
  .ko_compare(.ko_layer(), "one point a location")
})

test_that("replicate means at 1e15 are those of the values at 0", {
  skip_if_not_installed("gstat")
  L <- .ko_layer(replicates = TRUE)
  .ko_compare(L, "replicates")
})

test_that("a response whose raw values show their spread is not centred", {
  # Only responses that lose their spread are taken about a centre.  At 0 and
  # at 1e6 the raw values show it, so `mean` is the plain mean of the raw
  # response bit for bit, and kr_shift is (kr_pred - mean) over its standard
  # deviation bit for bit: centring every response would round both.
  skip_if_not_installed("gstat")
  L <- .ko_layer()
  cells <- .ko_cells(L)
  for (off in c(0, 1e6)) {
    inf <- paste("offset", off)
    p <- L
    p$zz <- off + p$zg
    a <- suppressWarnings(assign_features_to_polygons(p, cells))
    expect_identical(.regression_centre(a$zz, bounded = FALSE), 0, info = inf)
    d <- sf::st_drop_geometry(suppressWarnings(kriging_adequacy(a, "zz", cells)))
    # The plain mean of the raw values, summed in the order the function puts
    # the rows in (coordinates, cell ID, response), which does not depend on
    # the order given.
    xy <- sf::st_coordinates(a)
    o <- .canonical_order(xy[, 1], xy[, 2], .id_as_character(a$poly_id), a$zz)
    expect_identical(d$mean, as.numeric(tapply(a$zz[o], factor(a$poly_id[o],
                                                               levels = cells$poly_id),
                                               mean)), info = inf)
    shift_var <- d$plain_var - d$kr_var
    expect_identical(d$kr_shift,
                     ifelse(is.finite(shift_var) & shift_var > 1e-4 * d$plain_var,
                            (d$kr_pred - d$mean) / sqrt(pmax(shift_var, 0)), NA_real_),
                     info = inf)
    expect_true(sum(is.finite(d$kr_shift)) > 30L, info = inf)
  }
})
