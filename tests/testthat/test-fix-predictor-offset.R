# A predictor at a large offset is fitted, not dropped.
# The QR decomposition behind lm.fit() and lm() drops a column whose part
# orthogonal to the earlier columns is under 1e-7 of its norm, so a predictor
# whose spread is under about 1e-7 of its size (1e8 plus a spread of 1) was
# taken as collinear with the intercept and left out, in silence, of the
# profile's OLS fit, of the cell-mean regression of .morans_i_for_k() and of
# estimate_sac_range()'s OLS detrending.  rss, cp, Moran's I and z, the
# nugget and the range then read the residuals of the intercept alone, and
# moved with the predictor's origin.

# A smooth field, a predictor w and noise on 600 points; w as stored at 1e15,
# so that 1e8 + wq and 1e15 + wq hold the same spread.
.po_layer <- function(n = 600, seed = 3) {
  withr::with_seed(seed, {
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4)
    w <- rnorm(n)
    z <- 2 * w + sin(x / 2000) + cos(y / 1500) + rnorm(n, sd = 0.5)
  })
  wq <- (1e15 + w) - 1e15
  list(xy = cbind(x, y), wq = wq, z = z,
       sf = sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"),
                         crs = 32615))
}

# 2,000 points in 40 clusters with structure shorter than the clusters.
.po_clustered <- function(seed = 5) {
  withr::with_seed(seed, {
    nc <- 40; per <- 50
    cx <- runif(nc, 0, 1e4); cy <- runif(nc, 0, 1e4)
    x <- rep(cx, each = per) + rnorm(nc * per, sd = 150)
    y <- rep(cy, each = per) + rnorm(nc * per, sd = 150)
    w <- rnorm(nc * per)
    z <- 2 * w + sin(x / 400) + cos(y / 300) + rnorm(nc * per, sd = 0.3)
  })
  wq <- (1e15 + w) - 1e15
  list(wq = wq,
       sf = sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"),
                         crs = 32615))
}

.po_profile <- function(d, sac) {
  suppressWarnings(
    resolution_profile(d, "z", predictor_vars = "p", n_levels = 4, seed = 1, sac = sac))
}

.po_cols <- c("rss", "cp", "moran_i", "moran_z")

test_that("lm.fit() drops the predictor at 1e8 on this layer, so the fits below are reached", {
  lay <- .po_layer()
  expect_identical(stats::lm.fit(cbind(1, lay$wq), lay$z)$rank, 2L)
  expect_identical(stats::lm.fit(cbind(1, 1e8 + lay$wq), lay$z)$rank, 1L)
  expect_identical(stats::lm.fit(cbind(1, 1e15 + lay$wq), lay$z)$rank, 1L)
})

test_that("resolution_profile() gives one rss, cp, moran_i and moran_z at every offset of a predictor", {
  lay <- .po_layer()
  for (sac in list(2000, NULL)) {
    d <- lay$sf; d$p <- lay$wq
    ref <- .po_profile(d, sac)
    expect_identical(attr(ref, "variable"), "residuals")
    expect_true(all(is.finite(ref$rss)) && all(is.finite(ref$cp)))
    expect_true(any(is.finite(ref$moran_z)))
    for (off in c(1e4, 1e6, 1e7, 1e8, 1e10, 1e15)) {
      d$p <- off + lay$wq
      got <- .po_profile(d, sac)
      info <- sprintf("sac = %s, offset %g", if (is.null(sac)) "NULL" else format(sac), off)
      expect_identical(got$levels, ref$levels, info = info)
      for (cl in .po_cols)
        expect_equal(got[[cl]], ref[[cl]], tolerance = 1e-8, info = paste(info, cl))
    }
  }
})

test_that(".morans_i_for_k() regresses on a predictor at a large offset, loop and grouped", {
  lay <- .po_layer()
  for (k in c(20L, 40L)) {
    cl <- withr::with_seed(1, stats::kmeans(lay$xy, k, nstart = 1, iter.max = 50)$cluster)
    for (grouped in c(FALSE, TRUE)) {
      ref <- .morans_i_for_k(lay$xy, lay$z, matrix(lay$wq, ncol = 1L), cl, grouped = grouped)
      expect_true(all(is.finite(ref)))
      # The intercept alone gives another I: the predictor is read.
      alone <- .morans_i_for_k(lay$xy, lay$z, matrix(0, nrow(lay$xy), 0L), cl,
                               grouped = grouped)
      expect_gt(abs(alone[["I"]] - ref[["I"]]), 1e-3)
      for (off in c(1e4, 1e7, 1e8, 1e10, 1e15)) {
        got <- .morans_i_for_k(lay$xy, lay$z, matrix(off + lay$wq, ncol = 1L), cl,
                               grouped = grouped)
        expect_equal(got, ref, tolerance = 1e-8,
                     info = sprintf("k = %d, grouped = %s, offset %g", k, grouped, off))
      }
      # A second predictor at an offset of its own, beside one at 0, its
      # spread as stored at 1e12.
      p2 <- cbind(lay$wq, ((1e12 + lay$xy[, 1] / 1e4) - 1e12))
      ref2 <- .morans_i_for_k(lay$xy, lay$z, p2, cl, grouped = grouped)
      got2 <- .morans_i_for_k(lay$xy, lay$z, cbind(p2[, 1], 1e12 + p2[, 2]), cl,
                              grouped = grouped)
      expect_equal(got2, ref2, tolerance = 1e-8,
                   info = sprintf("k = %d, grouped = %s, two predictors", k, grouped))
    }
  }
})

test_that("estimate_sac_range() detrends on a predictor at a large offset", {
  skip_if_not_installed("gstat")
  for (lay in list(.po_layer(), .po_clustered())) {
    d <- lay$sf; d$p <- lay$wq
    ref <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
    expect_true(is.finite(attr(ref, "nugget")))
    for (off in c(1e8, 1e15)) {
      d$p <- off + lay$wq
      got <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
      expect_identical(is.na(as.numeric(got)), is.na(as.numeric(ref)), info = format(off))
      expect_equal(as.numeric(got), as.numeric(ref), tolerance = 1e-8, info = format(off))
      expect_equal(attr(got, "nugget"), attr(ref, "nugget"), tolerance = 1e-8,
                   info = format(off))
      expect_identical(attr(got, "rejected_reason"), attr(ref, "rejected_reason"),
                       info = format(off))
      expect_identical(attr(got, "detrend_method"), "ols")
    }
  }
})

test_that("determine_optimal_levels(criterion = \"morans_i\") gives one ranking and one moran_z at every offset of a predictor", {
  lay <- .po_clustered()
  run <- function(p) {
    d <- lay$sf; d$p <- p
    suppressWarnings(determine_optimal_levels(
      d, response_var = "z", predictor_vars = "p", criterion = "morans_i", max_levels = 60))
  }
  ref <- run(lay$wq)
  dg <- attr(ref, "diagnostics")
  expect_false(is.null(dg))
  expect_true(all(is.finite(dg$moran_z[dg$eval_ks])))
  for (off in c(1e8, 1e15)) {
    got <- run(off + lay$wq)
    gd <- attr(got, "diagnostics")
    expect_identical(as.integer(got), as.integer(ref), info = format(off))
    expect_identical(gd$eval_ks, dg$eval_ks, info = format(off))
    expect_equal(gd$moran_z, dg$moran_z, tolerance = 1e-8, info = format(off))
    expect_equal(gd$moran_i, dg$moran_i, tolerance = 1e-8, info = format(off))
  }
})

test_that("a constant predictor, at any value, gives the fit of the intercept alone", {
  lay <- .po_layer()
  d <- lay$sf
  base <- suppressWarnings(
    resolution_profile(d, "z", n_levels = 4, seed = 1, sac = 2000))
  cl <- withr::with_seed(1, stats::kmeans(lay$xy, 40L, nstart = 1, iter.max = 50)$cluster)
  alone <- .morans_i_for_k(lay$xy, lay$z, matrix(0, nrow(lay$xy), 0L), cl)
  d$p <- 0
  sac0 <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
  for (v in c(0, 5, 1e8, 1e15)) {
    d$p <- v
    got <- .po_profile(d, 2000)
    for (cn in c("rss", "moran_i", "moran_z"))
      expect_equal(got[[cn]], base[[cn]], tolerance = 1e-8, info = paste(v, cn))
    expect_equal(.morans_i_for_k(lay$xy, lay$z, matrix(v, nrow(lay$xy), 1L), cl), alone,
                 tolerance = 1e-8, info = format(v))
    e <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
    expect_equal(attr(e, "nugget"), attr(sac0, "nugget"), tolerance = 1e-8, info = format(v))
    expect_identical(attr(e, "rejected_reason"), attr(sac0, "rejected_reason"),
                     info = format(v))
  }
})

test_that("a predictor collinear in fact with another, at an offset, is still dropped", {
  lay <- .po_layer()
  cl <- withr::with_seed(1, stats::kmeans(lay$xy, 40L, nstart = 1, iter.max = 50)$cluster)
  one <- .morans_i_for_k(lay$xy, lay$z, matrix(lay$wq, ncol = 1L), cl)
  # Both at an offset, so the raw fit drops both and the fit about the
  # centres keeps one.
  two <- .morans_i_for_k(lay$xy, lay$z, cbind(1e8 + lay$wq, 3e8 + 2 * lay$wq), cl)
  expect_equal(two, one, tolerance = 1e-8)
  d <- lay$sf; d$p <- lay$wq
  ref <- suppressWarnings(
    resolution_profile(d, "z", predictor_vars = "p", n_levels = 4, seed = 1, sac = 2000))
  d$p <- 1e8 + lay$wq
  d$q <- 3e8 + 2 * lay$wq
  got <- suppressWarnings(
    resolution_profile(d, "z", predictor_vars = c("p", "q"), n_levels = 4, seed = 1,
                       sac = 2000))
  for (cn in .po_cols) expect_equal(got[[cn]], ref[[cn]], tolerance = 1e-8, info = cn)
  cly <- .po_clustered()
  d <- cly$sf; d$p <- cly$wq
  e1 <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
  d$p <- 1e8 + cly$wq
  d$q <- 3e8 + 2 * cly$wq
  e2 <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = c("p", "q")))
  expect_equal(attr(e2, "nugget"), attr(e1, "nugget"), tolerance = 1e-8)
  expect_identical(attr(e2, "rejected_reason"), attr(e1, "rejected_reason"))
})

# A little below the offset where lm.fit() drops the predictor, the raw fit
# keeps it but works on the rounding of the offset: with a predictor that
# explains most of the response, rss, cp, Moran's I and z and the nugget
# moved with the origin by up to 1e-4 of their size.
.po_strong <- function(seed = 3, n = 600) {
  withr::with_seed(seed, {
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4); w <- rnorm(n)
  })
  wq <- (1e15 + w) - 1e15
  z <- withr::with_seed(seed + 1L,
    2e4 * wq + sin(x / 2000) + cos(y / 1500) + rnorm(n, sd = 0.5))
  list(xy = cbind(x, y), wq = wq, z = z,
       sf = sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"),
                         crs = 32615))
}

test_that("a predictor the raw fit keeps at a large offset gives the numbers of the predictor at 0", {
  lay <- .po_strong()
  # The raw fit keeps the predictor at these offsets, so the gate is not
  # the dropped column.
  for (off in c(1e6, 3e6, 9.9e6))
    expect_identical(stats::lm.fit(cbind(1, off + lay$wq), lay$z)$rank, 2L)
  d <- lay$sf; d$p <- lay$wq
  ref <- .po_profile(d, 2000)
  ref_e <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
  cl <- withr::with_seed(1, stats::kmeans(lay$xy, 40L, nstart = 1, iter.max = 50)$cluster)
  ref_m <- .morans_i_for_k(lay$xy, lay$z, matrix(lay$wq, ncol = 1L), cl)
  expect_true(all(is.finite(ref_m)))
  for (off in c(1e6, 3e6, 9.9e6)) {
    d$p <- off + lay$wq
    got <- .po_profile(d, 2000)
    for (cn in .po_cols)
      expect_equal(got[[cn]], ref[[cn]], tolerance = 1e-8, info = paste(off, cn))
    e <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
    expect_equal(attr(e, "nugget"), attr(ref_e, "nugget"), tolerance = 1e-8, info = format(off))
    expect_identical(attr(e, "rejected_reason"), attr(ref_e, "rejected_reason"),
                     info = format(off))
    expect_equal(.morans_i_for_k(lay$xy, lay$z, matrix(off + lay$wq, ncol = 1L), cl), ref_m,
                 tolerance = 1e-8, info = format(off))
  }
})

test_that("a predictor near its centre is fitted on its raw values, to the last bit", {
  lay <- .po_layer()
  # The frozen predictors' shape: a centre a few times the spread.
  P <- cbind(lay$wq, 5 + 1.3 * lay$wq, 2020 + round(lay$wq))
  expect_identical(.centre_columns(P), P)
  X <- cbind(1, P)
  kept <- .lm_fit_centred(X, lay$z)
  expect_identical(kept$X, X)
  expect_identical(kept$fit, stats::lm.fit(X, lay$z))
  df <- data.frame(a = P[, 1], b = P[, 2], f = factor(rep(c("u", "v"), 300)))
  expect_identical(.centre_columns(df), df)
  # A constant column is left as it is, at any value.
  expect_identical(.centre_columns(cbind(rep(1e8, 5))), cbind(rep(1e8, 5)))
})

test_that("estimate_sac_range() detrends on a timestamp in seconds as on the same numbers", {
  skip_if_not_installed("gstat")
  lay <- .po_clustered()
  d <- lay$sf
  # Tenths of a second, as stored at 1.7e9.
  secs <- (1.7e9 + round(lay$wq * 10) / 10) - 1.7e9
  d$p <- secs
  ref <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
  for (p in list(1.7e9 + secs,
                 as.POSIXct(1.7e9 + secs, origin = "1970-01-01", tz = "UTC"),
                 as.difftime(1.7e9 + secs, units = "secs"))) {
    d$p <- p
    got <- suppressWarnings(estimate_sac_range(d, "z", predictor_vars = "p"))
    info <- class(p)[1]
    expect_identical(attr(got, "detrend_method"), "ols", info = info)
    expect_identical(is.na(as.numeric(got)), is.na(as.numeric(ref)), info = info)
    expect_equal(as.numeric(got), as.numeric(ref), tolerance = 1e-8, info = info)
    expect_equal(attr(got, "nugget"), attr(ref, "nugget"), tolerance = 1e-8, info = info)
    expect_identical(attr(got, "rejected_reason"), attr(ref, "rejected_reason"), info = info)
  }
})

test_that("an exact fit on a northing is still constant where the northing is fitted about its centre", {
  skip_if_not_installed("gstat")
  # The northing 4.5e6 plus a span of 4,000 is centred before the fits, and
  # a response computed from it carries the rounding of 4.317 times 4.5e6: the
  # test for an exact fit reads the terms of the predictor as given
  # (.rounding_design()), as it did before the centring.
  lay <- .po_layer()
  N <- 4.5e6 + lay$xy[, 2] * 0.4
  expect_false(is.null(attr(.centre_columns(cbind(N)), "centres")))
  v <- 1.7 - 4.317 * N + 4.317 * mean(N)
  cl <- withr::with_seed(1, stats::kmeans(lay$xy, 40L, nstart = 1, iter.max = 50)$cluster)
  for (grouped in c(FALSE, TRUE)) {
    m <- .morans_i_for_k(lay$xy, v, cbind(N), cl, grouped = grouped)
    expect_true(all(is.na(m)), info = format(grouped))
  }
  d <- lay$sf; d$N <- N; d$v <- v
  s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = "N"))
  expect_identical(attr(s, "rejected_reason"),
                   "the residuals on predictor_vars are constant: the predictors explain the response exactly")
  r <- suppressWarnings(resolution_profile(d, "v", predictor_vars = "N", n_levels = 3,
                                           nstart = 1, seed = 1))
  expect_true(attr(r, "one_cell")$constant)
  expect_true(all(is.na(r$cp)))
  # Beside it, a response the northing does not fit exactly is scored.
  d$v <- v + sin(lay$xy[, 1] / 2000)
  r <- suppressWarnings(resolution_profile(d, "v", predictor_vars = "N", n_levels = 3,
                                           nstart = 1, seed = 1, sac = 2000))
  expect_true(any(is.finite(r$cp)))
})
