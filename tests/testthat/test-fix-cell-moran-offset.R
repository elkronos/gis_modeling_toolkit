# tests/testthat/test-fix-cell-moran-offset.R
# ---------------------------------------------------------------------------
# Moran's I of the cell means tests its residuals against the rounding of
# cell means: .POOL_ULPS units in the last place of the raw response over
# the square root of each cell's rows, plus the rounding of the sums the
# means are formed from, on the residual degrees of freedom, with the bound
# on the regression's rounding taken on the cell design.  Rows whose OLS
# residuals are within the rounding of a row give NA first.
# At a large offset it scored the residuals of an exact fit, which are
# rounding, where the response's raw values lost their spread (1e12 and
# 1e15), and it refused cell means with real spread at 1e15 by the bound of
# a single row, so determine_optimal_levels(criterion = "morans_i") gave
# other counts at offsets 0 and 1e15 (BUG-65).
# ---------------------------------------------------------------------------

.cm_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# 800 points over 4 km with three standard normal columns.
.cm_layer <- function() {
  .cm_seed(5)
  n <- 800
  x <- stats::runif(n, 0, 4000); y <- stats::runif(n, 0, 4000)
  list(xy = cbind(x, y), w = stats::rnorm(n), w2 = stats::rnorm(n), f = stats::rnorm(n))
}

.cm_cells <- function(xy, k) {
  .cm_seed(k)
  stats::kmeans(xy, k, nstart = 1)$cluster
}

test_that("Moran's I of the cell means is NA on an exact fit at every offset", {
  lay <- .cm_layer()
  P <- cbind(w = lay$w, w2 = lay$w2)
  for (k in c(20L, 40L, 80L)) {
    cl <- .cm_cells(lay$xy, k)
    .cm_seed(100 + k)
    cf <- matrix(round(stats::rnorm(90, sd = 5), 2), ncol = 3)
    for (off in c(0, 1e6, 1e12, 1e15)) {
      scored <- 0L
      for (i in seq_len(nrow(cf))) {
        y <- off + cf[i, 1] + cf[i, 2] * lay$w + cf[i, 3] * lay$w2
        for (grouped in c(TRUE, FALSE)) {
          m <- .morans_i_for_k(lay$xy, y, P, cl, grouped = grouped)
          if (!is.na(m[["I"]]) || !is.na(m[["z"]])) scored <- scored + 1L
        }
      }
      expect_identical(scored, 0L, info = sprintf("k = %d, offset %g", k, off))
    }
  }
})

test_that("Moran's I of cell means with spread is the same at offsets 0 and 1e15", {
  lay <- .cm_layer()
  P <- matrix(lay$w, ncol = 1L)
  for (k in c(20L, 40L, 80L, 200L)) {
    cl <- .cm_cells(lay$xy, k)
    for (s in c(1, 3)) {
      v <- (1e15 + 1e5 * lay$w + s * lay$f) - 1e15
      for (grouped in c(TRUE, FALSE)) {
        ref <- .morans_i_for_k(lay$xy, v, P, cl, grouped = grouped)
        got <- .morans_i_for_k(lay$xy, 1e15 + v, P, cl, grouped = grouped)
        info <- sprintf("k = %d, s = %g, grouped = %s", k, s, grouped)
        expect_true(all(is.finite(ref)), info = info)
        expect_equal(got, ref, tolerance = 1e-8, info = info)
      }
    }
    # A spread of 0.3 is 2.4 units in the last place of 1e15 a row, which
    # the bound of a row calls rounding, and so does the bound of a cell.
    m <- .morans_i_for_k(lay$xy, 1e15 + 1e5 * lay$w + 0.3 * lay$f, P, cl, grouped = TRUE)
    expect_true(is.na(m[["I"]]) && is.na(m[["z"]]), info = sprintf("k = %d, s = 0.3", k))
  }
})

test_that("the cell-mean bound shrinks with the rows of each cell", {
  # 1e15 plus f on 40 cells of about 20 rows: the residuals are the centred
  # cell means, of standard deviation about 1 / sqrt(20), under the bound of
  # a row (three units of 0.125) and above that of a cell mean.
  lay <- .cm_layer()
  cl <- .cm_cells(lay$xy, 40L)
  none <- matrix(numeric(0), length(lay$f), 0L)
  v <- (1e15 + lay$f) - 1e15
  for (grouped in c(TRUE, FALSE)) {
    ref <- .morans_i_for_k(lay$xy, v, none, cl, grouped = grouped)
    expect_true(all(is.finite(ref)))
    expect_equal(.morans_i_for_k(lay$xy, 1e15 + v, none, cl, grouped = grouped), ref,
                 tolerance = 1e-8)
  }
  # The values' bound of .residuals_rounding() with the rows of each mean.
  cm <- as.numeric(tapply(v, cl, mean))
  nj <- as.numeric(table(cl))
  X <- matrix(1, length(cm), 1L)
  e <- cm - mean(cm)
  expect_false(.residuals_rounding(e, 1e15 + v, cm, X, mean(cm), n = nj))
  expect_true(.residuals_rounding(e * 0.05, 1e15 + v, cm, X, mean(cm), n = nj))
  # Without the rows it is the bound of a row, as before.
  expect_true(.residuals_rounding(e, 1e15 + v, cm, X, mean(cm)))
})

test_that("determine_optimal_levels(criterion = 'morans_i') gives the same counts at 0 and 1e15", {
  lay <- .cm_layer()
  cases <- list(
    exact = list(v = function(off) off - 5.98 - 4.76 * lay$w + 1.39 * lay$w2,
                 preds = c("w", "w2"), scored = FALSE),
    spread = list(v = function(off) off + ((1e15 + 1e5 * lay$w + lay$f) - 1e15),
                  preds = "w", scored = TRUE))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    res <- lapply(c(0, 1e15), function(off) {
      d <- data.frame(x = lay$xy[, 1], y = lay$xy[, 2], w = lay$w, w2 = lay$w2,
                      v = cs$v(off))
      d <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)
      suppressWarnings(determine_optimal_levels(d, response_var = "v",
                                                predictor_vars = cs$preds,
                                                criterion = "morans_i", max_levels = 40,
                                                on_no_elbow = "warn"))
    })
    expect_identical(as.integer(res[[1]]), as.integer(res[[2]]), info = nm)
    dg <- lapply(res, attr, "diagnostics")
    if (cs$scored) {
      z0 <- dg[[1]]$moran_z; z1 <- dg[[2]]$moran_z
      expect_true(any(is.finite(z0)), info = nm)
      expect_identical(is.finite(z0), is.finite(z1), info = nm)
      expect_equal(z1[is.finite(z1)], z0[is.finite(z0)], tolerance = 1e-8, info = nm)
    } else {
      expect_true(is.null(dg[[1]]) || !any(is.finite(dg[[1]]$moran_z)), info = nm)
      expect_true(is.null(dg[[2]]) || !any(is.finite(dg[[2]]$moran_z)), info = nm)
    }
  }
})

# 20,000 points over 4 km: ws the easting scaled to about -1.7 to 1.7, w
# and f standard normal.
.cm_big <- function() {
  .cm_seed(11)
  n <- 20000
  x <- stats::runif(n, 0, 4000); y <- stats::runif(n, 0, 4000)
  list(xy = cbind(x, y), ws = (x - 2000) / 1155, w = stats::rnorm(n), f = stats::rnorm(n))
}

.cm_big_cells <- function(xy, k) {
  .cm_seed(k)
  suppressWarnings(stats::kmeans(xy, k, nstart = 1, iter.max = 50)$cluster)
}

test_that("Moran's I of the cell means is NA on exact fits with many rows per cell", {
  # About 1,200 to 2,200 rows a cell.  The rounding of a cell mean's values
  # shrinks with its rows, but that of the sums it is formed from grows with
  # them, and a bound on the values alone scored exact fits at offset 0 on
  # the loop path.
  lay <- .cm_big()
  for (k in c(10L, 12L, 15L)) {
    cl <- .cm_big_cells(lay$xy, k)
    for (pr in c("w", "ws")) {
      P <- matrix(lay[[pr]], ncol = 1L)
      scored <- character(0)
      for (a in c(2, 7, 100, 1e4)) for (b in c(0.05, 0.3, 1, 5) * a / 100) {
        for (grouped in c(FALSE, TRUE)) {
          m <- .morans_i_for_k(lay$xy, a + b * lay[[pr]], P, cl, grouped = grouped)
          if (!is.na(m[["I"]]) || !is.na(m[["z"]]))
            scored <- c(scored, sprintf("a %g b %g grouped %s", a, b, grouped))
        }
      }
      expect_identical(scored, character(0), info = sprintf("k = %d on %s", k, pr))
    }
  }
})

test_that("Moran's I of the cell means is NA where the sums' rounding is the spread", {
  # 2 + 0.05 ws plus 30 units in the last place of 2 times f, on 1,200 to
  # 2,200 rows a cell: the cell means' residuals hold the spread over
  # sqrt(1,600), about 0.75 of a unit, beside the rounding of sums of 1,600
  # values, and z was scored up to 0.77 from that of the spread alone.  A
  # spread of 1e4 units is scored, with the I and z of the spread alone.
  lay <- .cm_big()
  u <- 2 * .Machine$double.eps
  P <- matrix(lay$ws, ncol = 1L)
  for (k in c(10L, 12L, 15L)) {
    cl <- .cm_big_cells(lay$xy, k)
    for (grouped in c(FALSE, TRUE)) {
      info <- sprintf("k = %d, grouped = %s", k, grouped)
      m <- .morans_i_for_k(lay$xy, 2 + 0.05 * lay$ws + 30 * u * lay$f, P, cl, grouped = grouped)
      expect_true(is.na(m[["I"]]) && is.na(m[["z"]]), info = info)
      ref <- .morans_i_for_k(lay$xy, 1e4 * u * lay$f, P, cl, grouped = grouped)
      got <- .morans_i_for_k(lay$xy, 2 + 0.05 * lay$ws + 1e4 * u * lay$f, P, cl,
                             grouped = grouped)
      expect_true(all(is.finite(ref)), info = info)
      expect_equal(got, ref, tolerance = 0.01, info = info)
    }
  }
})

test_that("Moran's I of the cell means is NA on an exact fit whose rounding a cell shares", {
  # 1e15 + 2 + b ws, b 0.05 and 0.1: the response changes by about a unit
  # in the last place of 1e15 (0.125) over the layer, so its values hold
  # three numbers and their rounding is the same over most rows of a cell.
  # The rows' residuals are within the rounding of a row, as at offset 0.
  lay <- .cm_big()
  P <- matrix(lay$ws, ncol = 1L)
  for (k in c(40L, 80L)) {
    cl <- .cm_big_cells(lay$xy, k)
    for (b in c(0.05, 0.1)) for (off in c(0, 1e15)) {
      m <- .morans_i_for_k(lay$xy, off + 2 + b * lay$ws, P, cl, grouped = TRUE)
      expect_true(is.na(m[["I"]]) && is.na(m[["z"]]),
                  info = sprintf("k = %d, b = %g, offset %g", k, b, off))
    }
  }
})

test_that("the cell-mean bound weights each cell by its rows and reads the design's rank", {
  # Values at 1e15, a unit in the last place 0.125: three units squared
  # are 0.1406.
  y_of <- function(yd) 1e15 + yd
  # Five cells of one row beside five of 400: the bound is the sum of
  # u2 / n_j on 9 degrees of freedom, 0.634, where u2 over the mean of the
  # rows would give 0.0063.
  .cm_seed(3)
  e <- stats::rnorm(10L)
  e <- (e - mean(e)) * sqrt(0.1 / sum((e - mean(e))^2))
  nj <- rep(c(1, 400), each = 5L)
  expect_true(.residuals_rounding(e, y_of(e), e, n = nj))
  expect_false(.residuals_rounding(e * sqrt(10), y_of(e), e * sqrt(10), n = nj))
  expect_false(.residuals_rounding(e, y_of(e), e, n = rep(200.5, 10L)))
  # The same with a design of two columns, on 8 degrees of freedom.
  X <- cbind(1, seq(-1, 1, length.out = 10L))
  ex <- stats::lm.fit(X, e)$residuals
  ex <- ex * sqrt(0.1 / sum(ex^2))
  yd <- drop(X %*% c(0.5, 1)) + ex
  expect_true(.residuals_rounding(ex, y_of(yd), yd, X, c(0.5, 1), n = nj))
  expect_false(.residuals_rounding(ex, y_of(yd), yd, X, c(0.5, 1), n = rep(200.5, 10L)))
  # Six means on a design of rank 4: 2 degrees of freedom, a bound of
  # 2 u2 = 0.28 with one row a mean, where 5 would give 0.70.
  .cm_seed(4)
  X4 <- cbind(1, matrix(stats::rnorm(18L), 6L))
  e4 <- stats::lm.fit(X4, stats::rnorm(6L))$residuals
  e4 <- e4 * sqrt(0.5 / sum(e4^2))
  yd4 <- drop(X4 %*% c(0.5, 1, -1, 0.25)) + e4
  expect_false(.residuals_rounding(e4, y_of(yd4), yd4, X4, c(0.5, 1, -1, 0.25), n = rep(1, 6L)))
  expect_true(.residuals_rounding(e4 * sqrt(0.5), y_of(yd4), yd4, X4, c(0.5, 1, -1, 0.25),
                                  n = rep(1, 6L)))
  # The rows' bound reads n - 1 degrees of freedom whatever the design.
  expect_true(.residuals_rounding(e4, y_of(yd4), yd4, X4, c(0.5, 1, -1, 0.25)))
  # The rounding of forming the means (`m`) adds to that of their values.
  e10 <- e * sqrt(10)
  expect_false(.residuals_rounding(e10, y_of(e10), e10, n = rep(100, 10L)))
  expect_true(.residuals_rounding(e10, y_of(e10), e10, n = rep(100, 10L), m = rep(0.5, 10L)))
})
