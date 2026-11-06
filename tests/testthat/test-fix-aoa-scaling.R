# tests/testthat/test-fix-aoa-scaling.R
# ---------------------------------------------------------------------------
# The variance test and the dropped-column slack of area_of_applicability()
# (BUG-88, BUG-97).
#
# .aoa_scaling() kept a predictor only where sd > sqrt(eps) * max(abs(x)), a
# bound on the column's distance from 0, not on its spread: a predictor with
# sd 1 at 1e8 (a timestamp, a coordinate-like column) was dropped as having
# "no variance", and its extrapolated rows were marked Inf or read as inside.
# A column the bound drops is now kept where .constant_column() finds
# spread.  A column still dropped is compared with new rows on 3 units in
# the last place of its values, where the slack was 1.5e-8 of their size
# (25 units for a constant timestamp at 1.7e9).
# ---------------------------------------------------------------------------

aoa_offset_layout <- function() {
  set.seed(1)
  n <- 120
  list(a = stats::rnorm(n), b = stats::rnorm(n),
       an = c(stats::rnorm(20), stats::rnorm(20, 4)), bn = stats::rnorm(40))
}

aoa_run <- function(tr, nw, vars) {
  res <- NULL
  log <- capture_spatialkit_log(
    res <- suppressWarnings(area_of_applicability(
      nw, train_sf = tr, predictor_vars = vars, weights = NULL, coords = "no")),
    level = logger::WARN)
  res$log <- log
  res
}


test_that("a predictor at an offset is kept and gives the index it gives at 0", {
  L <- aoa_offset_layout()
  for (o in c(1e8, 1.7e9, 1e15)) {
    # The control is the same stored spread at offset 0.
    tr  <- data.frame(a = o + L$a, b = L$b)
    nw  <- data.frame(a = o + L$an, b = L$bn)
    trc <- data.frame(a = (o + L$a) - o, b = L$b)
    nwc <- data.frame(a = (o + L$an) - o, b = L$bn)
    r  <- aoa_run(tr, nw, c("a", "b"))
    rc <- aoa_run(trc, nwc, c("a", "b"))
    expect_identical(r$dropped_vars, character(0))
    expect_identical(r$predictor_vars, c("a", "b"))
    expect_false(any(grepl("no variance", r$log)))
    expect_equal(r$aoa$DI, rc$aoa$DI, tolerance = 1e-8)
    expect_equal(r$threshold, rc$threshold, tolerance = 1e-8)
    expect_identical(r$aoa$AOA, rc$aoa$AOA)
    expect_identical(r$n_outside, rc$n_outside)
    # The rows 4 SD out on a are outside, as at 0.
    expect_identical(sum(!r$aoa$AOA[21:40]), sum(!rc$aoa$AOA[21:40]))
    expect_gte(sum(!r$aoa$AOA[21:40]), 15L)
    # The scaling records the training mean and sd.
    expect_equal(unname(r$scaling$center[["a"]]), mean(o + L$a), tolerance = 1e-12)
    expect_equal(unname(r$scaling$scale[["a"]]), stats::sd((o + L$a) - o),
                 tolerance = 1e-8)
  }
})

test_that("an extrapolation along a predictor at an offset is outside (BUG-97)", {
  set.seed(3); w <- stats::rnorm(300); wq <- (1e15 + w) - 1e15; q <- stats::rnorm(300)
  set.seed(9); nd_w <- c(stats::rnorm(50), stats::rnorm(20, 6)); nd_wq <- (1e15 + nd_w) - 1e15
  set.seed(10); nd_q <- stats::rnorm(70)
  r0 <- aoa_run(data.frame(p = wq, q = q), data.frame(p = nd_wq, q = nd_q), c("p", "q"))
  expect_identical(sum(!r0$aoa$AOA[51:70]), 20L)
  for (o in c(1e8, 1e15)) {
    r <- aoa_run(data.frame(p = o + wq, q = q), data.frame(p = o + nd_wq, q = nd_q),
                 c("p", "q"))
    expect_identical(r$dropped_vars, character(0))
    expect_false(any(grepl("no variance|never has", r$log)))
    expect_true(all(is.finite(r$aoa$DI)))
    expect_equal(r$aoa$DI, r0$aoa$DI, tolerance = 1e-8)
    expect_equal(r$threshold, r0$threshold, tolerance = 1e-8)
    expect_identical(r$aoa$AOA, r0$aoa$AOA)
  }
})

test_that(".aoa_scaling() keeps a column with spread and drops one without", {
  set.seed(4)
  z <- stats::rnorm(50)
  for (o in c(1e8, 1.7e9, 1e15)) {
    sc <- .aoa_scaling(cbind(a = o + z, b = z))
    expect_identical(unname(sc$keep), c(TRUE, TRUE))
  }
  # Units whose squares vanish or overflow (sd() gives 0 and Inf): still
  # spread, with the standard deviation of the column.
  for (s in c(1e-170, 1e200)) {
    sc <- .aoa_scaling(cbind(a = z * s, b = z))
    expect_true(sc$keep[["a"]])
    expect_equal(sc$scale[["a"]] / s, stats::sd(z), tolerance = 1e-12)
    expect_equal(sc$center[["a"]] / s, mean(z), tolerance = 1e-12)
  }
  # Constant columns, at 0, at an offset, and differing by a unit in the last
  # place of 1e8 (rounding), are dropped.
  u <- 2^26 * .Machine$double.eps
  sc <- .aoa_scaling(cbind(z = rep(0, 50), c8 = rep(1e8, 50),
                           r8 = 1e8 + rep(c(0, u), 25), a = z))
  expect_identical(unname(sc$keep), c(FALSE, FALSE, FALSE, TRUE))
  # A column the first bound keeps has the mean and sd it had.
  X <- cbind(a = z, b = 5 + 2 * z)
  sc <- .aoa_scaling(X)
  expect_identical(sc$center, colMeans(X))
  expect_identical(sc$scale, apply(X, 2L, stats::sd))
})

test_that("a dropped column reads new rows on the rounding of its values (BUG-88)", {
  L <- aoa_offset_layout()
  # A constant timestamp at 1.7e9; 20 new rows 25 seconds away.
  tr <- data.frame(a = rep(1.7e9, 120), b = L$b)
  nw <- data.frame(a = 1.7e9 + c(rep(0, 20), rep(25, 20)), b = L$bn)
  r <- aoa_run(tr, nw, c("a", "b"))
  expect_identical(r$dropped_vars, "a")
  expect_true(any(grepl("no variance in the training data: a", r$log)))
  expect_identical(r$aoa$DI[21:40], rep(Inf, 20))
  expect_true(all(is.finite(r$aoa$DI[1:20])))
  # Rounding of the training value is not a different value: a new row one
  # unit in the last place away is judged on b.
  u <- 2^30 * .Machine$double.eps
  nw2 <- data.frame(a = 1.7e9 + c(u, -u, 2 * u), b = L$bn[1:3])
  r2 <- aoa_run(tr, nw2, c("a", "b"))
  expect_true(all(is.finite(r2$aoa$DI)))
  # A column of zeros reads rounding of a value of 1 as 0, and a dummy of 1 as
  # different.
  tr0 <- data.frame(a = rep(0, 120), b = L$b)
  nw0 <- data.frame(a = c(1e-17, -1e-17, 1, 1e-9), b = L$bn[1:4])
  r0 <- aoa_run(tr0, nw0, c("a", "b"))
  expect_true(all(is.finite(r0$aoa$DI[1:2])))
  expect_identical(r0$aoa$DI[3:4], c(Inf, Inf))
})

test_that("the slack of a dropped column is 3 units in the last place", {
  L <- aoa_offset_layout()
  # At 1.7e9 a unit in the last place is 2^-22: 3 of them on either side are
  # rounding, 4 and 5 are a different value.
  u  <- 2^30 * .Machine$double.eps
  tr <- data.frame(a = rep(1.7e9, 120), b = L$b)
  nw <- data.frame(a = 1.7e9 + c(3, -3, 4, -4, 5, -5) * u, b = L$bn[1:6])
  r  <- aoa_run(tr, nw, c("a", "b"))
  expect_identical(r$dropped_vars, "a")
  expect_true(all(is.finite(r$aoa$DI[1:2])))
  expect_identical(r$aoa$DI[3:6], rep(Inf, 4))
  # A column of zeros: 3 units in the last place of 1 are 6.7e-16.
  tr0 <- data.frame(a = rep(0, 120), b = L$b)
  nw0 <- data.frame(a = c(6e-16, -6e-16, 1e-15, -1e-15), b = L$bn[1:4])
  r0  <- aoa_run(tr0, nw0, c("a", "b"))
  expect_true(all(is.finite(r0$aoa$DI[1:2])))
  expect_identical(r0$aoa$DI[3:4], c(Inf, Inf))
})

test_that("a dropped column at a subnormal value is read on its own rounding", {
  L <- aoa_offset_layout()
  # .ulp_of() underflows to 0 below about 1e-292; the rounding of such a value
  # is the smallest positive double, not the rounding of 1.  The same column
  # at 1 and at 1e-300 marks rows at 0, 2c and -c outside.
  tiny <- 2^-1074
  for (c0 in c(1, 1e-300, 1e-310)) {
    tr <- data.frame(a = rep(c0, 120), b = L$b)
    nw <- data.frame(a = c(0, 2, -1) * c0, b = L$bn[1:3])
    r  <- aoa_run(tr, nw, c("a", "b"))
    expect_identical(r$dropped_vars, "a")
    expect_identical(r$aoa$DI, rep(Inf, 3))
  }
  c0 <- 1e-310
  tr <- data.frame(a = rep(c0, 120), b = L$b)
  nw <- data.frame(a = c0 + c(0, 3, -3, 4, -4) * tiny, b = L$bn[1:5])
  r  <- aoa_run(tr, nw, c("a", "b"))
  expect_true(all(is.finite(r$aoa$DI[1:3])))
  expect_identical(r$aoa$DI[4:5], rep(Inf, 2))
})
