# tests/testthat/test-core-pool-reference.R
# ---------------------------------------------------------------------------
# The engine of R/pool.R against the reference functions of
# helper-pool-reference.R, bit for bit.
#
# On the thirteen layouts of the helper (9 to 3,000 rows, 2 to 1,500 regions,
# 18 to 2 rows per region, uniform and clustered, up to 536 regions of one
# row) every number the engine and the reference both compute is held to
# identical(), with FNN and with the exact searches of base R.  The layouts
# have no two candidate neighbours at one distance and no repeated location;
# test-core-pool-properties.R covers those.
#
# The reference functions call FNN, so every test here but the fifth skips
# without it.
# ---------------------------------------------------------------------------

# One layout with what the reference returns for it, computed once.
.pr_store <- new.env(parent = emptyenv())
.pr_case <- function(i) {
  key <- as.character(i)
  if (is.null(.pr_store[[key]])) {
    L   <- do.call(pool_layout, POOL_LAYOUTS[[i]])
    nn1 <- pool_nearest_other(L$xy)
    L$nn1 <- nn1
    L$R   <- region_estimates(L$y, L$gi, L$K, L$xy, want = c("P", "L"))
    L$Fr  <- fit_regions(L$gi, L$gi[1:3], L$y, L$xy, L$xy[1:3, , drop = FALSE], nn1)
    L$tag <- sprintf("layout %d: %d rows in %d regions", i, length(L$y), L$K)
    .pr_store[[key]] <- L
  }
  .pr_store[[key]]
}

# The two searches: FNN's kd-tree and the exact one.
.pr_paths <- c(TRUE, FALSE)
.pr_tag <- function(L, use_fnn) paste0(L$tag, if (use_fnn) ", with FNN" else ", exact search")


test_that("1. the pooled value and the plain mean are region_estimates()'s", {
  skip_if_not_installed("FNN")
  for (i in seq_along(POOL_LAYOUTS)) for (f in .pr_paths) {
    L <- .pr_case(i); tag <- .pr_tag(L, f)
    A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f)
    expect_identical(A$value, unname(L$R$est[, "L"]), info = tag)
    expect_identical(A$m, unname(L$R$est[, "P"]), info = tag)
    # The one-way parts the estimates are built from.
    expect_identical(A$n, L$R$a$n, info = tag)
    expect_identical(A$S, L$R$a$S, info = tag)
    expect_identical(A$N, L$R$a$N, info = tag)
    expect_identical(A$K, L$K, info = tag)
    expect_identical(A$status, "ok", info = tag)
  }
})

test_that("2. the parts of the pooled value and both residual vectors are fit_regions()'s", {
  skip_if_not_installed("FNN")
  for (i in seq_along(POOL_LAYOUTS)) for (f in .pr_paths) {
    L <- .pr_case(i); tag <- .pr_tag(L, f); Fr <- L$Fr
    A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f)
    expect_identical(A$sigma2, Fr$sigma2, info = tag)
    expect_identical(A$tau2, Fr$tau2, info = tag)
    expect_identical(A$w, Fr$w, info = tag)
    expect_identical(A$target, Fr$t, info = tag)
    expect_identical(A$Q, Fr$Q, info = tag)
    expect_identical(A$nb, Fr$nb, info = tag)
    expect_identical(A$resid_pooled, Fr$eL, info = tag)
    expect_identical(A$resid_mean, Fr$eP, info = tag)
    expect_identical(A$loo_pooled, Fr$loo_L, info = tag)
    expect_identical(A$loo_pooled_se, Fr$loo_L_se, info = tag)
    # The centroids are summed in coordinate order where the reference sums
    # them in row order, so they agree to rounding and not to the last bit.
    expect_equal(unname(A$cen), unname(rowsum(L$xy, L$gi) / A$n), tolerance = 1e-12,
                 info = tag)
  }
})

test_that("3. the plain mean's leave-one-out error is loo_plain()'s and fit_regions()'s", {
  skip_if_not_installed("FNN")
  for (i in seq_along(POOL_LAYOUTS)) for (f in .pr_paths) {
    L <- .pr_case(i); tag <- .pr_tag(L, f)
    A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f)
    expect_identical(A$loo_mean, loo_plain(L$y, L$gi, L$R$a, L$nn1), info = tag)
    expect_identical(A$loo_mean, L$Fr$loo_P, info = tag)
    expect_identical(A$loo_mean_se, stats::sd(L$Fr$eP^2) / sqrt(length(L$y)), info = tag)
    # The nearest other row the fit found for the rows alone in their region
    # is the reference's, and a fit handed the reference's rows is the same.
    alone <- which(A$n[L$gi] == 1L)
    expect_identical(.nearest_other(L$xy, L$gi, alone, f), L$nn1[alone], info = tag)
    B <- .pool_fit(L$y, L$gi, L$K, L$xy, nn1 = L$nn1, use_fnn = f)
    expect_identical(B, A, info = tag)
  }
  # The layouts hold regions of one row, which is where the two differ from
  # the residual formula of the regions that hold more.
  expect_gt(sum(.pr_case(7L)$R$a$n == 1L), 500L)
})

test_that("4. the interval's offsets are the 5th and 95th percentiles, and the ends loo_int()'s", {
  skip_if_not_installed("FNN")
  for (i in seq_along(POOL_LAYOUTS)) for (f in .pr_paths) {
    L <- .pr_case(i); tag <- .pr_tag(L, f)
    A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f)
    e <- A$resid_pooled
    off <- .obs_offsets(e, 0.90)
    expect_identical(off, stats::quantile(e, c(0.05, 0.95), names = FALSE, type = 7),
                     info = tag)
    expect_identical(cbind(A$value + off[1L], A$value + off[2L]),
                     loo_int(A$value, e, 1, 1), info = tag)
    # And against the reference's own value and residuals.
    expect_identical(cbind(A$value + off[1L], A$value + off[2L]),
                     loo_int(L$Fr$est[, "L"], L$Fr$eL, 1, 1), info = tag)
  }
})

test_that("5. the count at a size is rung_k()'s", {
  # rung_k() is arithmetic alone, so this one runs without FNN.
  expect_identical(.OBS_PER_CELL_LADDER, c(18, 12, 9, 6, 4.5, 3, 2))
  expect_identical(.OBS_PER_CELL_DEFAULT, 4.5)
  n <- 9:10000
  halves <- 0L
  for (m in .OBS_PER_CELL_LADDER) {
    want <- vapply(n, function(v) rung_k(v, m), integer(1))
    expect_identical(.count_at(n, m), want, info = sprintf("size %s", m))
    expect_identical(vapply(n, function(v) .count_at(v, m), integer(1)), want,
                     info = sprintf("size %s, one count at a time", m))
    halves <- halves + sum((n / m) %% 1 == 0.5)
  }
  # The range holds quotients that end in a half, which R rounds to the even
  # number: 4.5 to 4, 13.5 to 14, 2.5 to 2, 3.5 to 4.
  expect_gt(halves, 1000L)
  expect_identical(.count_at(9, 2), 4L)
  expect_identical(.count_at(27, 2), 14L)
  expect_identical(.count_at(45, 18), 2L)
  expect_identical(.count_at(63, 18), 4L)
  # Never fewer than two regions.
  expect_identical(.count_at(12, 18), 2L)
  expect_identical(.count_at(1, 4.5), 2L)
})

test_that("6. with one value beyond one per region the fit is region_estimates()'s", {
  skip_if_not_installed("FNN")
  # Forty regions of one row and a forty-first row in region 17: N - K = 1,
  # the fewest values sigma2 can rest on.
  cleanup <- pool_seed(61)
  xy <- cbind(stats::runif(40, 0, 1e5), stats::runif(40, 0, 1e5))
  xy <- rbind(xy, xy[17, ] + c(350, -210))
  y  <- sin(xy[, 1] / 2e4) + cos(xy[, 2] / 1.5e4) + stats::rnorm(41, 0, 0.6)
  cleanup()
  gi <- c(1:40, 17L); K <- 40L
  nn1 <- pool_nearest_other(xy)
  R  <- region_estimates(y, gi, K, xy, want = c("P", "L"))
  Fr <- fit_regions(gi, gi[1:3], y, xy, xy[1:3, , drop = FALSE], nn1)
  expect_identical(R$a$N - R$a$K, 1L)
  for (f in .pr_paths) {
    tag <- if (f) "with FNN" else "exact search"
    A <- .pool_fit(y, gi, K, xy, use_fnn = f)
    expect_identical(A$value, unname(R$est[, "L"]), info = tag)
    expect_identical(A$m, unname(R$est[, "P"]), info = tag)
    expect_identical(A$sigma2, R$a$sigma2, info = tag)
    expect_identical(A$sigma2, sum((y - A$m[gi])^2) / 1, info = tag)
    expect_identical(A$tau2, Fr$tau2, info = tag)
    expect_identical(A$w, Fr$w, info = tag)
    expect_identical(A$target, Fr$t, info = tag)
    expect_identical(A$nb, Fr$nb, info = tag)
    expect_identical(A$resid_pooled, Fr$eL, info = tag)
    expect_identical(A$resid_mean, Fr$eP, info = tag)
    expect_identical(A$loo_pooled, Fr$loo_L, info = tag)
    expect_identical(A$loo_mean, Fr$loo_P, info = tag)
    expect_true(A$status %in% c("ok", "no signal between cells"), info = tag)
  }
  # One row fewer and no region holds two: the reference stops, and the
  # engine returns the status without a fit.
  expect_error(region_estimates(y[1:40], 1:40, 40L, xy[1:40, ], want = c("P", "L")),
               "no group holds two")
  expect_identical(.pool_fit(y[1:40], 1:40, 40L, xy[1:40, ])$status,
                   "no cell holds two values")
})
