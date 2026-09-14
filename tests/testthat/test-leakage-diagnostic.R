# tests/testthat/test-leakage-diagnostic.R
# ---------------------------------------------------------------------------
# make_folds(block_kfold): the "blocks smaller than the autocorrelation range"
# warning must be reachable on the DEFAULT path.
#
# Before this file existed, `sac_range` inside make_folds() was assigned only
# under auto_range = TRUE, so the two leakage warnings -- both gated on
# is.finite(sac_range) -- could fire only in the one configuration that had
# already sized the blocks from the range and so never needed them.  On every
# default call the range stayed NA and the user learned nothing.  The fix
# estimates the range for the diagnostic alone whenever a response column is
# available (which it always is when a cv_*() function built the folds); it
# sizes nothing, so the folds are unchanged.
# ---------------------------------------------------------------------------

# A field with a short, well-identified range: a sum of sinusoids of period
# ~250 units on a 1000-unit extent.  estimate_sac_range() returns ~209 on it,
# which sits between the default block edge at k = 5 (1000 / 4 = 250, no
# leakage) and at k = 10 (1000 / 6 ~ 167, leakage).
leak_test_points <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  z <- sin(x / 40) + cos(y / 40) + rnorm(n, sd = 0.3)
  sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}

collect_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

test_that("the leakage warning fires on the default path when blocks are smaller than the range", {
  skip_if_not_installed("gstat")
  pts <- leak_test_points()
  r <- estimate_sac_range(pts, "z")
  expect_true(is.finite(r))
  expect_gt(r, 170); expect_lt(r, 250)   # the geometry below depends on this

  # k = 10 -> 30 target blocks -> 5 x 6 grid -> 167-unit edge < range.
  res <- collect_warnings(make_folds(pts, k = 10, method = "block_kfold",
                                     seed = 7, response_var = "z"))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "block dimension .* < autocorrelation range")
  expect_match(res$warnings, "auto_range = TRUE")

  # k = 5 -> 15 target blocks -> 4 x 4 grid -> 250-unit edge > range: silent.
  res5 <- collect_warnings(make_folds(pts, k = 5, method = "block_kfold",
                                      seed = 7, response_var = "z"))
  expect_length(res5$warnings, 0L)
})

test_that("the diagnostic sizes nothing: folds are identical with and without a response", {
  skip_if_not_installed("gstat")
  pts <- leak_test_points()
  with_resp <- suppressWarnings(
    make_folds(pts, k = 10, method = "block_kfold", seed = 7, response_var = "z"))
  without   <- make_folds(pts, k = 10, method = "block_kfold", seed = 7)
  expect_identical(with_resp$assignment$fold, without$assignment$fold)
  expect_identical(with_resp$params$grid_nx, without$params$grid_nx)
  expect_identical(with_resp$params$grid_ny, without$params$grid_ny)
  expect_identical(with_resp$params$block_size, without$params$block_size)
})

test_that("a hand-set block_size below the range warns too, without auto_range", {
  skip_if_not_installed("gstat")
  pts <- leak_test_points()
  res <- collect_warnings(make_folds(pts, k = 5, method = "block_kfold",
                                     seed = 7, response_var = "z",
                                     block_size = 100))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "block_size \\(100\\.0\\) < estimated autocorrelation range")
  # ... and a block_size above it does not.
  res2 <- collect_warnings(make_folds(pts, k = 5, method = "block_kfold",
                                      seed = 7, response_var = "z",
                                      block_size = 300))
  expect_length(res2$warnings, 0L)
})

test_that("without a response there is nothing to compare against, and no warning", {
  pts <- leak_test_points()
  res <- collect_warnings(make_folds(pts, k = 10, method = "block_kfold", seed = 7))
  expect_length(res$warnings, 0L)
})

test_that("fewer than 30 points skips the estimate rather than warning about it", {
  pts <- leak_test_points(n = 20)
  res <- collect_warnings(make_folds(pts, k = 3, method = "block_kfold",
                                     seed = 7, response_var = "z"))
  expect_length(res$warnings, 0L)
  expect_identical(spatialkit:::.sac_range_for_diagnostic(pts, "z", NULL, 1.0, 123L),
                   NA_real_)
})

test_that("the diagnostic estimate is unclassed and agrees with the public estimator", {
  skip_if_not_installed("gstat")
  pts <- leak_test_points()
  d <- spatialkit:::.sac_range_for_diagnostic(pts, "z", NULL, 1.0, 123L)
  expect_false(inherits(d, "sac_range"))
  expect_null(attributes(d))
  expect_equal(d, as.numeric(estimate_sac_range(pts, "z")))
})

test_that("the diagnostic reaches a cv_*() call that builds default folds", {
  skip_if_not_installed("gstat")
  pts <- leak_test_points()
  # lm_spatial_fit() is the suite's minimal spatial_fit backend (helper-lmfit.R).
  res <- collect_warnings(suppressMessages(
    cv_spatial(pts, "z", "w", fit_fn = function(tr) lm_spatial_fit(tr, "z", "w"),
               k = 10, seed = 7)))
  expect_true(any(grepl("block dimension .* < autocorrelation range", res$warnings)))
  expect_identical(res$value$n_folds_succeeded, res$value$n_folds_attempted)
})
