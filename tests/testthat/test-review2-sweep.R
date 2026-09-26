# tests/testthat/test-review2-sweep.R
# ---------------------------------------------------------------------------
# Second review round: cv_block_size_sweep().  Each test fails on the code before
# the fix.  Helpers are in helper-review2-folds.R.
# ---------------------------------------------------------------------------

# ---- cv_block_size_sweep() -------------------------------------------------

r2_sweep_layer <- function(x, y, seed = 1) {
  set.seed(seed)
  p <- r2_pts(x, y, a = rnorm(length(x)))
  p$z <- sin(sf::st_coordinates(p)[, 1] / 500) + 0.5 * p$a + rnorm(length(x), 0, 0.2)
  p
}

test_that("the sweep runs on points along one axis-parallel line", {
  set.seed(1)
  p <- r2_sweep_layer(5e5 + runif(60, 0, 6000), rep(5e6, 60))
  sw <- r2_quiet(cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 3, n_sizes = 4,
                                     sac = NA, quiet = TRUE))       # "no extent" before
  expect_gt(sum(sw$method == "block_kfold"), 0L)
  expect_true(all(sw$k == 3L))
})

test_that("the default ladder on a thin transect skips grids make_folds() would refuse", {
  set.seed(1)
  # A 6000 x 3 m transect.  The ladder's lower rung, 0.12 m, is a
  # 50000 x 25 grid of 1.25 million blocks, which used to abort the sweep with
  # an error about the units of a `block_size` nobody passed.
  p <- r2_sweep_layer(5e5 + c(0, 6000, runif(58, 0, 6000)), 5e6 + c(0, 3, runif(58, 0, 3)))
  sw <- r2_quiet(cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 5, n_sizes = 2,
                                     sac = NA, quiet = TRUE))
  blk <- sw[sw$method == "block_kfold", ]
  expect_equal(nrow(blk), 1L)
  expect_equal(blk$block_size, 1.5, tolerance = 1e-6)
  # Sizes the caller chose are refused up front, naming the argument.
  expect_error(cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, sac = NA, quiet = TRUE,
                                   block_sizes = c(0.001, 100)),
               "`block_sizes` 0.001 would each need a grid of more than")
})

test_that("the sweep warns when the default ladder cannot reach the range", {
  # A 5 km x 200 m corridor: the ladder stops at 100 m, while 1000 m blocks
  # would still give 5 along its length.
  set.seed(1)
  p <- r2_sweep_layer(5e5 + runif(60, 0, 5000), 5e6 + runif(60, 0, 200))
  expect_warning(
    r2_quiet(cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 4, n_sizes = 2,
                                 sac = 1000, quiet = TRUE)),
    "every block size in the default ladder .* is below the estimated autocorrelation range")
  # Not on a square, where no longer block would still give k blocks.
  sq <- r2_sweep_layer(runif(60, 0, 1000), runif(60, 0, 1000))
  expect_no_warning(
    r2_quiet(cv_block_size_sweep(sq, "z", "a", fit_fn = r2_fit, k = 4, n_sizes = 3,
                                 sac = 900, quiet = TRUE)))
})

test_that("a sac estimated in another CRS is converted to the sweep's units", {
  set.seed(1)
  d <- r2_pts(7e5 + runif(60, 0, 2000), 3.95e6 + runif(60, 0, 2000), crs = 32617,
              a = rnorm(60))
  d$z <- d$a + rnorm(60)
  sac_m <- structure(600, class = "sac_range", crs = sf::st_crs(32617))
  expect_warning(
    sw <- r2_quiet(cv_block_size_sweep(sf::st_transform(d, 2264), "z", "a",
                                       fit_fn = r2_fit, k = 3, n_sizes = 2,
                                       sac = sac_m, quiet = TRUE)),
    "`sac` was estimated in EPSG:32617, not in EPSG:2264")
  # 600 m is 1968.5 US survey feet.
  expect_equal(attr(sw, "sac_range"), 600 / 0.3048006, tolerance = 1e-3)
  # The same CRS, or no CRS recorded: used as given.
  sw2 <- r2_quiet(cv_block_size_sweep(d, "z", "a", fit_fn = r2_fit, k = 3, n_sizes = 2,
                                      sac = sac_m, quiet = TRUE))
  expect_equal(attr(sw2, "sac_range"), 600)
  expect_error(cv_block_size_sweep(d, "z", "a", fit_fn = r2_fit, k = 3, n_sizes = 2,
                                   sac = units::set_units(600, m), quiet = TRUE),
               "`sac` must be a plain number")
})

test_that("the plot reads every coverage column by closeness to its nominal level", {
  skip_if_not_installed("ggplot2")
  # cv_bayes() names coverage columns at full precision (coverage_97.5), and
  # a `metrics` function can add any; only coverage_50/80/95 were
  # recognised, and as higher-is-better, so every other level was captioned
  # "lower is better".  Over-coverage is miscalibration too.
  set.seed(1)
  p <- r2_sweep_layer(runif(60, 0, 1000), runif(60, 0, 1000))
  cov <- function(y, yhat) c(coverage_97.5 = mean(abs(y - yhat) < 2 * stats::sd(y)),
                             coverage_95 = mean(abs(y - yhat) < 1.96 * stats::sd(y)))
  caption <- function(metric) {
    sw <- r2_quiet(cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 3, n_sizes = 2,
                                       sac = NA, quiet = TRUE, metric = metric,
                                       metrics = cov))
    plot(sw)$labels$caption
  }
  expect_match(caption("coverage_97.5"), "closer to 0.975 is better")
  expect_match(caption("coverage_95"), "closer to 0.95 is better")
  expect_match(caption("RMSE"), "lower is better")
  expect_match(caption("R2"), "higher is better")
})
