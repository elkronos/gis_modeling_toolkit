# tests/testthat/test-review3-S5-folds.R
# ---------------------------------------------------------------------------
# Third review round: make_folds() and cv_block_size_sweep().  Each test
# fails on the code before the fix, except the one that pins the ladder
# where it did not change.  Helpers are in helper-review2-folds.R.
# ---------------------------------------------------------------------------

r3_collect_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

r3_sweep_layer <- function(x, y, seed = 1) {
  set.seed(seed)
  p <- r2_pts(x, y, a = rnorm(length(x)))
  p$z <- sin(sf::st_coordinates(p)[, 1] / 500) + 0.5 * p$a + rnorm(length(x), 0, 0.2)
  p
}

r3_cells <- function(bb, b) {
  d <- spatialkit:::.block_dims_from_size(bb, b)
  as.numeric(d$nx) * as.numeric(d$ny)
}


# ---- cv_block_size_sweep(): the default ladder's top -----------------------

test_that("the default ladder at k = 5 keeps all its sizes on a square and reaches side / 3", {
  # The ladder used to top out at half the shorter side, a 2 x 2 grid of four
  # cells, which k = 5 always dropped: five cross-validations for n_sizes = 6,
  # the last at ~0.30 of the side.  A range between that and side / 3 (a
  # 3 x 3 grid) then drew a warning that called the 0.30 rung "half the
  # shorter side" and said only blocks up to the longer side / k -- smaller
  # than rungs already run -- still gave k blocks.
  set.seed(4)
  p <- r3_sweep_layer(5e5 + runif(60, 0, 1000), 5e6 + runif(60, 0, 1000))
  bb <- sf::st_bbox(p)
  side <- min(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"])
  res <- r3_collect_warnings(r2_quiet(
    cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 5, n_sizes = 6,
                        sac = 0.31 * as.numeric(side), quiet = TRUE)))
  expect_length(res$warnings, 0L)
  sw <- res$value
  blk <- sw[sw$method == "block_kfold", ]
  expect_equal(nrow(blk), 6L)
  expect_equal(attr(sw, "n_fits"), 35L)
  top <- max(blk$block_size)
  expect_gt(top, as.numeric(side) / 3.05)
  expect_lte(top, as.numeric(side) / 2)
  expect_gte(r3_cells(bb, top), 5)
  expect_true(all(blk$k == 5L))
})

test_that("the ladder's top is unchanged where half the side already gives k blocks", {
  set.seed(4)
  p <- r3_sweep_layer(5e5 + runif(60, 0, 1000), 5e6 + runif(60, 0, 1000))
  bb <- sf::st_bbox(p)
  side <- as.numeric(min(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]))
  sw <- r2_quiet(cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 4, n_sizes = 3,
                                     include_random = FALSE, sac = NA, quiet = TRUE))
  expect_equal(max(sw$block_size), side / 2 * (1 - 1e-9))
  expect_equal(min(sw$block_size), side / 25)
})

test_that("the ladder warning on a corridor names the rung run and a size that gives k blocks", {
  set.seed(1)
  p <- r3_sweep_layer(5e5 + runif(60, 0, 5000), 5e6 + runif(60, 0, 200))
  bb <- sf::st_bbox(p)
  res <- r3_collect_warnings(r2_quiet(
    cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 4, n_sizes = 2,
                        sac = 1000, quiet = TRUE)))
  expect_length(res$warnings, 1L)
  msg <- res$warnings
  top_run <- max(res$value$block_size, na.rm = TRUE)
  expect_match(msg, sprintf("default ladder \\(up to %s, over the ",
                            format(signif(top_run, 3))))
  expect_no_match(msg, "half the shorter side")
  b <- as.numeric(sub(".*Blocks up to ([0-9.e+]+) still give k = 4 blocks.*", "\\1", msg))
  expect_true(is.finite(b))
  expect_gt(b, 1000)                       # past the range, so the advice works
  expect_gte(r3_cells(bb, b), 4)           # and the size it names gives k blocks
  expect_lt(r3_cells(bb, b * 1.05), 4)     # and is close to the largest that does
})


# ---- cv_block_size_sweep(): sizes the caller chose -------------------------

test_that("user block_sizes that give fewer than k blocks are warned about, naming a size that works", {
  set.seed(4)
  p <- r3_sweep_layer(5e5 + runif(60, 0, 1000), 5e6 + runif(60, 0, 1000))
  bb <- sf::st_bbox(p)
  res <- r3_collect_warnings(r2_quiet(
    cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 5,
                        block_sizes = c(100, 200, 400, 600), include_random = FALSE,
                        sac = NA, quiet = TRUE)))
  expect_equal(res$value$block_size, c(100, 200))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "`block_sizes` 400, 600 give fewer than k = 5 blocks")
  b <- as.numeric(sub(".*the largest size whose grid holds k blocks is ([0-9.]+).*", "\\1",
                      res$warnings))
  expect_gte(r3_cells(bb, b), 5)
  expect_lt(r3_cells(bb, b * 1.05), 5)
})

test_that("a units object for block_sizes is refused by name", {
  set.seed(4)
  p <- r3_sweep_layer(5e5 + runif(60, 0, 1000), 5e6 + runif(60, 0, 1000))
  expect_error(cv_block_size_sweep(p, "z", "a", fit_fn = r2_fit, k = 3,
                                   block_sizes = units::set_units(c(100, 200), m),
                                   sac = NA, quiet = TRUE),
               "`block_sizes` must be positive numbers, given as plain numbers in EPSG:32632 units")
})


# ---- make_folds(): the leakage diagnostic on a single row of blocks --------

test_that("a single row of blocks is compared with the range along the row only", {
  skip_if_not_installed("gstat")
  # Points on a 10 km line.  The one row spans the whole (zero) height and
  # borders no other block across it; min(w / nx, h / ny) compared the range
  # with 0 and warned on every call, and its advice (block_size = range)
  # shortened the blocks along the line.
  set.seed(11)
  n <- 300
  x <- sort(runif(n, 0, 10000))
  z <- as.numeric(t(chol(exp(-as.matrix(stats::dist(x)) / 50) + diag(1e-6, n))) %*% rnorm(n))
  p <- r2_pts(5e5 + x, rep(5e6, n), z = z)
  r <- suppressWarnings(as.numeric(estimate_sac_range(p, "z")))
  skip_if_not(is.finite(r) && r > 260 && r < 650, "range outside the fixture's window")
  # Automatic grid: 15 x 1, blocks 667 m long, longer than the range.
  res <- r3_collect_warnings(r2_quiet(
    make_folds(p, k = 5, method = "block_kfold", response_var = "z", seed = 1)))
  expect_equal(c(res$value$params$grid_nx, res$value$params$grid_ny), c(15, 1))
  expect_length(res$warnings, 0L)
  # block_nx: 1000 m blocks are silent, 250 m blocks warn.
  res <- r3_collect_warnings(r2_quiet(
    make_folds(p, k = 5, method = "block_kfold", response_var = "z", block_nx = 10, seed = 1)))
  expect_length(res$warnings, 0L)
  res <- r3_collect_warnings(r2_quiet(
    make_folds(p, k = 5, method = "block_kfold", response_var = "z", block_nx = 40, seed = 1)))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "block_nx/block_ny yield blocks smaller than autocorrelation range")

  # A 10 km x 100 m corridor: the 15 x 1 grid is compared along its length.
  set.seed(11)
  y <- runif(n, 0, 100)
  z2 <- as.numeric(t(chol(exp(-as.matrix(stats::dist(cbind(x, y))) / 50) +
                            diag(1e-6, n))) %*% rnorm(n))
  cp <- r2_pts(5e5 + x, 5e6 + y, z = z2)
  res <- r3_collect_warnings(r2_quiet(
    make_folds(cp, k = 5, method = "block_kfold", response_var = "z", seed = 1)))
  expect_equal(c(res$value$params$grid_nx, res$value$params$grid_ny), c(15, 1))
  expect_false(any(grepl("block dimension", res$warnings)))
})


# ---- make_folds(): a single block ------------------------------------------

test_that("a block size that leaves one block stops without announcing a lowered k", {
  set.seed(2)
  p <- r2_pts(5e5 + runif(60, 0, 1000), 5e6 + runif(60, 0, 1000))
  lines <- capture_spatialkit_log(
    expect_error(make_folds(p, k = 5, method = "block_kfold", block_size = 600, seed = 1),
                 "the block size \\(600\\) produces a single block"),
    level = logger::INFO)
  expect_false(log_has(lines, "Reducing k to match"))
  # Two blocks or more still lower k, and say so.
  lines <- capture_spatialkit_log(
    f <- make_folds(p, k = 5, method = "block_kfold", block_size = 400, seed = 1),
    level = logger::INFO)
  expect_true(log_has(lines, "block_size produces only 4 blocks \\(< k = 5\\). Reducing k"))
  expect_equal(f$k, 4L)
})


# ---- make_folds(): the 1,000,000-block guard -------------------------------

test_that("the grid guard names the argument that produced the grid, and refuses an overflow", {
  set.seed(1)
  p <- r2_pts(5e5 + runif(100, 0, 1000), 5e6 + runif(100, 0, 1000))
  e <- expect_error(make_folds(p, k = 5, method = "block_kfold", block_nx = 2000,
                               block_ny = 1000),
                    "2000 x 1000 = 2,000,000 cells, above the 1,000,000 this function will build")
  expect_match(conditionMessage(e), "`block_nx`/`block_ny`")
  expect_no_match(conditionMessage(e), "block_size|unset")
  e <- expect_error(make_folds(p, k = 5, method = "block_kfold", block_multiplier = 1e6),
                    "above the 1,000,000 this function will build")
  expect_match(conditionMessage(e), "`block_multiplier` \\(1e\\+06\\) x k \\(5\\)")
  expect_no_match(conditionMessage(e), "unset")
  # nx * ny overflowed to Inf and slipped past the guard into st_make_grid().
  expect_error(make_folds(p, k = 5, method = "block_kfold", block_size = 1e-200),
               "= more than 1e308 cells, above the 1,000,000 this function will build. Check that `block_size` \\(1e-200\\)")
})


# ---- argument validation ---------------------------------------------------

test_that("phi, min_train and block_multiplier are refused by name", {
  set.seed(1)
  p <- r2_pts(5e5 + runif(100, 0, 2000), 5e6 + runif(100, 0, 1000))
  pp <- sf::st_as_sf(sf::st_make_grid(p, n = c(10, 5), what = "centers"))
  expect_error(make_folds(p, method = "nndm", prediction_points = pp,
                          phi = units::set_units(100, m)),
               "`phi` must be a single non-negative number")
  expect_error(make_folds(p, method = "nndm", prediction_points = pp,
                          min_train = units::set_units(0.5, 1)),
               "`min_train` must be a single number in \\(0, 1\\)")
  for (bm in list(NA, c(1, 3), units::set_units(3, 1), "3", -1, 0, Inf, numeric(0)))
    expect_error(make_folds(p, k = 5, method = "block_kfold", block_multiplier = bm),
                 "`block_multiplier` must be a single positive number",
                 info = paste(format(bm), collapse = ","))
  # A valid fractional multiplier is still honoured.
  expect_equal(make_folds(p, k = 4, method = "block_kfold",
                          block_multiplier = 2.5)$params$block_multiplier, 2.5)
})

test_that("a misspelt response_var is an error for block_kfold whether or not auto_range is set", {
  set.seed(1)
  p <- r2_pts(5e5 + runif(60, 0, 1000), 5e6 + runif(60, 0, 1000), z = rnorm(60))
  expect_error(make_folds(p, k = 5, method = "block_kfold", response_var = "nope", seed = 1),
               "`response_var` 'nope' is not a column of `points_sf`")
  expect_error(make_folds(p, k = 5, method = "block_kfold", response_var = "nope",
                          auto_range = TRUE, seed = 1),
               "`response_var` 'nope' is not a column of `points_sf`")
  expect_error(make_folds(p, k = 5, method = "block_kfold", response_var = c("z", "z"), seed = 1),
               "`response_var` must be a single column name")
  # The methods that never read it still ignore it.
  expect_equal(make_folds(p, k = 5, method = "random_kfold", response_var = "nope", seed = 1)$k, 5L)
})


# ---- make_folds(auto_range = TRUE): the fallback warning's reason ----------

test_that("the auto_range fallback warning says why a bare NA came back", {
  skip_if_not_installed("gstat")
  set.seed(1)
  p20 <- r2_pts(5e5 + runif(20, 0, 1000), 5e6 + runif(20, 0, 1000), z = rnorm(20))
  expect_warning(r2_quiet(make_folds(p20, k = 3, method = "block_kfold", auto_range = TRUE,
                                     response_var = "z", seed = 1)),
                 "no autocorrelation range was identified \\(20 points, fewer than the 30")
  pc <- r2_pts(5e5 + runif(60, 0, 1000), 5e6 + runif(60, 0, 1000), z = rep(1, 60))
  expect_warning(r2_quiet(make_folds(pc, k = 3, method = "block_kfold", auto_range = TRUE,
                                     response_var = "z", seed = 1)),
                 "no autocorrelation range was identified \\(the response is constant\\)")
})
