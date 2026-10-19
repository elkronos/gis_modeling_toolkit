# tests/testthat/test-fix-folds-k-reduced.R
# ---------------------------------------------------------------------------
# make_folds() returning fewer folds than asked.  It was a log line only, at
# four sites in three methods, and the result kept no record of the k asked
# for.  Each test here fails on the code before the fix.
# ---------------------------------------------------------------------------

kr_pts <- function(n = 300, side = 1e4, seed = 42) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, side), y = 5e6 + runif(n, 0, side),
                          z = rnorm(n), a = rnorm(n),
                          g = sample(letters[1:3], n, TRUE)),
               coords = c("x", "y"), crs = 32632)
}

# Every warning a call raises, with its classes and the condition itself.
kr_collect <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = val, conditions = w,
       messages = vapply(w, conditionMessage, ""),
       reduced = Filter(function(c) inherits(c, "spatialkit_k_reduced"), w))
}

test_that("each of the four k-reduction sites raises the classed warning", {
  u <- kr_pts()
  calls <- list(
    block_size = function() make_folds(u, k = 5, method = "block_kfold",
                                       block_size = 4000, seed = 1),
    block_grid = function() make_folds(u, k = 5, method = "block_kfold",
                                       block_nx = 3, block_ny = 1, seed = 1),
    groups     = function() make_folds(u, k = 5, method = "leave_location_out",
                                       group_var = "g", seed = 1),
    rows       = function() make_folds(u[1:7, ], k = 10, method = "random_kfold",
                                       seed = 1))
  got <- c(block_size = 4L, block_grid = 3L, groups = 3L, rows = 7L)
  asked <- c(block_size = 5L, block_grid = 5L, groups = 5L, rows = 10L)
  for (nm in names(calls)) {
    res <- r2_quiet(kr_collect(calls[[nm]]()))
    expect_length(res$reduced, 1L)
    cnd <- res$reduced[[1L]]
    expect_s3_class(cnd, "warning")
    # The counts travel on the condition, so a handler need not read text.
    expect_identical(cnd$k, got[[nm]], info = nm)
    expect_identical(cnd$k_requested, asked[[nm]], info = nm)
    expect_match(conditionMessage(cnd),
                 sprintf("returning %d folds, not the k = %d asked for", got[[nm]], asked[[nm]]),
                 info = nm)
    # ... and on the result, which used to keep no record of the k asked for.
    expect_identical(res$value$k, got[[nm]], info = nm)
    expect_identical(res$value$params$k_requested, asked[[nm]], info = nm)
    expect_length(res$value$folds, got[[nm]])
  }
  # tryCatch() by class sees it; it used to be invisible to any handler.
  expect_identical(
    r2_quiet(tryCatch(calls$block_size(), spatialkit_k_reduced = function(w) "caught")),
    "caught")
})

test_that("a k that is honoured raises nothing and is recorded as asked", {
  u <- kr_pts()
  for (f in list(
    r2_quiet(expect_no_warning(make_folds(u, k = 5, method = "random_kfold", seed = 1))),
    r2_quiet(expect_no_warning(make_folds(u, k = 5, method = "block_kfold", seed = 1))),
    r2_quiet(expect_no_warning(make_folds(u, k = 3, method = "leave_location_out",
                                          group_var = "g", seed = 1))))) {
    expect_identical(f$k, f$params$k_requested)
  }
  # The leave-one-out methods take no k; NA says none was given.
  loo <- r2_quiet(suppressWarnings(make_folds(u[1:40, ], method = "buffered_loo", buffer = 300)))
  expect_true(is.na(loo$params$k_requested))
})

test_that("empty cells lower k once, with one warning that counts the blocks holding points", {
  # Four clusters under a 3 x 3 grid: nine cells, four of them occupied.  The
  # cell count alone promised five folds; four came back and nothing warned.
  set.seed(5)
  cx <- c(500, 9500, 500, 5000); cy <- c(500, 9500, 9500, 5000); g <- rep(1:4, each = 100)
  cl <- sf::st_as_sf(data.frame(x = 5e5 + cx[g] + rnorm(400, 0, 150),
                                y = 5e6 + cy[g] + rnorm(400, 0, 150)),
                     coords = c("x", "y"), crs = 32632)
  bs <- as.numeric(diff(sf::st_bbox(cl)[c("xmin", "xmax")])) / 3.2
  res <- r2_quiet(kr_collect(make_folds(cl, k = 5, method = "block_kfold",
                                        block_size = bs, seed = 1)))
  expect_identical(c(res$value$params$grid_nx, res$value$params$grid_ny), c(3, 3))
  expect_identical(res$value$k, 4L)
  expect_length(res$reduced, 1L)
  expect_match(res$messages, "3 x 3 grid .* 4 of its blocks hold points")
})

test_that("the size the k warning names gives the k it promises", {
  u <- kr_pts()
  res <- r2_quiet(kr_collect(make_folds(u, k = 5, method = "block_kfold",
                                        block_size = 4000, seed = 1)))
  msg <- res$messages[[1L]]
  expect_match(msg, "For 5 folds pass block_size = [0-9.]+ \\(the largest that gives 5\\)")
  s <- as.numeric(sub(".*For 5 folds pass block_size = ([0-9.]+) .*", "\\1", msg))
  expect_lt(s, 4000)
  f <- r2_quiet(expect_no_warning(make_folds(u, k = 5, method = "block_kfold",
                                             block_size = s, seed = 1)))
  expect_identical(f$k, 5L)
  # "or pass k = 4 to keep these blocks": that works too, quietly.
  f4 <- r2_quiet(expect_no_warning(make_folds(u, k = 4, method = "block_kfold",
                                              block_size = 4000, seed = 1)))
  expect_identical(f4$k, 4L)
})

test_that("the single-block error names a block size that works", {
  u <- kr_pts()
  err <- tryCatch(r2_quiet(make_folds(u, k = 5, method = "block_kfold",
                                      block_size = 5100, seed = 1)),
                  error = function(e) conditionMessage(e))
  expect_match(err, "the block size \\(5100\\) produces a single block covering the whole extent")
  # It used to say "set auto_range = FALSE", which leads back to the default
  # grid and its leakage warning.
  expect_no_match(err, "auto_range = FALSE")
  expect_match(err, "the largest that gives 5 folds is [0-9.]+")
  s <- as.numeric(sub(".*the largest that gives 5 folds is ([0-9.]+)\\.$", "\\1", err))
  f <- r2_quiet(expect_no_warning(make_folds(u, k = 5, method = "block_kfold",
                                             block_size = s, seed = 1)))
  expect_identical(f$k, 5L)
})

test_that("params records the block size built, not only the size asked", {
  u <- kr_pts()
  for (bs in c(2600, 3000, 4000)) {
    f <- r2_quiet(suppressWarnings(make_folds(u, k = 4, method = "block_kfold",
                                              block_size = bs, seed = 1)))
    boxes <- lapply(sf::st_geometry(f$params$blocks), sf::st_bbox)
    wd <- vapply(boxes, function(b) b[["xmax"]] - b[["xmin"]], numeric(1))
    ht <- vapply(boxes, function(b) b[["ymax"]] - b[["ymin"]], numeric(1))
    expect_equal(rep(f$params$block_width, length(wd)), wd)
    expect_equal(rep(f$params$block_height, length(ht)), ht)
    # A minimum edge: built at least as large, and under twice as large.
    expect_identical(f$params$block_size, bs)
    expect_gte(f$params$block_width, bs); expect_lt(f$params$block_width, 2 * bs)
    expect_gte(f$params$block_height, bs); expect_lt(f$params$block_height, 2 * bs)
  }
  # Supplied blocks have no single width.
  cells <- sf::st_as_sf(sf::st_make_grid(u, n = c(3, 3)))
  fs <- r2_quiet(make_folds(u, k = 4, method = "block_kfold", blocks = cells, seed = 1))
  expect_true(is.na(fs$params$block_width) && is.na(fs$params$block_height))
})

test_that("an axis shorter than block_size is spanned by one row, shorter than the size", {
  # A corridor just under 40 km long (39.8 km between its outermost points)
  # and 5 km wide, blocks of 10 km.  Three fit along
  # it; across it there is no room for one, so a single row spans the width
  # and each block is as tall as the extent: half the size asked.  "At least
  # block_size" holds on the long axis only.
  set.seed(4)
  cor <- sf::st_as_sf(data.frame(x = 5e5 + runif(300, 0, 40000), y = 5e6 + runif(300, 0, 5000),
                                 a = rnorm(300)),
                      coords = c("x", "y"), crs = 32632)
  bb <- sf::st_bbox(cor)
  w <- bb[["xmax"]] - bb[["xmin"]]; h <- bb[["ymax"]] - bb[["ymin"]]
  bs <- 10000
  expect_lt(h, bs); expect_gt(w, 3 * bs)
  f <- r2_quiet(expect_no_warning(make_folds(cor, k = 3, method = "block_kfold",
                                             block_size = bs, seed = 1)))
  p <- f$params
  expect_identical(c(p$grid_nx, p$grid_ny), c(3, 1))
  # The long axis: one to under two times the size.
  expect_equal(p$block_width, w / 3)
  expect_gte(p$block_width, bs); expect_lt(p$block_width, 2 * bs)
  # The short axis: the extent itself, which is below the size.
  expect_equal(p$block_height, h)
  expect_lt(p$block_height, bs)
  # params describes the polygons returned.
  boxes <- lapply(sf::st_geometry(p$blocks), sf::st_bbox)
  expect_equal(vapply(boxes, function(b) b[["xmax"]] - b[["xmin"]], numeric(1)),
               rep(p$block_width, 3))
  expect_equal(vapply(boxes, function(b) b[["ymax"]] - b[["ymin"]], numeric(1)),
               rep(p$block_height, 3))

  # A response with a range of about 9 km, wider than the corridor.
  D <- as.matrix(stats::dist(sf::st_coordinates(cor)))
  cor$z <- as.numeric(t(chol(exp(-3 * D / 9000) + diag(1e-8, 300))) %*% rnorm(300)) +
    rnorm(300, 0, 0.3)

  # The sweep's table reports the same cells.
  sw <- r2_quiet(suppressWarnings(suppressMessages(
    cv_block_size_sweep(cor, "z", "a", fit_fn = r2_fit, block_sizes = bs, k = 3,
                        include_random = FALSE, sac = 6000, quiet = TRUE))))
  expect_equal(sw$block_width, w / 3)
  expect_equal(sw$block_height, h)
  expect_lt(sw$block_height, sw$block_size)

  # auto_range on the same corridor: a range above the corridor's width gets
  # one row of blocks as tall as the corridor, so shorter than the range,
  # and the message says what is built.
  msgs <- character(0)
  a <- r2_quiet(suppressWarnings(withCallingHandlers(
    make_folds(cor, k = 4, method = "block_kfold", auto_range = TRUE,
               response_var = "z", seed = 1),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })))
  rng <- as.numeric(a$params$sac_range)
  # The estimate is 6044 where this was written; the check needs only that
  # it is above the corridor's width (4969).
  skip_if(!(is.finite(rng) && rng > h),
          "the estimated range is not above the corridor's width on this platform")
  expect_identical(a$params$grid_ny, 1)
  expect_equal(a$params$block_height, h)
  expect_lt(a$params$block_height, rng)
  expect_gte(a$params$block_width, rng); expect_lt(a$params$block_width, 2 * rng)
  expect_match(msgs[grepl("estimated spatial autocorrelation range", msgs)],
               sprintf("a %s x 1 grid of blocks %s wide and %s tall",
                       format(a$params$grid_nx), format(signif(a$params$block_width, 4)),
                       format(signif(h, 4))))
})

test_that("a cv_*() call that builds fewer folds passes the classed warning on", {
  u <- kr_pts()
  res <- r2_quiet(kr_collect(suppressMessages(
    cv_spatial(u, "z", "a", fit_fn = r2_fit, k = 5, block_size = 4000, seed = 1))))
  expect_identical(res$value$n_folds_attempted, 4L)
  expect_length(res$reduced, 1L)
})

test_that("cv_block_size_sweep() gathers the k reductions into one warning", {
  # Four clusters: at the coarser sizes fewer than k blocks hold points.
  set.seed(5)
  cx <- c(500, 9500, 500, 5000); cy <- c(500, 9500, 9500, 5000); g <- rep(1:4, each = 30)
  cl <- sf::st_as_sf(data.frame(x = 5e5 + cx[g] + rnorm(120, 0, 400),
                                y = 5e6 + cy[g] + rnorm(120, 0, 400),
                                a = rnorm(120)),
                     coords = c("x", "y"), crs = 32632)
  cl$z <- cl$a + rnorm(120)
  res <- r2_quiet(kr_collect(suppressMessages(
    cv_block_size_sweep(cl, "z", "a", fit_fn = r2_fit, block_sizes = c(600, 3000),
                        k = 5, include_random = FALSE, sac = 1000, quiet = TRUE))))
  sw <- res$value
  expect_true(any(sw$k < 5L))
  # One warning from the sweep, none passed through from make_folds().
  expect_length(res$reduced, 0L)
  hit <- grepl("cv_block_size_sweep\\(\\): fewer than the k = 5 folds asked for in 1 of the 2 cross-validations \\(block size 3000: 4 folds\\)",
               res$messages)
  expect_identical(sum(hit), 1L)
  # The table says what was built at each size.
  expect_true(all(c("block_width", "block_height") %in% names(sw)))
  expect_true(all(sw$block_width >= sw$block_size & sw$block_width < 2 * sw$block_size))
})
