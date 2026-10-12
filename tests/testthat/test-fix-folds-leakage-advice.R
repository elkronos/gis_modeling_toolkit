# tests/testthat/test-fix-folds-leakage-advice.R
# ---------------------------------------------------------------------------
# The leakage warning of make_folds(method = "block_kfold") ended "Pass
# block_size = <range> or auto_range = TRUE" whatever the layout.  On a
# square extent at k = 5 that gives five folds only for a range up to a third
# of the side: up to a half it returned four folds without a warning, and
# above it both settings stopped with the single-block error.  The advice is
# now worked out for the layout, and these tests do what it says.
# ---------------------------------------------------------------------------

la_layout <- function(kind, n = 400, w = 1e5, h = 1e5, seed = 1) {
  set.seed(seed)
  xy <- if (kind == "uniform") cbind(runif(n, 0, w), runif(n, 0, h)) else {
    c0 <- cbind(runif(6, 0.1 * w, 0.9 * w), runif(6, 0.1 * h, 0.9 * h))
    i <- sample(6, n, TRUE)
    c0[i, ] + cbind(rnorm(n, 0, 0.03 * w), rnorm(n, 0, 0.03 * h))
  }
  sf::st_as_sf(data.frame(x = 4e5 + xy[, 1], y = 4.5e6 + xy[, 2]),
               coords = c("x", "y"), crs = 32617)
}

# Do what one piece of advice says and return what comes back.
la_follow <- function(pts, k, ...) {
  reduced <- FALSE
  val <- tryCatch(
    withCallingHandlers(
      r2_quiet(make_folds(pts, k = k, method = "block_kfold", seed = 1, ...)),
      spatialkit_k_reduced = function(w) { reduced <<- TRUE; invokeRestart("muffleWarning") },
      warning = function(w) invokeRestart("muffleWarning"),
      message = function(m) invokeRestart("muffleMessage")),
    error = function(e) conditionMessage(e))
  list(value = val, reduced = reduced)
}

# Check every claim an advice text makes about a layout.
la_check <- function(adv, pts, k, label) {
  size <- as.numeric(sub(".*?block_size = ([0-9.]+).*", "\\1", adv, perl = TRUE))
  got  <- la_follow(pts, k, block_size = size)
  if (grepl("enough for k = [0-9]+ folds", adv)) {
    # "Pass block_size = X: ... enough for k = K folds."
    expect_false(is.character(got$value), info = label)
    expect_identical(got$value$k, as.integer(k), info = label)
    expect_false(got$reduced, info = label)
    return(invisible("k"))
  }
  # Otherwise it names the largest size that still gives k folds, and that
  # size is below the range the warning is about.
  expect_match(adv, "The largest block_size that gives [0-9]+ folds is [0-9.]+", info = label)
  alt <- as.numeric(sub(".*The largest block_size that gives [0-9]+ folds is ([0-9.]+),.*", "\\1", adv))
  expect_lt(alt, size)
  alt_got <- la_follow(pts, k, block_size = alt)
  expect_false(is.character(alt_got$value), info = label)
  expect_identical(alt_got$value$k, as.integer(k), info = label)
  expect_false(alt_got$reduced, info = label)
  expect_match(adv, "method = \"nndm\"", fixed = TRUE, info = label)
  if (grepl("stops with an error", adv)) {
    # "block_size = X stops with an error": it does, with the single-block error.
    expect_true(is.character(got$value), info = label)
    expect_match(got$value, "single block", info = label)
    return(invisible("error"))
  }
  # "... returns B folds, not K": B folds come back, with the classed warning.
  expect_match(adv, "returns [0-9]+ folds, not [0-9]+", info = label)
  b <- as.integer(sub(".*returns ([0-9]+) folds, not [0-9]+.*", "\\1", adv))
  expect_false(is.character(got$value), info = label)
  expect_identical(got$value$k, b, info = label)
  expect_lt(b, k)
  expect_true(got$reduced, info = label)
  invisible("fewer")
}

test_that("the advice does what it says for ranges of 0.1 to 0.9 of the side", {
  outcomes <- character(0)
  for (kind in c("uniform", "clustered")) for (shape in c("square", "wide")) {
    w <- if (shape == "square") 1e5 else 2.4e5
    pts <- la_layout(kind, w = w, seed = if (kind == "uniform") 3 else 4)
    bb  <- sf::st_bbox(pts)
    reg <- sf::st_sf(geometry = sf::st_as_sfc(bb))
    for (frac in seq(0.1, 0.9, by = 0.1)) for (k in c(3L, 5L)) {
      rng <- frac * 1e5 + 0.37           # not a round number, as an estimate is not
      adv <- r2_quiet(spatialkit:::.block_range_advice(pts, reg, bb, rng, k))
      # The size named is at least the range it stands for.
      size <- as.numeric(sub(".*?block_size = ([0-9.]+).*", "\\1", adv, perl = TRUE))
      expect_gte(size, rng)
      outcomes <- c(outcomes, la_check(adv, pts, k,
                                       sprintf("%s %s range %.1f k %d", kind, shape, frac, k)))
    }
  }
  # The sweep covers all three outcomes; a test that only ever saw one would
  # not be testing the other two.
  expect_setequal(unique(outcomes), c("k", "fewer", "error"))
})

test_that("a range below 10 CRS units is quoted as itself, not rounded up to 1", {
  # ceiling() turned a range of 0.3712 km into "Pass block_size = 1".
  set.seed(2)
  km <- sf::st_as_sf(data.frame(x = runif(200, 0, 10), y = runif(200, 0, 10)),
                     coords = c("x", "y"), crs = NA)
  bb <- sf::st_bbox(km); reg <- sf::st_sf(geometry = sf::st_as_sfc(bb))
  adv <- spatialkit:::.block_range_advice(km, reg, bb, 0.3712, 5L)
  expect_match(adv, "block_size = 0.372", fixed = TRUE)
  expect_identical(spatialkit:::.block_size_at_least(0.3712), 0.372)
  expect_identical(spatialkit:::.block_size_at_least(37128.4), 37129)
})

test_that("the default call's warning can be followed, end to end", {
  skip_if_not_installed("gstat")
  # A field whose range is known, on a 100 km square.  What decides the
  # advice is the range estimate_sac_range() returns for it, 43.5 km and
  # 59.8 km for true effective ranges of 0.4 and 0.5 of the side: blocks of
  # the first leave a 2 x 2 grid, blocks of the second a single block.  (The
  # iterative variogram fit this was first written against returned 37 km
  # and 53 km at 0.3 and 0.4; the estimate at 0.3 is now 28 km, under a third
  # of the side, for which five range-sized blocks do exist.)
  sim <- function(eff_range, seed) {
    set.seed(seed); n <- 700
    x <- runif(n, 0, 1e5); y <- runif(n, 0, 1e5)
    d <- as.matrix(stats::dist(cbind(x, y)))
    z <- as.numeric(t(chol(exp(-3 * d / eff_range) + diag(1e-8, n))) %*% rnorm(n)) +
      rnorm(n, 0, sqrt(0.2))
    sf::st_as_sf(data.frame(x = 4e5 + x, y = 4.5e6 + y, z = z),
                 coords = c("x", "y"), crs = 32617)
  }
  seen <- character(0)
  for (frac in c(0.4, 0.5)) {
    dat <- sim(frac * 1e5, seed = 12)
    w <- character(0)
    f <- withCallingHandlers(
      r2_quiet(make_folds(dat, k = 5, method = "block_kfold", response_var = "z", seed = 1)),
      warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
    adv <- grep("block dimension .* < autocorrelation range", w, value = TRUE)
    # The default 4 x 4 grid has 25 km blocks and both ranges are above it.
    expect_length(adv, 1L)
    if (length(adv) != 1L) next
    expect_identical(f$k, 5L)                 # the diagnostic sizes nothing
    # The range the advice is about is above a third of the side.
    expect_gt(as.numeric(f$params$sac_range), 1e5 / 3)
    seen <- c(seen, la_check(adv, dat, 5L, sprintf("range %.1f", frac)))
    # Where the text offers auto_range = TRUE as well, it gives the same
    # outcome as the block size beside it.
    if (grepl("or auto_range = TRUE", adv, fixed = TRUE)) {
      size <- as.numeric(sub(".*?block_size = ([0-9.]+).*", "\\1", adv, perl = TRUE))
      a <- la_follow(dat, 5L, response_var = "z", auto_range = TRUE)
      b <- la_follow(dat, 5L, block_size = size)
      expect_identical(is.character(a$value), is.character(b$value))
      if (!is.character(a$value)) expect_identical(a$value$k, b$value$k)
      expect_identical(a$reduced, b$reduced)
    }
  }
  # Neither of these ranges fits five range-sized blocks, so neither message
  # may promise them: one returns fewer folds and the other a single block.
  expect_false("k" %in% seen)
  expect_setequal(seen, c("fewer", "error"))
})

test_that("auto_range says what it built, and warns when that is fewer folds", {
  skip_if_not_installed("gstat")
  set.seed(12); n <- 400
  x <- runif(n, 0, 1e5); y <- runif(n, 0, 1e5)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-3 * d / 3e4) + diag(1e-8, n))) %*% rnorm(n)) + rnorm(n, 0, sqrt(0.2))
  dat <- sf::st_as_sf(data.frame(x = 4e5 + x, y = 4.5e6 + y, z = z),
                      coords = c("x", "y"), crs = 32617)
  msgs <- character(0); reduced <- 0L
  f <- tryCatch(withCallingHandlers(
    r2_quiet(make_folds(dat, k = 5, method = "block_kfold", response_var = "z",
                        auto_range = TRUE, seed = 1)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
    spatialkit_k_reduced = function(w) { reduced <<- reduced + 1L; invokeRestart("muffleWarning") }),
    error = function(e) NULL)
  skip_if(is.null(f) || !isTRUE(f$params$block_size > 0), "no range-sized grid was built")
  expect_true(any(grepl("using as minimum block size\\. Over the .* extent that is a [0-9]+ x [0-9]+ grid of blocks",
                        msgs)))
  # Fewer folds than asked never comes back without the warning.
  expect_identical(f$k < 5L, reduced == 1L)
  expect_identical(f$params$k_requested, 5L)
})
