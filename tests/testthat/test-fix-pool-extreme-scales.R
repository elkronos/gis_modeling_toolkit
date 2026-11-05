# tests/testthat/test-fix-pool-extreme-scales.R
# ---------------------------------------------------------------------------
# The test for no spread reads the raw values and the column centred on its
# median and scaled by its largest distance from it, so a spread far from 0
# is found, while values that differ only by their rounding (a standard
# deviation of 3 units in the last place or less) still have none; where
# the raw values lose their spread to rounding the fit runs on that column
# and is put back in the column's units; values whose squares the
# arithmetic cannot hold get NA pooled columns with a reason; and the
# profile says what each of its two tests of a constant found (BUG-21).
# ---------------------------------------------------------------------------

.xs_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# 800 rows in the 64 cells of an 8 x 8 km grid, with a smooth response.
.xs_layer <- function() {
  .xs_seed(2121)
  x <- stats::runif(800, 0, 8000); y <- stats::runif(800, 0, 8000)
  z <- sin(x / 2000) + cos(y / 1500) + stats::rnorm(800, sd = 0.5)
  L <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, z = z), coords = c("x", "y"),
                    crs = 32632)
  L$cell_id <- as.integer(8L * floor(y / 1000) + floor(x / 1000) + 1L)
  L
}

# The same rows as the engine reads them.
.xs_engine <- function(L) {
  u <- sort(unique(L$cell_id))
  gi <- match(L$cell_id, u)
  list(gi = gi, K = length(u), n = tabulate(gi, length(u)), xy = sf::st_coordinates(L))
}

.xs_catch <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w, classes = vapply(w, function(c) class(c)[1L], ""))
}

.xs_bounds <- paste0("which is values up to 1e150 in absolute value whose largest distance ",
                     "from their median is at least 1e-150")

# One unit in the last place of a positive value: the power of two at or
# below it times the machine epsilon.  log2() of a value just below a power
# of two can round up to it.
.xs_ulp <- function(v) {
  p <- 2^min(floor(log2(v)), 1023)
  if (p > v) p <- p / 2
  p * .Machine$double.eps
}

test_that("1e15 plus a spread of 1 is pooled as a column with spread", {
  L <- .xs_layer()
  L$offset <- 1e15 + L$z
  L$shift  <- L$offset - 1e15        # exact: the values lie within a factor 2 of 1e15
  a <- .xs_catch(summarize_by_cell(L, "offset", id_col = "cell_id", estimate = "pooled"))
  b <- .xs_catch(summarize_by_cell(L, "shift", id_col = "cell_id", estimate = "pooled"))
  expect_false("spatialkit_pool_degenerate" %in% a$classes)
  expect_identical(a$classes, b$classes)
  ra <- attr(a$value, "pooled")$columns$resp_offset
  rb <- attr(b$value, "pooled")$columns$resp_shift
  expect_identical(ra$status, rb$status)
  expect_identical(ra$status, "ok")
  expect_false(ra$constant)
  # The weights and every number in squared units are those of the spread
  # alone; the pooled values are, to the half unit in the last place that
  # storing them at 1e15 leaves (2^-4).
  expect_equal(a$value$..pooled_w_resp_offset, b$value$..pooled_w_resp_shift, tolerance = 1e-12)
  for (f in c("sigma2", "tau2", "own_weight", "loo", "loo_se", "loo_one_cell", "structure_p"))
    expect_equal(ra[[f]], rb[[f]], tolerance = 1e-12, info = f)
  expect_identical(ra$structured, rb$structured)
  expect_lte(max(abs((a$value$..pooled_resp_offset - 1e15) - b$value$..pooled_resp_shift)),
             2^-4 + 1e-12)
  # The plain columns are the plain means of the values themselves.
  plain <- summarize_by_cell(L, "offset", id_col = "cell_id", estimate = "mean")
  expect_identical(a$value$resp_mean_offset, plain$resp_mean_offset)
})

test_that("1e6 plus 1e-9 x noise has spread, and the engine's fit is that of the noise", {
  L <- .xs_layer(); g <- .xs_engine(L)
  y <- 1e6 + 1e-9 * L$z
  # By the raw values alone the column has no spread.
  expect_true(.no_spread(sum((y - mean(y))^2), y))
  A <- .pool_fit(y, g$gi, g$K, g$xy, use_fnn = FALSE)
  E <- .pool_fit((y - 1e6) * 1e9, g$gi, g$K, g$xy, use_fnn = FALSE)
  expect_false(A$constant)
  expect_false(A$flat)
  expect_identical(A$status, E$status)
  expect_equal(A$w, E$w, tolerance = 1e-12)
  for (f in c("sigma2", "tau2", "Q", "loo_pooled", "loo_pooled_se", "loo_mean", "loo_mean_se",
              "loo_one", "loo_one_se"))
    expect_equal(A[[f]] * 1e18, E[[f]], tolerance = 1e-12, info = f)
  expect_equal(A$f_p, E$f_p, tolerance = 1e-10)
  expect_identical(A$n, E$n)
  expect_identical(A$cen, E$cen)
  expect_identical(A$nb, E$nb)
  # The sums and plain means are those of the column.
  expect_identical(A$S, as.numeric(rowsum(y, g$gi)))
  expect_identical(A$m, A$S / A$n)
  # One unit in the last place of 1e6 is 1.16e-10, 0.116 in units of the noise.
  expect_lt(max(abs((A$value - 1e6) * 1e9 - E$value)), 0.116)
  expect_lt(max(abs((A$target - 1e6) * 1e9 - E$target)), 0.116)
  expect_lt(max(abs(A$resid_pooled * 1e9 - E$resid_pooled)), 0.116)
  expect_lt(max(abs(A$resid_mean * 1e9 - E$resid_mean)), 0.116)
  # The same through the column's summary.
  P <- .pool_column(y, L$cell_id, g$xy, obs_level = 0.9, want = "pooled", use_fnn = FALSE)
  expect_false(P$record$constant)
  expect_identical(P$record$status, E$status)
  expect_true(all(is.finite(P$record$offsets)))
  # A column of one value is constant at any size, and is not refused.
  for (v in c(0.1, 1e15, 1e200, 1e-170)) {
    C <- .pool_fit(rep(v, length(g$gi)), g$gi, g$K, g$xy, use_fnn = FALSE)
    expect_true(C$constant, info = format(v))
    expect_identical(C$status, "no spread inside any cell", info = format(v))
    expect_identical(C$value, C$m, info = format(v))
  }
})

test_that("a column far from 0 that varies only between cells is not called constant", {
  L <- .xs_layer(); g <- .xs_engine(L)
  # Four levels 1 apart, each held by whole cells: exact at 1e15, where a
  # unit in the last place is 0.125, so the levels are 8 such units apart.
  y <- 1e15 + (g$gi %% 4)
  expect_true(.no_spread(sum((y - mean(y))^2), y))
  A <- .pool_fit(y, g$gi, g$K, g$xy, use_fnn = FALSE)
  E <- .pool_fit(y - 1e15, g$gi, g$K, g$xy, use_fnn = FALSE)
  expect_false(A$constant)
  expect_true(A$flat)
  expect_identical(A$status, "no spread inside any cell")
  expect_identical(A$status, E$status)
  expect_identical(A$value, A$m)
  expect_identical(A$w, rep(1, g$K))
  expect_identical(A$sigma2, 0)
  # Between the cells the column has spread, and the fit measures it.
  expect_gt(A$tau2, 0)
  for (f in c("tau2", "Q", "loo_pooled", "loo_mean", "loo_one", "loo_one_se"))
    expect_equal(A[[f]], E[[f]], tolerance = 1e-12, info = f)
  expect_lte(max(abs((A$target - 1e15) - E$target)), 2^-4)
  # Levels one unit in the last place apart are rounding: a standard
  # deviation of about one such unit, under the bound of three.
  y <- 1e15 + 0.125 * (g$gi %% 4)
  expect_lt(stats::sd(y - 1e15) / .xs_ulp(1e15), .POOL_ULPS)
  A <- .pool_fit(y, g$gi, g$K, g$xy, use_fnn = FALSE)
  expect_true(A$constant)
  expect_identical(A$status, "no spread inside any cell")
  expect_identical(A$w, rep(1, g$K))
  expect_identical(A$tau2, 0)
  expect_identical(A$loo_pooled, NA_real_)
})

test_that("the scaled reading never takes away a spread the raw values show", {
  L <- .xs_layer(); g <- .xs_engine(L)
  # Cells of +1 and -1 with a spread of 12 units in the last place inside
  # them.  Centred on their median, 1, the values are larger than they are,
  # and there that spread is rounding; on the raw values it is not.
  .xs_seed(5)
  sgn <- ifelse(seq_len(g$K) %% 9 < 5, 1, -1)
  y <- sgn[g$gi] * (1 + 12 * .Machine$double.eps * stats::rnorm(length(g$gi)))
  res <- y - (as.numeric(rowsum(y, g$gi)) / g$n)[g$gi]
  expect_false(.no_spread(sum(res^2), y) || .no_spread(.ss_within(res, g$gi, g$n), y))
  sc <- .spread_scaled(y)
  expect_equal(sc$centre, 1, tolerance = 1e-14)
  expect_true(.spread_readings(sc, g$gi, g$n)$flat)
  A <- .pool_fit(y, g$gi, g$K, g$xy, use_fnn = FALSE)
  expect_false(A$flat)
  expect_false(A$constant)
  expect_true(A$status %in% c("ok", "no signal between cells"))
  expect_gt(A$sigma2, 0)
})

test_that("the column of .spread_scaled() is exact, and its largest value lies in [1, 2)", {
  L <- .xs_layer()
  for (v in list(1e15 + L$z, 1e6 + 1e-9 * L$z, L$z * 1e200, L$z * 1e-170, L$z,
                 c(.Machine$double.xmax, -.Machine$double.xmax / 2))) {
    s <- .spread_scaled(v)
    expect_true(all(is.finite(s$y)))
    expect_gte(max(abs(s$y)), 1)
    expect_lt(max(abs(s$y)), 2)
    expect_identical(s$spread, max(abs(s$y)) * s$scale)
    expect_identical(log2(s$scale), round(log2(s$scale)))   # a power of two
    expect_equal(s$centre + s$scale * s$y, v, tolerance = 1e-14)
    expect_identical(s$centre, stats::median(v))
    # `unit` is one unit in the last place of the largest absolute value,
    # in the units of the scaled column.
    expect_identical(s$unit * s$scale, .xs_ulp(max(abs(v))))
  }
  # Within a factor of two of the median the scaled column holds the
  # differences exactly.
  v <- 1e15 + L$z
  s <- .spread_scaled(v)
  expect_identical(s$centre + s$scale * s$y, v)
  expect_identical(s$y * s$scale, v - stats::median(v))
  # A column of one value gives zeros, at any size.
  for (c0 in c(0.1, 1e200, 1e-170, 0)) {
    s <- .spread_scaled(rep(c0, 5))
    expect_identical(s$y, numeric(5), info = format(c0))
    expect_identical(s$spread, 0, info = format(c0))
    expect_identical(s$unit, .Machine$double.eps, info = format(c0))
  }
  # A largest value just below a power of two, where log2() rounds up to
  # it: the first divisor stays at or below it, so `unit` is a unit in the
  # last place of that value (512 times the machine epsilon), not of 1024.
  v <- c(1024 * (1 - .Machine$double.eps / 2), 1000)
  expect_identical(log2(v[1L]), 10)
  s <- .spread_scaled(v)
  expect_gte(max(abs(s$y)), 1)
  expect_lt(max(abs(s$y)), 2)
  expect_identical(s$unit * s$scale, 512 * .Machine$double.eps)
  expect_identical(s$centre + s$scale * s$y, v)
})

test_that("values whose squares the arithmetic cannot hold get NA pooled columns and a reason", {
  expect_identical(.POOL_Y_MAX, 1e150)
  expect_identical(.POOL_Y_MIN, 1e-150)
  L <- .xs_layer()
  L$big  <- L$z * 1e200
  L$tiny <- L$z * 1e-170
  L$one  <- 1e200
  r <- .xs_catch(summarize_by_cell(L, "z", c("big", "tiny", "one"), id_col = "cell_id",
                                   estimate = "pooled", obs_level = 0.9))
  w <- Filter(function(c) inherits(c, "spatialkit_pool_unavailable"), r$warnings)
  expect_length(w, 1L)
  expect_identical(w[[1L]]$reason, "values out of range")
  expect_identical(w[[1L]]$columns, c("big", "tiny"))
  expect_identical(conditionMessage(w[[1L]]), paste0(
    "summarize_by_cell(): the pooled columns of 'big' and 'tiny' are NA: their values lie ",
    "beyond what the arithmetic of pooled values holds, ", .xs_bounds, "; beyond these ",
    "bounds the squares of many values can overflow or lose their precision. Rescaled by a ",
    "power of ten to within that, the columns are pooled. The plain columns are unaffected."))
  rec <- attr(r$value, "pooled")$columns
  for (v in c("pred_big", "pred_tiny")) {
    expect_identical(rec[[v]]$status, "values out of range", info = v)
    expect_false(rec[[v]]$constant, info = v)
    for (col in c("..pooled_", "..pooled_w_", "..pooled_target_", "..obs_lo_", "..obs_hi_")) {
      expect_true(paste0(col, v) %in% names(r$value), info = paste(col, v))
      expect_true(all(is.na(r$value[[paste0(col, v)]])), info = paste(col, v))
    }
    for (f in c("sigma2", "tau2", "own_weight", "loo", "loo_one_cell"))
      expect_identical(rec[[v]][[f]], NA_real_, info = paste(v, f))
    expect_identical(rec[[v]]$offsets, c(NA_real_, NA_real_), info = v)
  }
  # The print says what the record says.
  out <- utils::capture.output(print(attr(r$value, "pooled")))
  expect_identical(sum(grepl("status      : values out of range", out, fixed = TRUE)), 2L)
  # The plain columns are as they were.
  plain <- .xs_catch(summarize_by_cell(L, "z", c("big", "tiny", "one"), id_col = "cell_id",
                                       estimate = "mean"))$value
  for (v in c("big", "tiny"))
    expect_identical(r$value[[paste0("pred_mean_", v)]], plain[[paste0("pred_mean_", v)]], info = v)
  # A column of one value at 1e200 has no spread, and says so.
  expect_identical(rec$pred_one$status, "no spread inside any cell")
  expect_true(rec$pred_one$constant)
  expect_true("spatialkit_pool_degenerate" %in% r$classes)
  # The response beside them is pooled as it is alone.
  alone <- .xs_catch(summarize_by_cell(L, "z", id_col = "cell_id", estimate = "pooled",
                                       obs_level = 0.9))$value
  expect_identical(r$value$..pooled_resp_z, alone$..pooled_resp_z)
  # One column alone is "its", and the bounds hold at the edges.
  one <- .xs_catch(summarize_by_cell(L, "big", id_col = "cell_id", estimate = "pooled"))
  w1 <- Filter(function(c) inherits(c, "spatialkit_pool_unavailable"), one$warnings)
  expect_match(conditionMessage(w1[[1L]]), "of 'big' are NA: its values lie beyond", fixed = TRUE)
  expect_match(conditionMessage(w1[[1L]]), "the column is pooled", fixed = TRUE)
  g <- .xs_engine(L)
  at <- function(y) .pool_fit(y, g$gi, g$K, g$xy, use_fnn = FALSE)$status
  zs <- L$z / max(abs(L$z))
  expect_false(identical(at(zs * 1e150), "values out of range"))
  expect_identical(at(zs * 1.0000001e150), "values out of range")
  d <- L$z - stats::median(L$z)
  zt <- d / max(abs(d))
  expect_false(identical(at(zt * 1.0000001e-150), "values out of range"))
  expect_identical(at(zt * 0.9999999e-150), "values out of range")
  # A column beyond the bounds that does not vary inside its cells and does
  # between them is refused too: what it has between the cells is squares.
  flat_big <- 1e200 * (1 + g$gi %% 4)
  expect_identical(at(flat_big), "values out of range")
  F <- .pool_fit(flat_big, g$gi, g$K, g$xy, use_fnn = FALSE)
  expect_false(F$constant)
  expect_true(F$flat)
  expect_identical(F$tau2, NA_real_)
  # At the bounds themselves the values are held: an absolute value of 1e150
  # and a largest distance from the median of 1e-150.
  edge <- function(y) .pool_out_of_range(y, .spread_scaled(y))
  expect_false(edge(c(-1e150, 0, 1e150)))
  expect_true(edge(c(-1e150, 0, 1e150 * (1 + 1e-15))))
  expect_identical(.spread_scaled(c(-1e-150, 0, 1e-150))$spread, 1e-150)
  expect_false(edge(c(-1e-150, 0, 1e-150)))
  expect_true(edge(c(-1e-150, 0, 1e-150) * (1 - 1e-15)))
})

test_that("the profile's loo reads the engine's test of a constant", {
  L <- .xs_layer()
  L$z <- 1e6 + 1e-9 * L$z
  r <- .xs_catch(resolution_profile(L, "z", estimate = "pooled", obs_per_cell = c(18, 4.5),
                                    sac = NULL))
  p <- r$value
  e <- attr(p, "estimate")
  expect_false(e$constant)
  expect_true(all(is.finite(p$loo)))
  expect_true(all(is.finite(p$own_weight)))
  expect_true(is.finite(e$loo_one_cell))
  # The error of one mean is that of the noise, in the column's units.
  noise <- (L$z - 1e6) * 1e9
  e1 <- (noise - mean(noise)) * length(noise) / (length(noise) - 1)
  expect_equal(e$loo_one_cell * 1e18, mean(e1^2), tolerance = 1e-12)
  # Cp's test reads the column the same way where the raw values show no
  # spread, and scores it too (BUG-36).
  expect_false(attr(p, "one_cell")$constant)
  expect_false("spatialkit_constant_response" %in% r$classes)
  expect_true(all(is.finite(p$cp)))
  expect_false("spatialkit_pool_unavailable" %in% r$classes)
  out <- utils::capture.output(print(p))
  expect_false(any(grepl("one cell: loo NA", out, fixed = TRUE)))
})

test_that("where the two tests of a constant disagree, the profile says what each found", {
  L <- .xs_layer()
  # The two tests read the values the same way (BUG-36): 1e15 plus a spread
  # of 1 is constant to neither.
  L$off <- 1e15 + L$z
  r <- .xs_catch(resolution_profile(L, "off", estimate = "pooled", obs_per_cell = c(18, 4.5),
                                    sac = NULL))
  expect_false(attr(r$value, "one_cell")$constant)
  expect_false(attr(r$value, "estimate")$constant)
  expect_false("spatialkit_constant_response" %in% r$classes)
  expect_true(all(is.finite(r$value$cp)))
  # They read different rows: the test for `cp` the subsample's, the
  # estimate's the layer's.  A response of 5 with one row of 6, and a
  # subsample of 100 that misses that row.
  .xs_seed(11)
  d <- sf::st_as_sf(data.frame(x = stats::runif(600, 0, 4000), y = stats::runif(600, 0, 4000),
                               v = 5), coords = c("x", "y"), crs = 32615)
  d$v[17] <- 6
  r <- .xs_catch(resolution_profile(d, "v", estimate = "pooled", sample_n = 100, n_levels = 3,
                                    sac = NULL, seed = 1))
  p <- r$value
  expect_true(attr(p, "one_cell")$constant)
  expect_identical(attr(p, "one_cell")$n, 100L)
  expect_false(attr(p, "estimate")$constant)
  expect_identical(attr(p, "estimate")$n, 600L)
  expect_true(all(is.finite(p$loo)))
  # The warning no longer says that no cell count fits the response better
  # than another, beside a `loo` that picks one, and says why the two differ.
  k <- Filter(function(c) inherits(c, "spatialkit_constant_response"), r$warnings)
  expect_length(k, 1L)
  expect_identical(k[[1L]]$n, 100L)
  expect_identical(conditionMessage(k[[1L]]), paste0(
    "resolution_profile(): the test for `cp` finds the response constant over the 100 scored ",
    "row(s) of the subsample (its total sum of squares, 0, is rounding error at most on their ",
    "scale), so `cp`, `moran_i` and `moran_z` are NA at every level. The estimate reads the 600 ",
    "row(s) of the layer with a value, and finds spread there: `loo` is scored."))
  expect_false("spatialkit_pool_unavailable" %in% r$classes)
  out <- utils::capture.output(print(p))
  expect_true(any(grepl(paste0("one cell    : rss [^,]+, constant on the rows the test for cp ",
                               "reads \\(one mean"), out)))
  expect_false(any(grepl("a constant", out, fixed = TRUE)))
  expect_false(any(grepl("raw values", c(out, conditionMessage(k[[1L]])), fixed = TRUE)))
  expect_error(select_resolution(p, criterion = "cp"), paste0(
    "select_resolution(): `cp` is NA at every level because the test for `cp` finds the ",
    "response constant on the 100 scored row(s) it reads: its total sum of squares is rounding ",
    "error at most on their scale. The estimate reads 600 row(s), and finds spread there, so ",
    "`loo` is scored: select_resolution(criterion = \"loo\") reads it."), fixed = TRUE)
  expect_error(select_resolution(p, criterion = "moran_z", diagnostic = TRUE),
               "`moran_z` is NA at every level because the test for `cp`", fixed = TRUE)
  expect_true(select_resolution(p, criterion = "loo")$best %in% p$levels)
  expect_s3_class(summary(p), "resolution_summary")
  # A response both tests find constant keeps the words it had.
  set.seed(7)
  L$mix <- sample(c(0.1, 0.3 / 3), nrow(L), TRUE)
  r <- .xs_catch(resolution_profile(L, "mix", estimate = "pooled", obs_per_cell = c(18, 4.5),
                                    sac = NULL))
  k <- Filter(function(c) inherits(c, "spatialkit_constant_response"), r$warnings)
  expect_match(conditionMessage(k[[1L]]), paste0(
    "rounding error at most), so no cell count fits it better than another: `cp`, `moran_i`, ",
    "`moran_z`, `loo`, `loo_se` and `own_weight` are NA at every level."), fixed = TRUE)
  out <- utils::capture.output(print(r$value))
  expect_true(any(grepl(", a constant (one mean", out, fixed = TRUE)))
  expect_error(select_resolution(r$value, criterion = "cp"),
               "because the response is constant: its total sum of squares", fixed = TRUE)
})

test_that("a profile of values beyond the arithmetic says why its loo is NA", {
  L <- .xs_layer()
  L$z <- L$z * 1e200
  for (est in c("pooled", "mean")) {
    r <- .xs_catch(resolution_profile(L, "z", estimate = est, obs_per_cell = c(18, 4.5),
                                      sac = NULL))
    q <- r$value
    expect_true(all(is.na(q$loo)), info = est)
    expect_true(all(is.na(q$loo_se)), info = est)
    expect_true(all(is.na(q$own_weight)), info = est)
    e <- attr(q, "estimate")
    expect_false(e$constant, info = est)
    expect_identical(e$loo_one_cell, NA_real_, info = est)
    expect_identical(e$loo_one_cell_se, NA_real_, info = est)
    # The test for `cp` finds these values constant (their sum of squares
    # is Inf, and so is its bound).  One warning says why every column that
    # reads the values is NA, and none calls them constant.
    expect_true(attr(q, "one_cell")$constant, info = est)
    w <- Filter(function(c) inherits(c, "spatialkit_pool_unavailable"), r$warnings)
    expect_length(w, 1L)
    expect_identical(w[[1L]]$reason, "values out of range", info = est)
    expect_identical(w[[1L]]$columns, "z", info = est)
    expect_identical(conditionMessage(w[[1L]]), paste0(
      "resolution_profile(): the values of 'z' lie beyond what the arithmetic of the ",
      "profile holds, ", .xs_bounds, "; beyond these bounds the squares of many values can ",
      "overflow or lose their precision. So `cp`, `moran_i`, `moran_z`, ",
      if (est == "pooled") "`loo`, `loo_se` and `own_weight`" else "`loo` and `loo_se`",
      " are NA at every level. Rescaled by a power of ten to within that, the response is ",
      "scored."), info = est)
    expect_false("spatialkit_constant_response" %in% r$classes, info = est)
    out <- utils::capture.output(print(q))
    expect_true(any(grepl("one cell: loo NA (the values are out of range)", out, fixed = TRUE)),
                info = est)
    expect_false(any(grepl("constant", out, fixed = TRUE)), info = est)
    expect_error(select_resolution(q, criterion = "loo"), paste0(
      "select_resolution(): `loo` is NA at every level because the values of the response ",
      "lie beyond what the arithmetic of the estimate holds, ", .xs_bounds, ". Rescaled by a ",
      "power of ten to within that, the response is scored."), fixed = TRUE)
    expect_error(select_resolution(q, criterion = "cp"), paste0(
      "select_resolution(): `cp` is NA at every level because the values of the response ",
      "lie beyond what the arithmetic of the profile holds, ", .xs_bounds, ". Rescaled by a ",
      "power of ten to within that, the response is scored."), fixed = TRUE)
    msg <- tryCatch(summary(q), error = conditionMessage)
    expect_match(msg, paste0(
      "no count criterion is finite at any level. The values of the response lie beyond what ",
      "the arithmetic of the profile holds, ", .xs_bounds, ", so no criterion that reads them ",
      "is finite: cp, loo and moran_z are NA at every level, and reliability needs a usable ",
      "variogram of them;"), fixed = TRUE)
    expect_false(grepl("constant", msg, fixed = TRUE), info = est)
  }
  # At 1e-170 the sum of squares is 0, and the print's line for one cell
  # says why rather than calling the response a constant.
  L$z <- .xs_layer()$z * 1e-170
  q <- suppressWarnings(resolution_profile(L, "z", estimate = "pooled", obs_per_cell = 4.5,
                                           sac = NULL))
  expect_identical(attr(q, "one_cell")$rss, 0)
  out <- utils::capture.output(print(q))
  expect_true(any(grepl("one cell    : rss 0, values out of range (one mean", out, fixed = TRUE)))
  # Values beyond the bounds that the test for `cp` does not find constant:
  # the warning names the estimate's columns alone.
  L$z <- .xs_layer()$z * 1e151
  r <- .xs_catch(resolution_profile(L, "z", estimate = "pooled", obs_per_cell = 4.5,
                                    sac = NULL))
  expect_false(attr(r$value, "one_cell")$constant)
  w <- Filter(function(c) inherits(c, "spatialkit_pool_unavailable"), r$warnings)
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), paste0(
    "lie beyond what the arithmetic of the estimate holds, ", .xs_bounds, "; beyond these ",
    "bounds the squares of many values can overflow or lose their precision. So `loo`, ",
    "`loo_se` and `own_weight` are NA at every level."), fixed = TRUE)
  # A response the profile's estimate finds constant keeps its own words.
  L$z <- 0.1
  q <- suppressWarnings(resolution_profile(L, "z", estimate = "pooled", obs_per_cell = 4.5,
                                           sac = NULL))
  expect_true(attr(q, "estimate")$constant)
  expect_false(.estimate_beyond_range(attr(q, "estimate")))
  msg <- tryCatch(summary(q), error = conditionMessage)
  expect_match(msg, "cp, loo and moran_z are NA at every level", fixed = TRUE)
  expect_false(grepl("beyond what the arithmetic", msg, fixed = TRUE))
})

# ---- Rounding is not spread, on the scaled column either ---------------------

# Columns that are one number in truth, stored as values that differ in
# their last bits, and one whose spread is below a unit in the last place.
.xs_rounding <- function(L) {
  N <- nrow(L)
  .xs_seed(7)
  list(mix   = sample(c(0.1, 0.3 / 3), N, TRUE),
       sqrt2 = sample(c(2, sqrt(2)^2), N, TRUE),
       feet  = sample(c(5280 * 0.3048, 5280 / (1 / 0.3048)), N, TRUE),
       third = c(1 / 3 + .Machine$double.eps / 4, rep(1 / 3, N - 1L)),
       next1 = rep(c(0.1, 0.1 + .Machine$double.eps / 16), length.out = N),
       bits  = 0.1 * (1 + .Machine$double.eps * rep_len(c(-2, 1, 0, 2, -1), N)),
       below = 1e6 + 1e-10 * L$z)
}

test_that("values that differ only by their rounding have no spread on the scaled column", {
  L <- .xs_layer(); g <- .xs_engine(L)
  cols <- .xs_rounding(L)
  for (nm in names(cols)) {
    y <- cols[[nm]]
    expect_gt(length(unique(y)), 1L)
    sc <- .spread_scaled(y)
    # The scaled column holds them as differences of whole units in the
    # last place, and their standard deviation is under the bound.
    expect_lte(stats::sd(sc$y) / sc$unit, .POOL_ULPS)
    expect_true(.constant_column(y, sc), info = nm)
    A <- .pool_fit(y, g$gi, g$K, g$xy, use_fnn = FALSE)
    expect_true(A$constant, info = nm)
    expect_true(A$flat, info = nm)
    expect_identical(A$status, "no spread inside any cell", info = nm)
    expect_identical(A$w, rep(1, g$K), info = nm)
    expect_identical(c(A$sigma2, A$tau2), c(0, 0), info = nm)
    expect_identical(A$value, A$m, info = nm)
    expect_identical(A$loo_pooled, NA_real_, info = nm)
  }
  # Through the summary: one degenerate warning, every column constant,
  # weights of 1 and no test of structure on rounding.
  for (nm in names(cols)) L[[nm]] <- cols[[nm]]
  r <- .xs_catch(summarize_by_cell(L, "z", names(cols), id_col = "cell_id",
                                   estimate = "pooled", obs_level = 0.9))
  d <- Filter(function(c) inherits(c, "spatialkit_pool_degenerate"), r$warnings)
  expect_length(d, 1L)
  expect_identical(d[[1L]]$reason, "no spread inside any cell")
  expect_identical(d[[1L]]$columns, names(cols))
  expect_identical(d[[1L]]$constant, stats::setNames(rep(TRUE, length(cols)), names(cols)))
  ns <- Filter(function(c) inherits(c, "spatialkit_no_structure"), r$warnings)
  for (w in ns) expect_false(any(names(cols) %in% sub("^pred_", "", w$columns)))
  rec <- attr(r$value, "pooled")$columns
  for (nm in names(cols)) {
    expect_identical(r$value[[paste0("..pooled_w_pred_", nm)]], rep(1, nrow(r$value)), info = nm)
    expect_true(rec[[paste0("pred_", nm)]]$constant, info = nm)
    expect_identical(rec[[paste0("pred_", nm)]]$loo, NA_real_, info = nm)
  }
})

test_that("a profile of a response that differs only by its rounding has no loo", {
  L <- .xs_layer()
  L$mix <- .xs_rounding(L)$mix
  r <- .xs_catch(resolution_profile(L, "mix", estimate = "pooled", obs_per_cell = c(18, 4.5),
                                    sac = NULL))
  p <- r$value
  e <- attr(p, "estimate")
  expect_true(e$constant)
  expect_identical(p$loo, c(NA_real_, NA_real_))
  expect_identical(p$own_weight, c(NA_real_, NA_real_))
  expect_identical(e$loo_one_cell, NA_real_)
  expect_error(select_resolution(p, criterion = "loo"),
               "`loo` is NA at every level because the response is constant", fixed = TRUE)
})

test_that("a column with spread over its rows and rounding inside its cells is flat", {
  L <- .xs_layer(); g <- .xs_engine(L)
  # 1e15 + z / 2: 4 units in the last place over the rows, 2 inside the
  # cells.  The fit is made on the scaled column with the column's own
  # readings, so nothing is fitted on the rounding inside the cells.
  y <- 1e15 + L$z / 2
  sc <- .spread_scaled(y)
  expect_gt(stats::sd(sc$y) / sc$unit, .POOL_ULPS)
  A <- .pool_fit(y, g$gi, g$K, g$xy, use_fnn = FALSE)
  expect_false(A$constant)
  expect_true(A$flat)
  expect_identical(A$status, "no spread inside any cell")
  expect_identical(A$w, rep(1, g$K))
  expect_identical(A$sigma2, 0)
  expect_identical(A$value, A$m)
  expect_gt(A$tau2, 0)
  expect_identical(A$loo_pooled, 0)
  # The warning says the values do not vary to their precision, which is
  # what it means here and for cells at 1e15 and -1e15 with a spread of 1.
  L$half <- y
  L$pm   <- ifelse(g$gi %% 2 == 0, 1e15, -1e15) + L$z
  r <- .xs_catch(summarize_by_cell(L, "half", "pm", id_col = "cell_id", estimate = "pooled"))
  d <- Filter(function(c) inherits(c, "spatialkit_pool_degenerate"), r$warnings)
  expect_length(d, 1L)
  expect_identical(conditionMessage(d[[1L]]), paste0(
    "summarize_by_cell(): the values of 'half' and 'pm' do not vary inside any cell, to the ",
    "precision of the values, so every own weight is 1 and the pooled value is the plain ",
    "mean. attr(, \"pooled\") says which columns are constant."))
  expect_identical(d[[1L]]$constant, c(half = FALSE, pm = FALSE))
})

test_that("a column constant inside large regions is flat on the scaled reading too", {
  # Two regions of 500 rows at 0.1 and 1/3.  The sums of 500 such values are
  # off by rounding, so inside the regions the sum about the means as summed
  # is above both bounds, on the raw values and on the scaled column; only
  # the sum about the corrected means (.ss_within()) shows that there is no
  # spread there.
  gi <- rep(1:2, each = 500L); n <- c(500L, 500L); N <- 1000L
  xy <- local({ .xs_seed(3); cbind(stats::runif(N, 0, 1e4), stats::runif(N, 0, 1e4)) })
  y <- c(0.1, 1 / 3)[gi]
  sc <- .spread_scaled(y)
  res <- sc$y - (as.numeric(rowsum(sc$y, gi)) / n)[gi]
  expect_false(.scaled_rounding(sum(res^2), N - 2L, sc))
  expect_true(.scaled_rounding(.ss_within(res, gi, n), N - 2L, sc))
  expect_true(.spread_readings(sc, gi, n)$flat)
  A <- .pool_fit(y, gi, 2L, xy, use_fnn = FALSE)
  expect_false(A$constant)
  expect_true(A$flat)
  expect_identical(A$status, "no spread inside any cell")
  expect_identical(A$w, c(1, 1))
  expect_identical(A$sigma2, 0)
  expect_identical(A$value, A$m)
})


# ---- The fit of the scaled column, case by case --------------------------------

test_that("a column far from 0 with spread between and inside its cells is fitted on the scaled column", {
  L <- .xs_layer(); g <- .xs_engine(L)
  # Over the rows the raw values show their spread; inside the cells they
  # lose it (4 units in the last place against a raw bound of about 18).
  y <- 1e15 + 10 * (g$gi %% 4) + L$z
  e <- y - 1e15                      # exact: within a factor 2 of 1e15
  expect_false(.no_spread(sum((y - mean(y))^2), y))
  res <- y - (as.numeric(rowsum(y, g$gi)) / g$n)[g$gi]
  expect_true(.no_spread(sum(res^2), y))
  A <- .pool_fit(y, g$gi, g$K, g$xy, use_fnn = FALSE)
  E <- .pool_fit(e, g$gi, g$K, g$xy, use_fnn = FALSE)
  expect_false(A$constant)
  expect_false(A$flat)
  expect_identical(A$status, E$status)
  expect_equal(A$w, E$w, tolerance = 1e-12)
  for (f in c("sigma2", "tau2", "Q", "loo_pooled", "loo_pooled_se", "loo_one"))
    expect_equal(A[[f]], E[[f]], tolerance = 1e-12, info = f)
  expect_lte(max(abs((A$value - 1e15) - E$value)), 2^-4)
  # Through the summary.
  L$far <- y; L$near <- e
  a <- suppressWarnings(summarize_by_cell(L, "far", id_col = "cell_id", estimate = "pooled"))
  b <- suppressWarnings(summarize_by_cell(L, "near", id_col = "cell_id", estimate = "pooled"))
  ra <- attr(a, "pooled")$columns$resp_far; rb <- attr(b, "pooled")$columns$resp_near
  expect_false(ra$constant)
  expect_identical(ra$status, rb$status)
  expect_equal(a$..pooled_w_resp_far, b$..pooled_w_resp_near, tolerance = 1e-12)
  expect_equal(ra$sigma2, rb$sigma2, tolerance = 1e-12)
  expect_lte(max(abs((a$..pooled_resp_far - 1e15) - b$..pooled_resp_near)), 2^-4)
})

test_that("one cell of values far from 0 keeps the plain mean as its value", {
  L <- .xs_layer(); g <- .xs_engine(L)
  y <- 1e15 + L$z
  one <- rep(1L, length(y))
  A <- .pool_fit(y, one, 1L, g$xy, use_fnn = FALSE)
  E <- .pool_fit(y - 1e15, one, 1L, g$xy, use_fnn = FALSE)
  expect_identical(A$status, "one cell")
  expect_false(A$constant)
  expect_identical(A$value, A$m)
  expect_identical(A$m, sum(y) / length(y))
  expect_identical(A$w, 1)
  for (f in c("sigma2", "loo_pooled", "loo_one", "loo_one_se"))
    expect_equal(A[[f]], E[[f]], tolerance = 1e-12, info = f)
})

test_that("the profile's error of one mean on raw values with spread is the raw arithmetic's", {
  L <- .xs_layer()
  .xs_seed(9)
  L$big <- 1000 + L$z
  L$cnt <- stats::rpois(nrow(L), 4)
  xy <- sf::st_coordinates(L)
  for (v in c("big", "cnt", "z")) {
    p <- suppressWarnings(resolution_profile(L, v, estimate = "pooled", obs_per_cell = 4.5,
                                             sac = NULL))
    # The rows in the order the profile reads them: by coordinates, then
    # by the response.
    y <- as.numeric(L[[v]]); n <- length(y)
    y <- y[order(xy[, 1L], xy[, 2L], y)]
    e1 <- (y - mean(y)) * n / (n - 1)
    expect_identical(attr(p, "estimate")$loo_one_cell, mean(e1^2), info = v)
    expect_identical(attr(p, "estimate")$loo_one_cell_se, stats::sd(e1^2) / sqrt(n), info = v)
  }
})


# ---- The standard errors of the leave-one-out errors at any size ---------------

test_that("the standard error of a leave-one-out error is finite inside the bounds", {
  # Inside 2^-255 to 2^255 the numbers are the formula's, bit for bit.
  .xs_seed(4)
  for (e in list(stats::rnorm(50), stats::rnorm(50) * 1e70, stats::rnorm(50) * 1e-70,
                 c(0, 0, 1), stats::rpois(30, 3) - 3))
    expect_identical(.loo_mean_square(e), c(mean(e^2), stats::sd(e^2) / sqrt(length(e))))
  expect_identical(.loo_mean_square(c(0, 0)), c(0, 0))
  # Beyond them the squares of the squares overflow or underflow; the
  # numbers are those of the residuals scaled by a power of two.
  e <- stats::rnorm(50)
  ref <- c(mean(e^2), stats::sd(e^2) / sqrt(50))
  for (k in c(80, 140, -80, -140)) {
    got <- .loo_mean_square(e * 10^k)
    expect_true(all(is.finite(got)) && all(got > 0), info = k)
    expect_equal(got / 10^(2 * k), ref, tolerance = 1e-12, info = k)
  }
  # Through the summary and the profile, at 1e80 and 1e-100.
  L <- .xs_layer()
  L$up <- L$z * 1e80; L$down <- L$z * 1e-100
  r <- suppressWarnings(summarize_by_cell(L, "z", c("up", "down"), id_col = "cell_id",
                                          estimate = "pooled"))
  rec <- attr(r, "pooled")$columns
  expect_equal(rec$pred_up$loo_se / 1e160, rec$resp_z$loo_se, tolerance = 1e-12)
  expect_equal(rec$pred_down$loo_se / 1e-200, rec$resp_z$loo_se, tolerance = 1e-12)
  p  <- suppressWarnings(resolution_profile(L, "up", estimate = "pooled",
                                            obs_per_cell = c(18, 4.5), sac = NULL))
  p0 <- suppressWarnings(resolution_profile(L, "z", estimate = "pooled",
                                            obs_per_cell = c(18, 4.5), sac = NULL))
  expect_true(all(is.finite(p$loo_se)))
  expect_equal(p$loo_se / 1e160, p0$loo_se, tolerance = 1e-12)
  expect_equal(attr(p, "estimate")$loo_one_cell_se / 1e160,
               attr(p0, "estimate")$loo_one_cell_se, tolerance = 1e-12)
})
