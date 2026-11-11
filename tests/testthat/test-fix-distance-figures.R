# tests/testthat/test-fix-distance-figures.R
# ---------------------------------------------------------------------------
# A distance in CRS units is printed with four significant digits in the
# warnings, log lines, prints and plot captions that quote it, so the figure
# reads the same in any unit.  With a fixed number of decimals (%.0f, %.1f) a
# layer in a CRS whose unit is a large length read "fitted range (0) exceeds
# the largest lag the variogram was fitted over (0.0006598 = 1 x cutoff 0)"
# (BUG-256).  Each test runs one layer in metres and the same layer in UTM
# with +to_meter=1000000, and reads the figures back from the text.
# ---------------------------------------------------------------------------

.df_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.df_big_unit <- "+proj=utm +zone=17 +datum=WGS84 +to_meter=1000000 +no_defs"

# The layer in metres (EPSG:32617) or in units of 1000 km.
.df_layer <- function(xy, unit) {
  L <- sf::st_as_sf(transform(xy, x = 5e5 + x, y = 4e6 + y), coords = c("x", "y"),
                    crs = 32617)
  if (unit == "m") L else sf::st_transform(L, .df_big_unit)
}

# 300 points on a 1 km square with a trend and no sill.
.df_trend <- function() {
  .df_seed(5)
  xy <- data.frame(x = stats::runif(300, 0, 1000), y = stats::runif(300, 0, 1000))
  xy$t <- xy$x / 50 + stats::rnorm(300)
  xy
}

# 300 points of an exponential field (scale 60 m) plus 0.5 a plus noise.
.df_field <- function() {
  .df_seed(7)
  xy <- data.frame(x = stats::runif(300, 0, 1000), y = stats::runif(300, 0, 1000))
  D <- as.matrix(stats::dist(xy))
  xy$a <- stats::rnorm(300)
  xy$z <- as.numeric(t(chol(exp(-D / 60) + diag(1e-8, 300))) %*% stats::rnorm(300)) +
    0.5 * xy$a + stats::rnorm(300, sd = 0.2)
  xy
}

# A range given by hand, as resolution_profile() and the sweep take one.
.df_sac <- function(L, range) {
  structure(range, class = c("sac_range", "numeric"),
            variogram_model = data.frame(model = c("Nug", "Exp"), psill = c(0.3, 1),
                                         range = c(0, range / 3)),
            crs = sf::st_crs(L))
}

# The figures a text quotes after each fixed lead-in, as numbers.
.df_figures <- function(text, leads) {
  vapply(leads, function(lead) {
    m <- regmatches(text, regexpr(paste0(lead, "[-+.0-9e]+"), text, perl = TRUE))
    if (!length(m)) return(NA_real_)
    # A figure that ends a sentence keeps its full stop out.
    as.numeric(sub("\\.$", "", sub(lead, "", m, perl = TRUE)))
  }, numeric(1))
}

# Each figure is the value to four significant digits, and the figure in the
# unit of 1000 km is the figure in metres times 1e-6 to four significant
# digits.  The values come from variogram fits, whose last bits move with
# the platform, so they are read from the result rather than pinned; the
# relative tolerance of 1e-3 is one unit in the fourth digit, which a value
# lying on a rounding boundary can move in either unit.
.df_expect_four_digits <- function(fig, value, info) {
  expect_false(anyNA(fig), info = info)
  expect_equal(unname(fig), unname(signif(value, 4)), tolerance = 1e-9, info = info)
}
.df_expect_same <- function(fig_m, fig_u, info) {
  expect_equal(unname(fig_u) / 1e-6, unname(fig_m), tolerance = 1e-3, info = info)
}

test_that("estimate_sac_range()'s over-cutoff warning and plot() quote the range and the lag in any unit", {
  skip_if_not_installed("gstat")
  xy <- .df_trend()
  leads <- c(range = "fitted range \\(", lag = "fitted over \\(",
             cutoff = "x cutoff ")
  out <- list()
  for (unit in c("m", "1e6")) {
    L <- .df_layer(xy, unit)
    r <- NULL
    log <- capture_spatialkit_log(r <- estimate_sac_range(L, "t"), level = logger::WARN)
    line <- grep("estimate_sac_range(): fitted range", log, value = TRUE, fixed = TRUE)
    expect_length(line, 1L)
    expect_true(is.na(r))
    expect_identical(attr(r, "rejected_reason"), "fitted range exceeds the largest lag fitted")
    # The old text read "fitted range (0) ... (0.0006598 = 1 x cutoff 0)".
    expect_false(grepl("(0)", line, fixed = TRUE), info = unit)
    expect_false(grepl("cutoff 0)", line, fixed = TRUE), info = unit)
    fig <- .df_figures(line, leads)
    .df_expect_four_digits(fig, c(attr(r, "rejected_range"), attr(r, "cutoff_dist"),
                                  attr(r, "cutoff_dist")), info = unit)
    cap <- plot(r)$labels$subtitle
    cap <- gsub("\n", " ", cap, fixed = TRUE)
    fig_p <- .df_figures(cap, c(range = "the fitted range \\(", lag = "largest lag fitted \\("))
    .df_expect_four_digits(fig_p, c(attr(r, "rejected_range"), attr(r, "cutoff_dist")),
                           info = paste(unit, "plot"))
    out[[unit]] <- c(fig, fig_p)
  }
  .df_expect_same(out$m, out[["1e6"]], "over-cutoff")
  # The figures of the issue: 2639 and 0.002639, 659.8 and 0.0006598.
  expect_gt(out$m[["range"]], 1000); expect_lt(out$m[["range"]], 1e4)
  expect_lt(out[["1e6"]][["range"]], 0.01)
})

test_that("an identified range is quoted in any unit by estimate_sac_range(), plot() and make_folds()", {
  skip_if_not_installed("gstat")
  xy <- .df_field()
  out <- list()
  for (unit in c("m", "1e6")) {
    L <- .df_layer(xy, unit)
    r <- NULL
    log <- capture_spatialkit_log(r <- estimate_sac_range(L, "z"), level = logger::INFO)
    rng <- as.numeric(r)
    expect_true(is.finite(rng), info = unit)
    # The directional line: the all-directions estimate it names.
    dir_line <- grep("directional ranges vary", log, value = TRUE, fixed = TRUE)
    if (length(dir_line)) {
      f <- .df_figures(dir_line, c(all = "all-directions estimate \\("))
      .df_expect_four_digits(f, rng, info = paste(unit, "directional"))
      expect_false(grepl("= 0.0[,)]", dir_line), info = unit)
    }
    sub <- plot(r)$labels$subtitle
    f_plot <- .df_figures(sub, c(range = "Effective range = "))
    .df_expect_four_digits(f_plot, rng, info = paste(unit, "subtitle"))

    # make_folds(): the range it reads, a block_size below it, and blocks
    # too small for it.  The messages quote the range make_folds() estimated.
    msgs <- character(0); warns <- character(0)
    f <- withCallingHandlers(
      make_folds(L, k = 5, method = "block_kfold", response_var = "z",
                 auto_range = TRUE, block_size = rng / 2, seed = 1),
      message = function(c) { msgs <<- c(msgs, conditionMessage(c)); invokeRestart("muffleMessage") },
      warning = function(c) { warns <<- c(warns, conditionMessage(c)); invokeRestart("muffleWarning") })
    r_mf <- as.numeric(f$params$sac_range)
    expect_true(is.finite(r_mf), info = unit)
    m_line <- grep("estimated spatial autocorrelation range = ", msgs, value = TRUE, fixed = TRUE)
    expect_length(m_line, 1L)
    f_msg <- .df_figures(m_line, c(range = "autocorrelation range = "))
    .df_expect_four_digits(f_msg, r_mf, info = paste(unit, "message"))
    w_line <- grep("block_size (", warns, value = TRUE, fixed = TRUE)
    expect_length(w_line, 1L)
    f_warn <- .df_figures(w_line, c(size = "block_size \\(", range = "autocorrelation range \\("))
    .df_expect_four_digits(f_warn, c(rng / 2, r_mf), info = paste(unit, "block_size warning"))
    # The old text read "block_size (0.0) < estimated autocorrelation range (0.0)".
    expect_false(grepl("(0.0)", w_line, fixed = TRUE), info = unit)

    # Supplied blocks smaller than the range: their scale, the range, and a
    # size at least the range.
    warns <- character(0)
    blk <- sf::st_sf(geometry = sf::st_make_grid(L, n = c(12, 12)))
    f2 <- withCallingHandlers(
      suppressMessages(make_folds(L, k = 5, method = "block_kfold", response_var = "z",
                                  blocks = blk, seed = 1)),
      warning = function(c) { warns <<- c(warns, conditionMessage(c)); invokeRestart("muffleWarning") })
    b_line <- grep("median scale", warns, value = TRUE, fixed = TRUE)
    expect_length(b_line, 1L)
    r_b <- as.numeric(f2$params$sac_range)
    f_blk <- .df_figures(b_line, c(scale = "median scale ", range = "autocorrelation range \\(",
                                   across = "supply blocks at least "))
    .df_expect_four_digits(f_blk[c("range")], r_b, info = paste(unit, "blocks"))
    expect_gte(f_blk[["across"]], r_b)
    expect_lt(f_blk[["across"]], r_b * 1.01)
    expect_lt(f_blk[["scale"]], r_b)
    expect_gt(f_blk[["scale"]], 0)
    out[[unit]] <- c(rng = rng, plot = f_plot, msg = f_msg, warn = f_warn, blk = f_blk)
  }
  .df_expect_same(out$m, out[["1e6"]], "identified range")
})

test_that("resolution_profile(), cv_block_size_sweep() and kriging_adequacy() print the range in any unit", {
  xy <- .df_field()
  out <- list()
  for (unit in c("m", "1e6")) {
    s <- if (unit == "m") 1 else 1e-6
    L <- .df_layer(xy, unit)
    sac <- .df_sac(L, 240 * s)
    p <- suppressWarnings(suppressMessages(
      resolution_profile(L, "z", sac = sac, n_levels = 6, nstart = 2)))
    pr <- paste(utils::capture.output(print(p)), collapse = "\n")
    f_rp <- .df_figures(pr, c(ladder = "from range ", vg = "partial sill 1, range "))
    .df_expect_four_digits(f_rp, c(240, 240) * s, info = paste(unit, "resolution_profile"))
    # The old print read "floor 16 from range 0".
    expect_false(grepl("from range 0,", pr, fixed = TRUE), info = unit)

    # A range the data cannot support: the warning names it.
    sac20 <- .df_sac(L, 20 * s)
    w <- character(0)
    withCallingHandlers(
      suppressMessages(resolution_profile(L, "z", sac = sac20, n_levels = 6, nstart = 2)),
      warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
    w <- grep("cells no wider than the autocorrelation range", w, value = TRUE, fixed = TRUE)
    expect_length(w, 1L)
    f_ns <- .df_figures(w, c(range = "autocorrelation range \\("))
    .df_expect_four_digits(f_ns, 20 * s, info = paste(unit, "not supported"))

    sw <- suppressWarnings(suppressMessages(cv_block_size_sweep(
      L, "z", "a", fit_fn = function(tr) lm_spatial_fit(tr, "z", "a"), k = 3,
      n_sizes = 3, sac = sac, quiet = TRUE)))
    pr <- paste(utils::capture.output(print(sw)), collapse = "\n")
    f_sw <- .df_figures(pr, c(range = "estimated autocorrelation range: "))
    .df_expect_four_digits(f_sw, 240 * s, info = paste(unit, "sweep print"))
    f_swp <- .df_figures(plot(sw)$labels$subtitle, c(range = "autocorrelation range "))
    .df_expect_four_digits(f_swp, 240 * s, info = paste(unit, "sweep plot"))
    out[[unit]] <- c(rp = f_rp, ns = f_ns, sw = f_sw, swp = f_swp)
  }
  .df_expect_same(out$m, out[["1e6"]], "prints")
})

test_that("kriging_adequacy()'s print quotes an identified range in any unit", {
  skip_if_not_installed("gstat")
  xy <- .df_field()
  out <- list()
  for (unit in c("m", "1e6")) {
    L <- .df_layer(xy, unit)
    r <- suppressWarnings(estimate_sac_range(L, "z"))
    expect_true(is.finite(r), info = unit)
    cells <- sf::st_sf(poly_id = 1:36, geometry = sf::st_make_grid(L, n = c(6, 6)))
    a <- suppressWarnings(suppressMessages(assign_features_to_polygons(L, cells)))
    ka <- suppressWarnings(suppressMessages(kriging_adequacy(a, "z", cells, sac = r)))
    pr <- grep("^  variogram: ", withr::with_options(list(width = 10000),
                                                     utils::capture.output(print(ka))),
               value = TRUE)
    expect_length(pr, 1L)
    f <- .df_figures(pr, c(range = "%\\), range "))
    expect_true(isTRUE(attr(ka, "range_identified")), info = unit)
    .df_expect_four_digits(f, attr(ka, "range"), info = unit)
    out[[unit]] <- f
  }
  .df_expect_same(out$m, out[["1e6"]], "kriging_adequacy print")
})

# Every log line and warning a call emits, without the logger's prefix.
.df_text <- function(expr) {
  w <- character(0)
  lines <- capture_spatialkit_log(withCallingHandlers(expr,
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
    message = function(c) invokeRestart("muffleMessage")), level = logger::INFO)
  unique(c(sub("^[A-Z]+ \\[[^]]*\\] ", "", lines), w))
}

# The one text that contains `key`, and the figure after `lead` in it.
.df_one <- function(text, key, lead, info) {
  line <- grep(key, text, value = TRUE, fixed = TRUE)
  expect_length(line, 1L)
  if (length(line) != 1L) return(NA_real_)
  f <- unname(.df_figures(line, lead))
  expect_true(is.finite(f) && f > 0, info = info)
  f
}

test_that("the hole-effect warning and the REML line of estimate_sac_range() quote the range in any unit", {
  skip_if_not_installed("gstat")
  # A sine of period 250 m along x on a jittered 15 x 15 grid, without
  # random draws: its variogram falls over the shorter half of the lags by
  # 0.40 of its mean, well past the 0.15 that calls it decreasing.
  g <- expand.grid(i = 0:14, j = 0:14)
  id <- seq_len(nrow(g))
  xy <- data.frame(x = g$i * 1000 / 14 + ((id * 37) %% 11 - 5),
                   y = g$j * 1000 / 14 + ((id * 53) %% 13 - 6))
  xy$z <- sin(2 * pi * xy$x / 250) + 0.2 * (((id * 7919) %% 101) / 101 - 0.5)
  fld <- .df_field()
  out <- list()
  for (unit in c("m", "1e6")) {
    r <- NULL
    txt <- .df_text(r <- estimate_sac_range(.df_layer(xy, unit), "z", lags = "equal"))
    expect_identical(attr(r, "rejected_reason"), "empirical variogram decreases with distance",
                     info = unit)
    hole <- .df_one(txt, "the fit reported", "the fit reported ", paste(unit, "hole"))
    # The old text read "(the fit reported 0)" in the unit of 1000 km.
    .df_expect_four_digits(hole, attr(r, "rejected_range"), info = paste(unit, "hole"))
    reml <- .df_one(.df_text(estimate_sac_range(.df_layer(fld, unit), "z", predictor_vars = "a",
                                                detrend = "reml")),
                    "REML detrending", "range parameter ", paste(unit, "REML"))
    out[[unit]] <- c(hole = hole, reml = reml)
  }
  .df_expect_same(out$m, out[["1e6"]], "hole effect and REML")
})

test_that("make_folds() quotes block dimensions, distances outside blocks and nndm medians in any unit", {
  skip_if_not_installed("gstat")
  xy <- .df_field()
  out <- list()
  for (unit in c("m", "1e6")) {
    L <- .df_layer(xy, unit)
    f <- NULL
    txt <- .df_text(f <- make_folds(L, k = 5, method = "block_kfold", response_var = "z",
                                    block_multiplier = 20, seed = 1))
    dim_f <- .df_one(txt, "block dimension (", c(dim = "block dimension \\("), paste(unit, "dimension"))
    dim_r <- .df_one(txt, "block dimension (", c(range = "< autocorrelation range \\("),
                     paste(unit, "dimension range"))
    .df_expect_four_digits(dim_r, as.numeric(f$params$sac_range), info = paste(unit, "dimension"))
    # The old text read "block dimension (0.0) < autocorrelation range (0.0)".
    txt <- .df_text(f <- make_folds(L, k = 5, method = "block_kfold", response_var = "z",
                                    block_nx = 12, block_ny = 12, seed = 1))
    nx_r <- .df_one(txt, "block_nx/block_ny yield", c(range = "autocorrelation range \\("),
                    paste(unit, "block_nx"))
    .df_expect_four_digits(nx_r, as.numeric(f$params$sac_range), info = paste(unit, "block_nx"))
    # Blocks over half the points: the rest fall outside them, the farthest
    # by tens of metres.
    blk <- sf::st_sf(geometry = sf::st_make_grid(L[1:150, ], n = c(4, 4)))
    far <- .df_one(.df_text(make_folds(L, k = 5, method = "block_kfold", blocks = blk, seed = 1)),
                   "outside every supplied block", c(far = "the farthest by "), paste(unit, "outside"))
    nn <- .df_one(.df_text(make_folds(L, k = 5, method = "nndm", seed = 1,
                                      prediction_points = L[1:50, ])),
                  "folds. Target prediction-to-training", c(med = "Realised test-to-training distance: median "),
                  paste(unit, "nndm"))
    kn <- .df_one(.df_text(make_folds(L, k = 5, method = "knndm", seed = 1,
                                      prediction_points = L[1:50, ])),
                  "Held-out-to-training distance", c(med = "Held-out-to-training distance: median "),
                  paste(unit, "knndm"))
    out[[unit]] <- c(dim = dim_f, range = dim_r, nx = nx_r, far = far, nndm = nn, knndm = kn)
  }
  .df_expect_same(out$m, out[["1e6"]], "make_folds")
})

test_that("the GWR fallback bandwidth is quoted in any unit, and an adaptive one as a count", {
  skip_if_not_installed("sp")
  # The fallback is reached only when the bandwidth search fails, so the
  # search is made to fail.  A fixed fallback is a third of the bounding
  # box's diagonal in CRS units: 46.15 on this 100 m square in metres and
  # 4.615e-05 in units of 1000 km, where "%.4f" printed 0.0000 (BUG-273).
  local_mocked_bindings(.gwr_bw_search = function(...) stop("forced failure"),
                        .package = "spatialkit")
  .df_seed(1)
  xy <- data.frame(x = stats::runif(80, 0, 100), y = stats::runif(80, 0, 100),
                   w = stats::rnorm(80), v = stats::rnorm(80))
  xy$yy <- 1 + 2 * xy$w - xy$v + stats::rnorm(80, sd = 0.3)
  out <- list()
  for (unit in c("m", "1e6")) {
    L <- .df_layer(xy, unit)
    f <- NULL
    txt <- .df_text(f <- fit_gwr_model(L, "yy", c("w", "v"), adaptive = FALSE))
    expect_true(isTRUE(f$info$bandwidth_is_fallback))
    bw <- f$info$bandwidth
    fit_f <- .df_one(txt, "arbitrary fallback bandwidth", c(bw = "fallback bandwidth = "),
                     paste(unit, "fit"))
    cmp_f <- .df_one(.df_text(compare_models(list(G = f))), "used a fallback bandwidth",
                     c(bw = "used a fallback bandwidth \\("), paste(unit, "compare"))
    sel <- NULL
    txt_s <- .df_text(sel <- gwr_model_selection(L, "yy", c("w", "v"), adaptive = FALSE))
    sel_f <- .df_one(txt_s, "arbitrary fallback bandwidth", c(bw = "fallback bandwidth = "),
                     paste(unit, "selection"))
    expect_identical(sel$bandwidth_source, "fallback")
    # The bounding box's diagonal does not come from a fit, and 46.1534
    # lies far from a rounding boundary in its fourth digit, so the last
    # bits the projection moves between platforms leave the four digits
    # alone: the read-back is exact against signif().
    for (fig in list(fit_f, cmp_f, sel_f))
      .df_expect_four_digits(fig, bw, info = unit)
    .df_expect_four_digits(sel_f, sel$bandwidth, info = paste(unit, "selection value"))
    expect_false(any(grepl("0.0000", c(txt, txt_s), fixed = TRUE)), info = unit)
    out[[unit]] <- c(fit = fit_f, compare = cmp_f, selection = sel_f)
  }
  expect_equal(unname(out$m), rep(46.15, 3))
  expect_equal(unname(out[["1e6"]]), rep(4.615e-05, 3))
  .df_expect_same(out$m, out[["1e6"]], "fallback bandwidth")

  # Adaptive: a whole count of neighbours, 50 here, and on six rows the
  # count of 10 clamped to 6, with both counts whole.
  L <- .df_layer(xy, "m")
  f <- NULL
  txt <- .df_text(f <- fit_gwr_model(L, "yy", c("w", "v"), adaptive = TRUE))
  expect_equal(f$info$bandwidth, 50)
  line <- grep("arbitrary fallback bandwidth", txt, value = TRUE, fixed = TRUE)
  expect_match(line, "fallback bandwidth = 50. ", fixed = TRUE)
  expect_match(.df_text(compare_models(list(G = f))), "used a fallback bandwidth (50)",
               fixed = TRUE, all = FALSE)
  expect_match(.df_text(gwr_model_selection(L, "yy", c("w", "v"), adaptive = TRUE)),
               "fallback bandwidth = 50. ", fixed = TRUE, all = FALSE)
  f6 <- NULL
  txt <- .df_text(f6 <- fit_gwr_model(L[1:6, ], "yy", c("w", "v"), adaptive = TRUE))
  expect_equal(f6$info$bandwidth, 6)
  expect_match(txt, "fallback bandwidth = 6 (derived as 10, then clamped", fixed = TRUE,
               all = FALSE)
  expect_false(any(grepl(".0000", txt, fixed = TRUE)))
})

test_that(".fmt_bw() prints a whole adaptive count in full and a fixed bandwidth to four digits", {
  # A reachable adaptive fallback is below 10000, where four significant
  # digits print the same count, so a larger count pins the branch: as a
  # count of neighbours 12345 prints whole, as a length it rounds to 12340.
  expect_identical(.fmt_bw(12345, TRUE), "12345")
  expect_identical(.fmt_bw(12345, FALSE), "12340")
  expect_identical(.fmt_bw(123456789, TRUE), "123456789")
  expect_identical(.fmt_bw(4.61534e-05, FALSE), "4.615e-05")
  # A count that is not whole, or a missing flag, prints as a length.
  expect_identical(.fmt_bw(10.25, TRUE), "10.25")
  expect_identical(.fmt_bw(46.1534, NULL), "46.15")
})
