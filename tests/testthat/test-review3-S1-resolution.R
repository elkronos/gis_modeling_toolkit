# tests/testthat/test-review3-S1-resolution.R
# ---------------------------------------------------------------------------
# Regressions from the third review of the resolution step:
# determine_optimal_levels(), resolution_profile() and their messages.
# ---------------------------------------------------------------------------

# Every R warning an expression raises, muffled.
r3_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

# `reps` visits to each location in `xy` (metres, offset into UTM 32N).
r3_visits <- function(xy, reps, z = NULL) {
  i <- rep(seq_len(nrow(xy)), each = reps)
  d <- data.frame(x = 5e5 + xy[i, 1], y = 5e6 + xy[i, 2])
  if (!is.null(z)) d$z <- z
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}

r3_stations <- cbind(c(0, 1000, 5000, 9000, 9500), c(0, 8000, 2000, 500, 9000))

# A hand-made sac_range, so the variogram-based columns exist without gstat.
r3_sac <- function(range, nugget = 1, psill = 1, ...) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, if (is.finite(range)) range / 3 else 100),
                   stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(32632), ...)
}

r3_pts <- function(n = 200, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  d$z <- rnorm(n)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}


test_that("a WSS that falls to zero is the elbow when the rest has none", {
  eb <- spatialkit:::.elbow_from_wss
  # A curve like c / k to k = 4, then zero: one cell per location at k = 5.
  w <- c(1000 / (1:4), 0)
  out <- eb(w)
  expect_true(out$structured)
  expect_identical(out$knee_k, 5L)
  # Zero is relative: k-means leaves floating-point residue where the exact
  # value is 0, and the log of that dragged the line down as log(0) would.
  expect_identical(eb(c(1000 / (1:4), 7.8e-17))$knee_k, 5L)
  # Two levels are too few for a line, but not for a fall to zero.
  two <- eb(c(10, 0))
  expect_true(two$structured)
  expect_identical(two$knee_k, 2L)
  # Without a zero nothing changes: two levels, no elbow, the midpoint.
  expect_false(eb(c(5, 3))$structured)
  expect_equal(eb(c(5, 3))$knee_k, 1)
  # An elbow in the positive part keeps its place ahead of a later zero.
  bent <- eb(c(1000, 400, 150, 60 * 4 / (4:12), 0))
  expect_true(bent$structured)
  expect_identical(bent$knee_k, 4L)
})


test_that("repeat visits to five stations give one cell per station", {
  # Five stations visited thirty times each: the WSS reaches 0 at k = 5, and
  # the level used to be dropped from the elbow, which then said "no cluster
  # structure" and answered 3 (1 m of jitter answered 5).
  pts <- r3_visits(r3_stations, 30)
  out <- r3_warnings(determine_optimal_levels(pts, max_levels = 12))
  expect_identical(out$value[1L], 5L)
  expect_false(any(grepl("no elbow", out$warnings)))
  # The same with coordinates that are not whole metres.
  frac <- r3_visits(r3_stations + 0.123456, 30)
  expect_identical(suppressWarnings(determine_optimal_levels(frac, max_levels = 12))[1L], 5L)
  # Two stations visited thirty times each: 2, not the 1 a two-level ladder gave.
  expect_identical(suppressWarnings(determine_optimal_levels(r3_visits(r3_stations[1:2, ], 30)))[1L],
                   2L)
  # The profile's elbow column names the same count.
  prof <- resolution_profile(pts, min_cell_n = 1, n_levels = 4, nstart = 3)
  expect_identical(prof$levels, 2:5)
  expect_true(is.finite(prof$elbow[prof$levels == 5L]))
  expect_identical(select_resolution(prof, "elbow")$best, 5L)
})


test_that("two clusters of repeat-visited stations keep their elbow at 2", {
  # Twenty stations in two groups, ten visits each: the WSS at k = 20 is
  # floating-point residue (7.8e-17), which used to drag the whole log-log
  # line down and, at max_levels = 30, report no cluster structure.
  set.seed(2)
  st <- rbind(cbind(runif(10, 0, 500), runif(10, 0, 500)),
              cbind(runif(10, 8000, 8500), runif(10, 8000, 8500)))
  pts <- r3_visits(st, 10)
  for (ml in c(12, 30)) {
    out <- r3_warnings(determine_optimal_levels(pts, max_levels = ml))
    expect_identical(out$value[1L], 2L, info = ml)
    expect_false(any(grepl("no elbow", out$warnings)), info = ml)
  }
})


test_that("the no-elbow warnings name the bound that ended the ladder", {
  set.seed(1)
  two <- sf::st_as_sf(
    data.frame(x = 5e5 + c(runif(25, 0, 10), runif(25, 90, 100)),
               y = 5e6 + c(runif(25, 0, 10), runif(25, 90, 100))),
    coords = c("x", "y"), crs = 32632)
  # A two-level ladder has no line to test: it used to be told it fell in a
  # straight line "as it does for points with no cluster structure".
  out <- r3_warnings(determine_optimal_levels(two, max_levels = 2))
  expect_true(any(grepl("too short to read an elbow.*raise max_levels", out$warnings)))
  expect_false(any(grepl("no cluster structure", out$warnings)))
  out <- r3_warnings(determine_optimal_levels(two, max_levels = 1))
  expect_true(any(grepl("max_levels = 1, raised to 2", out$warnings)))
  # With three levels the elbow is there.
  expect_length(r3_warnings(determine_optimal_levels(two, max_levels = 3))$warnings, 0L)
  # A ladder ended by the points, not by max_levels, says so.
  set.seed(5)
  few <- sf::st_as_sf(data.frame(x = 5e5 + runif(12, 0, 1000), y = 5e6 + runif(12, 0, 1000)),
                      coords = c("x", "y"), crs = 32632)
  out <- r3_warnings(determine_optimal_levels(few, max_levels = 40))
  expect_true(any(grepl("k = 1 to 11, one short of the 12 points", out$warnings)))
  expect_false(any(grepl("max_levels\\)", out$warnings)))
})


test_that("combined returns the geometric ranking when the elbow is below ten cells", {
  # Six separated clusters: the elbow is 6, a count Moran's I cannot score
  # (the nine-cell floor).  Ranking the window anyway put the smallest count
  # it scores, ten, first whatever the response did.
  set.seed(3)
  K <- 6; n <- 360
  cen <- cbind(c(0, 6000, 12000, 0, 6000, 12000), c(0, 0, 0, 7000, 7000, 7000))
  cl  <- rep(seq_len(K), length.out = n)
  d <- data.frame(x = 5e5 + cen[cl, 1] + rnorm(n, 0, 250),
                  y = 5e6 + cen[cl, 2] + rnorm(n, 0, 250), p = rnorm(n))
  d$z <- d$p + rnorm(n)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  geo <- determine_optimal_levels(pts, max_levels = 30)
  expect_identical(geo[1L], 6L)
  lines <- capture_spatialkit_log(
    k <- suppressWarnings(determine_optimal_levels(pts, max_levels = 30,
                                                   response_var = "z", predictor_vars = "p")))
  dg <- attr(k, "diagnostics")
  expect_identical(dg$knee_k, 6L)
  expect_identical(as.integer(k), as.integer(geo))
  expect_identical(dg$criterion, "geometric")
  expect_match(dg$fallback, "below the ten cells Moran's I needs")
  expect_null(dg$combined_rank)
  expect_true(log_has(lines, "the WSS elbow is at k = 6, below the ten cells"))
})


test_that("combined ranks geometry on the log-log sag and unscored candidates last", {
  # Twelve separated clusters: the elbow is 12, its window 8 to 16, and 8
  # and 9 sit below the nine-cell floor.  The linear chord across the window
  # used to rank a k next to the window's own middle first.
  set.seed(1)
  K <- 12; n <- 720
  cen <- cbind(rep(c(0, 6000, 12000, 18000), 3), rep(c(0, 7000, 14000), each = 4))
  cl  <- rep(seq_len(K), length.out = n)
  d <- data.frame(x = 5e5 + cen[cl, 1] + rnorm(n, 0, 250),
                  y = 5e6 + cen[cl, 2] + rnorm(n, 0, 250), p = rnorm(n))
  d$z <- d$p + rnorm(n)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  expect_identical(determine_optimal_levels(pts, max_levels = 30)[1L], 12L)
  k  <- suppressWarnings(determine_optimal_levels(pts, max_levels = 30,
                                                  response_var = "z", predictor_vars = "p"))
  dg <- attr(k, "diagnostics")
  expect_identical(dg$knee_k, 12L)
  expect_identical(dg$criterion, "combined")
  cr <- dg$combined_rank
  ks <- as.integer(names(cr))
  scored <- is.finite(dg$moran_z[ks])
  expect_true(any(scored) && any(!scored))
  # Every unscored candidate takes the last place on the Moran axis, so a
  # scored one comes first ...
  expect_true(k[1L] %in% ks[scored])
  # ... and the unscored rank behind the elbow on the rank average.
  expect_true(all(cr[!scored] > cr[ks == 12L]))
})


test_that("combined on points with no elbow is ordered by Moran's z alone", {
  set.seed(7)
  n <- 800
  d <- data.frame(x = 5e5 + runif(n, 0, 1e4), y = 5e6 + runif(n, 0, 1e4), p = rnorm(n))
  d$z <- d$p + rnorm(n)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  k  <- suppressWarnings(determine_optimal_levels(pts, max_levels = 40,
                                                  response_var = "z", predictor_vars = "p"))
  dg <- attr(k, "diagnostics")
  cr <- dg$combined_rank
  ks <- as.integer(names(cr))
  scored <- is.finite(dg$moran_z[ks])
  expect_true(any(scored))
  # A flat geometric axis: every unscored candidate is tied ...
  expect_length(unique(cr[!scored]), 1L)
  # ... behind every scored one, which come in order of |z|.
  expect_true(max(cr[scored]) < min(cr[!scored]))
  sk <- ks[scored]
  expect_identical(k[seq_along(sk)], sk[order(abs(dg$moran_z[sk]))])
  # Ties go to the candidates nearest the elbow the window was drawn around.
  expect_identical(k[length(sk) + 1L], dg$knee_k)
})


test_that("resolution_profile() drops points with no coordinates with an R warning", {
  pts <- r3_pts(60)
  holed <- rbind(pts[1:10, ], sf::st_sf(z = 0, geometry = sf::st_sfc(sf::st_point(), crs = 32632)),
                 pts[11:60, ])
  out <- r3_warnings(resolution_profile(holed, n_levels = 3, nstart = 2))
  expect_true(any(grepl("resolution_profile\\(\\): dropping 1 point", out$warnings)))
  expect_identical(attr(out$value, "bounds")$n, 60L)
})


test_that("a range below the shortest lag gives cp no nugget and reliability nothing", {
  # That refusal means the structure cannot be told from a nugget, so the
  # nugget is not identified either; on white noise detrended by REML it was
  # 6e-7 on a sill of 0.99 and Cp ran to the support ceiling.  Cp takes its
  # noise variance from the finest level instead.
  pts <- r3_pts(200)
  rej <- r3_sac(NA_real_, nugget = 6e-7, psill = 0.99)
  attr(rej, "rejected_reason") <- "fitted range is below the shortest lag fitted"
  out <- r3_warnings(resolution_profile(pts, "z", sac = rej, n_levels = 4, nstart = 2))
  expect_true(any(grepl("below the shortest lag fitted.*not identified", out$warnings)))
  expect_true(all(is.na(out$value$reliability)))
  expect_true(all(is.finite(out$value$cp)))
  expect_identical(attr(out$value, "cp_noise")$source, "finest-level residual mean square")
  expect_null(attr(out$value, "variogram"))
  # The other refusals keep their rule: a range past the lags still gives Cp
  # its nugget.
  attr(rej, "rejected_reason") <- "fitted range exceeds the largest lag fitted"
  attr(rej, "variogram_model")$psill[1] <- 0.5
  out <- r3_warnings(resolution_profile(pts, "z", sac = rej, n_levels = 4, nstart = 2))
  expect_true(all(is.finite(out$value$cp)))
})


test_that("a nugget next to zero is warned about like a nugget of zero", {
  pts <- r3_pts(200)
  out <- r3_warnings(resolution_profile(pts, "z", sac = r3_sac(300, nugget = 6e-7, psill = 0.99),
                                        n_levels = 4, nstart = 2))
  expect_true(any(grepl("nugget is 0 \\(6e-07 against a sill of 0.99\\)", out$warnings)))
  # A nugget that is a real share of the sill raises nothing.
  expect_length(r3_warnings(resolution_profile(pts, "z", sac = r3_sac(300, nugget = 0.01),
                                               n_levels = 4, nstart = 2))$warnings, 0L)
})


test_that("warnings about the variogram estimated here do not blame a `sac` argument", {
  skip_if_not_installed("gstat")
  pts <- r3_pts(200)
  rej <- r3_sac(NA_real_, nugget = 0.5)
  attr(rej, "rejected_reason") <- "fitted range exceeds the largest lag fitted"
  local_mocked_bindings(estimate_sac_range = function(...) rej)
  out <- r3_warnings(resolution_profile(pts, "z", n_levels = 4, nstart = 2))
  expect_true(any(grepl("the variogram estimated from `response_var`", out$warnings,
                        fixed = TRUE)))
  expect_false(any(grepl("`sac` reports", out$warnings, fixed = TRUE)))
  # A supplied one is still called `sac`.
  out <- r3_warnings(resolution_profile(pts, "z", sac = rej, n_levels = 4, nstart = 2))
  expect_true(any(grepl("`sac` reports no usable range", out$warnings, fixed = TRUE)))
})


test_that("a ceiling one short of the points is not credited to the distinct locations", {
  pts <- r3_pts(30)
  prof <- resolution_profile(pts, min_cell_n = 1, n_levels = 4, nstart = 2)
  b <- attr(prof, "bounds")
  expect_identical(b$ceiling, 29L)
  expect_output(print(prof), "ceiling 29, one short of the 30 points")
  expect_match(spatialkit:::.ladder_edge(29L, prof$levels, prof$levels, b),
               "one short of the points")
  # With repeat visits the distinct locations do bind, and say so.
  rep5 <- resolution_profile(r3_visits(r3_stations, 30), min_cell_n = 1, n_levels = 4,
                             nstart = 2)
  expect_output(print(rep5), "ceiling 5 from 5 distinct locations")
})
