# tests/testthat/test-review2-resolution.R
# ---------------------------------------------------------------------------
# Regressions from the second review of the resolution step:
# resolution_profile() and determine_optimal_levels().
# ---------------------------------------------------------------------------

# A hand-made sac_range, so the variogram-based columns exist without gstat.
r2_sac <- function(range, nugget = 1, psill = 1, crs = sf::st_crs(32632), ...) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = crs, ...)
}

r2_pts <- function(n = 300, seed = 1, ext = 1000, x0 = 5e5, y0 = 5e6, crs = 32632) {
  set.seed(seed)
  d <- data.frame(x = x0 + runif(n, 0, ext), y = y0 + runif(n, 0, ext), w = rnorm(n))
  d$z <- d$w + rnorm(n)
  sf::st_as_sf(d, coords = c("x", "y"), crs = crs)
}

# Every R warning an expression raises, muffled.
r2_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}


test_that("one empty point does not turn determine_optimal_levels() into 1", {
  set.seed(1)
  pts <- sf::st_as_sf(
    data.frame(x = 5e5 + c(runif(25, 0, 10), runif(25, 90, 100)),
               y = 5e6 + c(runif(25, 0, 10), runif(25, 90, 100))),
    coords = c("x", "y"), crs = 32632)
  clean <- determine_optimal_levels(pts, max_levels = 6)
  holed <- rbind(pts[1:10, ], sf::st_sf(geometry = sf::st_sfc(sf::st_point(), crs = 32632)),
                 pts[11:50, ])
  out <- r2_warnings(determine_optimal_levels(holed, max_levels = 6))
  expect_identical(out$value, clean)
  expect_true(any(grepl("dropping 1 point", out$warnings)))
  # Under a split the positions index the layer as passed: the empty row is
  # in neither half, and every other row is in one.
  set.seed(2)
  d <- data.frame(x = runif(80, 0, 1000), y = runif(80, 0, 1000), w = rnorm(80))
  d$z <- d$w + rnorm(80)
  p80 <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  p81 <- rbind(p80[1:40, ], sf::st_sf(w = 0, z = 0,
                                      geometry = sf::st_sfc(sf::st_point(), crs = 32632)),
               p80[41:80, ])
  sp <- attr(suppressWarnings(determine_optimal_levels(p81, max_levels = 6,
                                                       select_on = "split")), "split")
  expect_setequal(c(sp$selection, sp$estimation), setdiff(seq_len(81), 41L))
})


test_that("a few missing predictor values do not take Moran's z away from determine_optimal_levels()", {
  # Three NAs in 400 rows used to make every cell mean holding one NA and the
  # whole window fall back to the geometric ranking.
  pts <- r2_pts(400, seed = 2)
  holed <- pts
  holed$w[c(5, 50, 150)] <- NA
  dol <- function(p) suppressWarnings(determine_optimal_levels(
    p, max_levels = 40, response_var = "z", predictor_vars = "w", criterion = "morans_i"))
  a <- dol(pts); b <- dol(holed)
  da <- attr(a, "diagnostics"); db <- attr(b, "diagnostics")
  expect_false(is.null(db))
  ks <- da$eval_ks[is.finite(da$moran_z[da$eval_ks])]
  expect_gt(length(ks), 0L)
  expect_true(all(is.finite(db$moran_z[ks])))
  expect_equal(db$moran_z[ks], da$moran_z[ks], tolerance = 0.2)
  lines <- capture_spatialkit_log(dol(holed))
  expect_true(log_has(lines, "3 of 400 row"))
})


test_that("a misspelt column is an error in determine_optimal_levels(), and a missing one a warning", {
  pts <- r2_pts(200)
  expect_error(determine_optimal_levels(pts, response_var = "Z", predictor_vars = "w"),
               "column 'Z' not found")
  expect_error(determine_optimal_levels(pts, response_var = "z", predictor_vars = c("w", "Elev")),
               "Elev.*not found")
  out <- r2_warnings(determine_optimal_levels(pts, max_levels = 6, response_var = "z",
                                              criterion = "morans_i"))
  expect_true(any(grepl("predictor_vars were not given", out$warnings)))
  expect_null(attr(out$value, "diagnostics"))
})


test_that("repeat visits let the ladder reach one cell per location", {
  # Five stations visited thirty times each: k-means can put a cell on every
  # station (WSS 0), which the cap one short of the locations used to drop.
  st <- data.frame(x = 5e5 + c(0, 1000, 0, 1000, 500), y = 5e6 + c(0, 0, 1000, 1000, 500))
  d  <- st[rep(1:5, each = 30), ]
  set.seed(3)
  d$z <- rnorm(150) + rep(1:5, each = 30)
  pp <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  prof <- resolution_profile(pp, "z", levels = 2:5, nstart = 3, sac = r2_sac(300))
  expect_identical(prof$levels, 2:5)
  expect_identical(prof$wss[4], 0)
  expect_true(is.finite(prof$cp[4]))
  b <- attr(resolution_profile(pp, min_cell_n = 1, n_levels = 4, nstart = 3), "bounds")
  expect_identical(b$ceiling, 5L)
  expect_identical(b$ceiling_from, "distinct locations")
  # Without repeats the cap stays one short of the points, which
  # stats::kmeans() needs.
  expect_identical(resolution_profile(r2_pts(20), levels = 2:20, nstart = 2)$levels, 2:19)
  # determine_optimal_levels() reaches k = 5 as well.
  lines <- capture_spatialkit_log(suppressWarnings(
    determine_optimal_levels(pp, max_levels = 12)))
  expect_true(log_has(lines, "k = 1 to 5"))
})


test_that("reliability does not move when the layer is rotated", {
  # A 3000 x 120 strip: its bounding box is 21 times the hull's area once it
  # is turned 45 degrees, and reliability's domain term used to be read off
  # the box (argmax 8 axis-aligned, 2 rotated).
  strip <- function(theta) {
    set.seed(1)
    u <- cbind(runif(600, -1500, 1500), runif(600, -60, 60))
    R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2)
    xy <- u %*% t(R)
    sf::st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2], z = rnorm(600)),
                 coords = c("x", "y"), crs = 32632)
  }
  sac <- r2_sac(90, nugget = 2, psill = 1)
  lv  <- c(2, 4, 8, 16, 32)
  a <- resolution_profile(strip(0), "z", sac = sac, levels = lv, nstart = 2)
  b <- resolution_profile(strip(pi / 4), "z", sac = sac, levels = lv, nstart = 2)
  expect_equal(b$reliability, a$reliability, tolerance = 0.02)
  expect_identical(select_resolution(a, "reliability")$best,
                   select_resolution(b, "reliability")$best)
})


test_that("a sac is read in its own CRS", {
  skip_if_not(!is.na(sf::st_crs(2263)))
  # The same range in metres (UTM 18N) and in US feet (NY Long Island state
  # plane) on the same points: the feet version put the floor at 2 where
  # the metre version put it at 8.
  pts <- r2_pts(300, seed = 4, ext = 5000, x0 = 583000, y0 = 4507000, crs = 32618)
  m  <- resolution_profile(pts, "z", sac = r2_sac(1800, crs = sf::st_crs(32618)),
                           n_levels = 4, nstart = 2)
  ft <- resolution_profile(pts, "z", sac = r2_sac(1800 / 0.3048006, crs = sf::st_crs(2263)),
                           n_levels = 4, nstart = 2)
  expect_identical(attr(ft, "bounds")$floor, attr(m, "bounds")$floor)
  expect_gt(attr(m, "bounds")$floor, 2L)
  expect_equal(ft$reliability, m$reliability, tolerance = 0.01)
})


test_that("a rejected range gives cp its nugget and reliability nothing, and says so", {
  pts <- r2_pts(300)
  rej <- r2_sac(NA_real_, nugget = 0.5)
  attr(rej, "variogram_model")$range[2] <- 3000
  attr(rej, "rejected_reason") <- "fitted range exceeds the largest lag fitted"
  out <- r2_warnings(resolution_profile(pts, "z", sac = rej, n_levels = 4, nstart = 2))
  prof <- out$value
  expect_true(any(grepl("fitted range exceeds the largest lag fitted", out$warnings)))
  expect_true(all(is.finite(prof$cp)))
  expect_true(all(is.na(prof$reliability)))
  expect_equal(attr(prof, "variogram")$nugget, 0.5)
  expect_error(select_resolution(prof, "reliability"), "identified range")
  # A model that did not converge gives neither.
  attr(rej, "rejected_reason") <- "variogram model did not converge"
  out <- r2_warnings(resolution_profile(pts, "z", sac = rej, n_levels = 4, nstart = 2))
  expect_true(any(grepl("did not converge", out$warnings)))
  expect_true(all(is.na(out$value$cp)) && all(is.na(out$value$reliability)))
})


test_that("a zero nugget, and a sac of the other variable, are warned about", {
  pts <- r2_pts(300)
  out <- r2_warnings(resolution_profile(pts, "z", sac = r2_sac(300, nugget = 0),
                                        n_levels = 4, nstart = 2))
  expect_true(any(grepl("nugget is 0", out$warnings)))
  # A residual variogram applied to the raw response, and the reverse.
  det <- r2_sac(300, detrended = TRUE)
  out <- r2_warnings(resolution_profile(pts, "z", sac = det, n_levels = 4, nstart = 2))
  expect_true(any(grepl("variogram of detrended residuals", out$warnings)))
  expect_true(isTRUE(attr(out$value, "variogram")$detrended))
  raw <- r2_sac(300, detrended = FALSE)
  out <- r2_warnings(resolution_profile(pts, "z", "w", sac = raw, n_levels = 4, nstart = 2))
  expect_true(any(grepl("variogram of the raw response", out$warnings)))
  # Matching, or unlabelled, raises nothing.
  expect_length(r2_warnings(resolution_profile(pts, "z", "w", sac = det, n_levels = 4,
                                               nstart = 2))$warnings, 0L)
  expect_length(r2_warnings(resolution_profile(pts, "z", sac = r2_sac(300), n_levels = 4,
                                               nstart = 2))$warnings, 0L)
})


test_that("a sac supplied under select_on = 'split' is flagged", {
  pts <- r2_pts(200)
  lines <- capture_spatialkit_log(
    resolution_profile(pts, "z", sac = r2_sac(300), n_levels = 4, nstart = 2,
                       select_on = "split"))
  expect_true(log_has(lines, "must be fitted on the selection half"))
})


test_that("range_floor = FALSE starts the ladder at 2 and still reports the floor", {
  pts <- r2_pts(300)
  sac <- r2_sac(250)
  on  <- resolution_profile(pts, "z", sac = sac, n_levels = 5, nstart = 2)
  off <- resolution_profile(pts, "z", sac = sac, n_levels = 5, nstart = 2, range_floor = FALSE)
  expect_gt(attr(on, "bounds")$floor, 2L)
  expect_identical(min(on$levels), attr(on, "bounds")$floor)
  expect_identical(min(off$levels), 2L)
  expect_identical(attr(off, "bounds")$floor, attr(on, "bounds")$floor)
  expect_false(attr(off, "bounds")$range_floor)
  expect_output(print(off), "not applied")
  expect_error(resolution_profile(pts, range_floor = NA), "TRUE or FALSE")
})
