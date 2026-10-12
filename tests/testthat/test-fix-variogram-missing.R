# tests/testthat/test-fix-variogram-missing.R
# ---------------------------------------------------------------------------
# Rows with a missing response are dropped before the n_max subsample is
# drawn, and the subsample does not depend on the order of the rows.
#
# The subsample used to be drawn over every row, so a sparse response in a
# large layer was estimated from about n_max / n of its usable rows: 300
# valued rows among 20,300 gave 239.8 from the 74 that fell in the draw
# (124.0 from the 300 alone), a different number for each seed (87.4 to
# 242.3), and 100 valued rows among 20,100 came back NA as "20 point(s) with
# a finite value to model".
# ---------------------------------------------------------------------------

# An exponential field (effective range R, nugget share nug) on a 1000 m
# square, with x and y kept as columns so layers can be stacked.
fvm_field <- function(n, seed, R = 300, nug = 0.2) {
  set.seed(seed)
  xy <- cbind(5e5 + runif(n, 0, 1000), 5e6 + runif(n, 0, 1000))
  z <- as.numeric(t(chol((1 - nug) * exp(-3 * as.matrix(stats::dist(xy)) / R) +
                           diag(nug, n))) %*% rnorm(n))
  sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z), coords = c("x", "y"),
               crs = 32632, remove = FALSE)
}
# m rows with no response, scattered over the same square.
fvm_pad <- function(m, seed) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = 5e5 + runif(m, 0, 1000), y = 5e6 + runif(m, 0, 1000),
                          z = NA_real_), coords = c("x", "y"), crs = 32632,
               remove = FALSE)
}

test_that("rows with no response do not use up the n_max subsample", {
  skip_if_not_installed("gstat")
  valued <- fvm_field(300, seed = 3)
  padded <- rbind(valued, fvm_pad(20000, seed = 4))
  alone <- estimate_sac_range(valued, "z")
  both  <- estimate_sac_range(padded, "z")
  expect_true(is.finite(alone))
  # The same estimate, to the last bit: the 20,000 empty rows change nothing.
  expect_identical(as.numeric(both), as.numeric(alone))
  expect_identical(attr(both, "variogram"), attr(alone, "variogram"))
  expect_identical(attr(both, "variogram_model"), attr(alone, "variogram_model"))
  # The counts say what happened to the rows.
  expect_identical(attr(alone, "n_used"), 300L)
  expect_identical(attr(alone, "n_missing"), 0L)
  expect_identical(attr(alone, "n_subsampled"), 0L)
  expect_identical(attr(both, "n_used"), 300L)
  expect_identical(attr(both, "n_missing"), 20000L)
  expect_identical(attr(both, "n_subsampled"), 0L)
  # No subsample is drawn, so the seed cannot matter.
  for (s in c(1L, 2L, 99L))
    expect_identical(as.numeric(estimate_sac_range(padded, "z", seed = s)),
                     as.numeric(alone))
  # The rows left out are logged with their count.
  lines <- capture_spatialkit_log(estimate_sac_range(padded, "z"))
  expect_true(log_has(lines, "20000 of 20300 row\\(s\\) have a missing or non-finite response"))
})

test_that("a sparse response is not refused for the rows the subsample missed", {
  skip_if_not_installed("gstat")
  valued <- fvm_field(100, seed = 5)
  padded <- rbind(valued, fvm_pad(20000, seed = 6))
  alone <- suppressWarnings(estimate_sac_range(valued, "z"))
  both  <- suppressWarnings(estimate_sac_range(padded, "z"))
  # It used to be NA with "20 point(s) with a finite value to model, fewer
  # than the 30 a variogram range is estimated from": 100 exist.
  expect_false(grepl("fewer than the 30", attr(both, "rejected_reason") %||% ""))
  expect_identical(as.numeric(both), as.numeric(alone))
  expect_identical(attr(both, "rejected_reason"), attr(alone, "rejected_reason"))
  expect_identical(attr(both, "n_used"), 100L)
  # When too few rows carry a value, the count in the reason is the layer's.
  few <- rbind(valued[1:12, ], fvm_pad(500, seed = 7))
  r <- suppressWarnings(estimate_sac_range(few, "z"))
  expect_true(is.na(r))
  expect_identical(attr(r, "rejected_reason"), paste0(
    "12 point(s) with a finite value to model, fewer than the 30 a variogram ",
    "range is estimated from"))
})

test_that("with more valued rows than n_max, the subsample is drawn from them alone", {
  skip_if_not_installed("gstat")
  valued <- fvm_field(400, seed = 8)
  padded <- rbind(valued, fvm_pad(1500, seed = 9))
  alone <- suppressWarnings(estimate_sac_range(valued, "z", n_max = 150))
  both  <- suppressWarnings(estimate_sac_range(padded, "z", n_max = 150))
  expect_identical(as.numeric(both), as.numeric(alone))
  expect_identical(attr(both, "variogram"), attr(alone, "variogram"))
  expect_identical(attr(both, "n_used"), 150L)
  expect_identical(attr(both, "n_subsampled"), 250L)
  expect_identical(attr(both, "n_missing"), 1500L)
  # The 150 points are all of the variogram's: every pair is one of theirs.
  expect_lte(sum(attr(both, "variogram")$np), choose(150, 2))
})

test_that("rows with no response take no part in choosing the projection", {
  skip_if_not_installed("gstat")
  # Lon/lat input is projected to a CRS picked from the layer's extent.  The
  # empty rows reach far beyond the valued ones here, and used to move that
  # choice (a UTM zone for the valued rows alone, an equal-area CRS padded).
  valued <- sf::st_transform(fvm_field(120, seed = 10), 4326)
  set.seed(11)
  far <- sf::st_as_sf(data.frame(x = runif(300, -30, 40), y = runif(300, 20, 65),
                                 z = NA_real_), coords = c("x", "y"), crs = 4326,
                      remove = FALSE)
  valued$x <- sf::st_coordinates(valued)[, 1]; valued$y <- sf::st_coordinates(valued)[, 2]
  padded <- rbind(valued, far)
  alone <- suppressWarnings(estimate_sac_range(valued, "z"))
  both  <- suppressWarnings(estimate_sac_range(padded, "z"))
  expect_true(attr(both, "crs") == attr(alone, "crs"))
  expect_identical(as.numeric(both), as.numeric(alone))
  expect_identical(attr(both, "max_dist"), attr(alone, "max_dist"))
})

test_that("a missing predictor leaves its row out before the subsample too", {
  skip_if_not_installed("gstat")
  pts <- fvm_field(200, seed = 12)
  set.seed(13); pts$w <- rnorm(nrow(pts))
  gaps <- pts
  gaps$w[c(3, 50, 120)] <- c(NA, Inf, NA)
  full <- suppressWarnings(estimate_sac_range(pts[-c(3, 50, 120), ], "z", "w"))
  r <- suppressWarnings(estimate_sac_range(gaps, "z", "w"))
  expect_identical(as.numeric(r), as.numeric(full))
  expect_identical(attr(r, "n_used"), 197L)
  expect_identical(attr(r, "n_missing"), 3L)
  expect_true(attr(r, "detrended"))
})

test_that("a predictor missing on every row falls back to the raw response with a warning", {
  skip_if_not_installed("gstat")
  # Leaving out every row would leave nothing to model, and the trend cannot
  # be fitted either way, so the documented fallback stands: the raw response,
  # an R warning, and detrended = FALSE.
  pts <- fvm_field(150, seed = 14)
  raw <- suppressWarnings(estimate_sac_range(pts, "z"))
  pts$bad <- NA_real_
  expect_warning(fb <- estimate_sac_range(pts, "z", predictor_vars = "bad"),
                 "OLS detrending on bad failed .* fitted to the RAW response")
  expect_false(attr(fb, "detrended"))
  expect_identical(as.numeric(fb), as.numeric(raw))
  expect_identical(attr(fb, "n_used"), 150L)
})

test_that("a predictor present on fewer than 30 rows gives the NA with their count, as it did", {
  skip_if_not_installed("gstat")
  # 240 rows with a response, the predictor on 20 of them.  A trend can be
  # fitted to 20 rows and a variogram cannot, so this is the unclassed NA
  # that says how many rows there were to model.  It is not the raw-response
  # fallback: that is for a trend that cannot be fitted at all (the test
  # above), and a range of the raw response is not what was asked for.
  pts <- fvm_field(240, seed = 15)
  set.seed(16); pts$w <- rnorm(nrow(pts))
  pts$w[21:240] <- NA
  expect_no_warning(r <- estimate_sac_range(pts, "z", predictor_vars = "w"))
  expect_true(is.na(r))
  expect_false(inherits(r, "sac_range"))
  expect_identical(names(attributes(r)), "rejected_reason")
  expect_identical(attr(r, "rejected_reason"), paste0(
    "20 point(s) with a finite value to model, fewer than the 30 a variogram ",
    "range is estimated from"))
  # One complete row is still a trend to fit and nothing to fit a variogram to.
  pts$w[2:240] <- NA
  expect_no_warning(r1 <- estimate_sac_range(pts, "z", predictor_vars = "w", detrend = "reml"))
  expect_identical(attr(r1, "rejected_reason"), paste0(
    "1 point(s) with a finite value to model, fewer than the 30 a variogram ",
    "range is estimated from"))
  # Thirty complete rows are enough, and the result is of the residuals.
  pts$w <- NA_real_
  set.seed(17); pts$w[1:40] <- rnorm(40)
  r40 <- suppressWarnings(estimate_sac_range(pts, "z", predictor_vars = "w",
                                             structure_alpha = 1))
  expect_s3_class(r40, "sac_range")
  expect_true(attr(r40, "detrended"))
  expect_identical(attr(r40, "n_used"), 40L)
  expect_identical(attr(r40, "n_missing"), 200L)
})

test_that("above n_max the answer does not depend on the order of the rows", {
  skip_if_not_installed("gstat")
  # The subsample indexes rows, so the same 500 points in another order gave
  # another 250 of them: 9 of 10 long-range layers changed their answer.
  for (s in 1:3) {
    pts <- fvm_field(500, seed = 300 + s)
    set.seed(40 + s); shuffled <- pts[sample(nrow(pts)), ]
    a <- suppressWarnings(estimate_sac_range(pts, "z", n_max = 250))
    b <- suppressWarnings(estimate_sac_range(shuffled, "z", n_max = 250))
    expect_identical(b, a, info = paste("layer", s))
    expect_identical(attr(a, "n_subsampled"), 250L)
  }
})
