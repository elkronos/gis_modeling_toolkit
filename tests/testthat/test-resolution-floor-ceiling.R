# tests/testthat/test-resolution-floor-ceiling.R
# ---------------------------------------------------------------------------
# A range floor above the ceiling of the classic ladder.
#
# On a layer of about 900 rows the floor read 2,506 cells on a ladder whose
# top was 99, and what the caller could read of that was uneven: the floor
# above the data's own ceiling raised spatialkit_not_supported, but the
# floor above the ceiling that sample_n sets was a log line alone.
# ---------------------------------------------------------------------------

# Points on a 100 km square in EPSG:5070 with a response that varies over
# `scale` metres.
.fc_layer <- function(n, scale, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, 0, 1e5), y = stats::runif(n, 0, 1e5))
  d$z <- sin(d$x / scale) + cos(d$y / scale) + stats::rnorm(n, sd = 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 5070)
}

# A sac_range of the given range, made by hand so that the floor is known:
# ceiling(1e10 / range^2) cells on the 100 km square.
.fc_sac <- function(pts, range) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.2, 1), range = c(0, range / 3),
                   stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(pts))
}

.fc_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws)
}

.fc_of <- function(r, class) r$warnings[vapply(r$warnings, inherits, logical(1), class)]

# 3000 points, fitted on a subsample of 100: the subsample fits at most 50
# cells, the layer supports 333 at min_cell_n = 9, and a range of 10 km asks
# for 100.
.fc_subsample <- function(range_floor = TRUE, ...) {
  pts <- .fc_layer(3000, 8000, 2)
  .fc_catch(resolution_profile(pts, "z", sac = .fc_sac(pts, 10000), n_levels = 8,
                               sample_n = 100, nstart = 2, range_floor = range_floor, ...))
}

# 300 points: the layer supports 33 cells, and a range of 2 km asks for 2500.
.fc_data <- function(range_floor = TRUE) {
  pts <- .fc_layer(300, 3000, 1)
  .fc_catch(resolution_profile(pts, "z", sac = .fc_sac(pts, 2000), n_levels = 8,
                               nstart = 2, range_floor = range_floor))
}


test_that("a floor above the ceiling sample_n sets is a classed warning, not a log line alone", {
  skip_if_not_installed("gstat")
  lines <- capture_spatialkit_log(r <- .fc_subsample(), level = logger::WARN)
  b <- attr(r$value, "bounds")
  expect_identical(b$ceiling_from, "sample_n")
  expect_identical(b$ceiling, 50L)
  expect_identical(b$floor, 100L)
  expect_true(b$supported)
  w <- .fc_of(r, "spatialkit_floor_above_ceiling")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$reason, "sample_n")
  expect_identical(w[[1L]]$floor, 100)
  expect_identical(w[[1L]]$ceiling, 50L)
  expect_identical(w[[1L]]$n, 3000L)
  expect_identical(w[[1L]]$n_sample, 100L)
  expect_true(startsWith(conditionMessage(w[[1L]]), paste0(
    "resolution_profile(): cells no wider than the autocorrelation range (10000) need at ",
    "least 100 of them, which the 3000 points support, but k-means on the 100-point ",
    "subsample can fit at most 50; raise `sample_n`. The profile runs from 2 to 50.")))
  # The layer supports the floor, so it is not spatialkit_not_supported.
  expect_length(.fc_of(r, "spatialkit_not_supported"), 0L)
  # Still logged, once.
  expect_identical(sum(grepl("raise `sample_n`", lines, fixed = TRUE)), 1L)
  # tryCatch() sees it, whatever range_floor is.
  pts <- .fc_layer(3000, 8000, 2)
  expect_identical(tryCatch(resolution_profile(pts, "z", sac = .fc_sac(pts, 10000),
                                               n_levels = 8, sample_n = 100, nstart = 2),
                            spatialkit_floor_above_ceiling = function(w) "caught"),
                   "caught")
  # A floor above the data's own ceiling is still spatialkit_not_supported,
  # and only that.
  d <- .fc_data()
  expect_length(.fc_of(d, "spatialkit_not_supported"), 1L)
  expect_length(.fc_of(d, "spatialkit_floor_above_ceiling"), 0L)
  expect_false(attr(d$value, "bounds")$supported)
})
