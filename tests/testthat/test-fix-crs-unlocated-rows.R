# tests/testthat/test-fix-crs-unlocated-rows.R
# ---------------------------------------------------------------------------
# One row of a lon/lat layer at a coordinate that is not finite made
# .geod_distance() stop ("NAs are not allowed in subscripted assignments"),
# .crs_distance_error() turned the stop into NA for every candidate, and a
# continental layer was put in one UTM zone with no comparison and no log
# line.  Such a row now takes no part in the score, as an empty row does,
# and .geod_distance() gives NA for a pair that is not finite.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
stations <- function() {
  set.seed(1)
  data.frame(lon = stats::runif(39, -120, -75), lat = stats::runif(39, 30, 47))
}
as_points <- function(d) sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)

test_that(".geod_distance() gives NA for a pair that is not finite and the distances of the others", {
  fin <- .geod_distance(c(-120, -100), c(30, 40), c(-100, -80), c(40, 45))
  for (bad in list(Inf, -Inf, NA_real_, NaN)) {
    # Without a warning: tan() of a latitude that is not finite gives NaN
    # with "NaNs produced".
    expect_no_warning(d <- .geod_distance(c(-120, bad, -100), c(30, 40, 40), c(-100, -100, -80),
                                          c(40, 40, 45)))
    expect_length(d, 3L)
    expect_true(is.na(d[2L]), label = paste("lon1 =", bad))
    expect_identical(d[-2L], fin, label = paste("lon1 =", bad))
  }
  expect_no_warning(d <- .geod_distance(-104, c(40, Inf), -103, 41))
  expect_true(is.na(d[2L]))
  expect_identical(d[1L], .geod_distance(-104, 40, -103, 41))
  expect_true(is.na(.geod_distance(NA_real_, 40, NA_real_, 40)))
  expect_no_warning(d <- .geod_distance(-104, c(-Inf, NaN, NA_real_), c(-103, 0, 1), c(Inf, 41, NaN)))
  expect_true(all(is.na(d)))
})

test_that("a row at a coordinate that is not finite is scored as an empty row is", {
  # The 39 stations score 0.050690 in UTM zone 14 and get an Albers CRS
  # with "0.61% against the zone's 5.07%"; with a 40th row at (Inf, 40)
  # every candidate scored NA, and UTM zone 14 was kept with no log line
  # and an NA in `crs_choice`.
  st <- stations()
  e39 <- lapply(c(32614, 5070), function(code) .crs_distance_error(as_points(st), sf::st_crs(code)))
  for (bad in list(c(Inf, 40), c(-100, -Inf))) {
    x <- as_points(rbind(st, data.frame(lon = bad[1], lat = bad[2])))
    expect_identical(lapply(c(32614, 5070), function(code) .crs_distance_error(x, sf::st_crs(code))),
                     e39, label = paste(bad, collapse = ", "))
  }
  # Above 40 rows the sample runs along the points sorted by coordinate, and
  # the row at Inf, sorted last, took the last place in it.
  set.seed(2)
  st60 <- data.frame(lon = stats::runif(60, -120, -75), lat = stats::runif(60, 30, 47))
  e60 <- lapply(c(32614, 5070), function(code) .crs_distance_error(as_points(st60), sf::st_crs(code)))
  x <- as_points(rbind(st60, data.frame(lon = c(Inf, -100), lat = c(40, Inf))))
  expect_identical(lapply(c(32614, 5070), function(code) .crs_distance_error(x, sf::st_crs(code))), e60)
  # The chooser reads the layer as the 39 stations (BUG-137): it kept UTM
  # zone 14 when its extent read the row.  The zone's score is the stations'
  # and not NA.
  x <- quiet(ensure_projected(as_points(rbind(st, data.frame(lon = Inf, lat = 40)))))
  ch <- attr(x, "crs_choice")
  expect_identical(ch$distance_error[ch$name == "UTM zone 14"], e39[[1L]])
  ref <- quiet(ensure_projected(as_points(st)))
  expect_identical(sf::st_crs(x), sf::st_crs(ref))
  expect_identical(ch, attr(ref, "crs_choice"))
  expect_false(identical(sf::st_crs(x), sf::st_crs(32614)))
})
