# tests/testthat/test-fix-resolution-range-floor.R
# ---------------------------------------------------------------------------
# The range floor is reported and not applied unless asked: the ladder of
# resolution_profile() starts at 2 by default, so an estimated (and on a
# response with no structure, spurious) range no longer sets the fewest cells
# a criterion can choose.  `reliability`, whose own optimum lies below the
# floor, is read at or above it by select_resolution().
# ---------------------------------------------------------------------------

frf_field <- function(n = 400, seed = 2, a = 100, nugget = 0.3) {
  set.seed(seed)
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / a) + diag(1e-8, n))) %*% rnorm(n)) +
    rnorm(n, sd = sqrt(nugget))
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
}
frf_sac <- function(pts, range = 240, nugget = 0.3, psill = 1) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(pts))
}
frf_quiet <- function(expr) suppressWarnings(suppressMessages(expr))


test_that("range_floor = FALSE is the default: the ladder starts at 2 and the floor is reported", {
  expect_identical(formals(resolution_profile)$range_floor, FALSE)
  pts <- frf_field()
  sac <- frf_sac(pts)
  p <- frf_quiet(resolution_profile(pts, "z", sac = sac, n_levels = 8, nstart = 3))
  b <- attr(p, "bounds")
  # A range of 240 on about a square kilometre: a floor of 16 or so, which
  # used to be the first level of the ladder.
  expect_gt(b$floor, 10L)
  expect_identical(b$floor, as.integer(ceiling(b$area / 240^2)))
  expect_true(b$supported)
  expect_false(b$range_floor)
  expect_identical(min(p$levels), 2L)
  expect_identical(max(p$levels), b$ceiling)
  expect_true(any(p$levels < b$floor))
  expect_output(print(p), sprintf("floor %d from range 240, not applied", b$floor))

  # Asked for, the floor starts the ladder as the default used to.
  on <- frf_quiet(resolution_profile(pts, "z", sac = sac, n_levels = 8, nstart = 3,
                                     range_floor = TRUE))
  expect_identical(min(on$levels), attr(on, "bounds")$floor)
  expect_true(attr(on, "bounds")$range_floor)
  # Cp over the ladder from 2 is free to choose below the floor; its edge
  # label never names a floor that was not applied.
  s <- select_resolution(p, "cp")
  expect_false(isTRUE(grepl("range floor", s$edge)))
})

test_that("a spurious range on white noise does not start the ladder above 2", {
  skip_if_not_installed("gstat")
  # White noise on a clustered layout, two of the seeds on which the profile's
  # own variogram returned ranges of 40 and 45 km and the ladder started at 98
  # and at 75 cells of 100.  Whatever range the variogram returns here, the
  # default ladder starts at 2.
  # The generator of the layers this was measured on, random stream and all
  # (the draws that are not used keep the stream in step with it).
  sim_xy <- function(n, towns) {
    k <- sample(nrow(towns), n, replace = TRUE)
    in_town <- runif(n) < 0.6
    x <- ifelse(in_town, towns[k, 1] + rnorm(n, 0, 12e3), runif(n, 0, 500e3))
    y <- ifelse(in_town, towns[k, 2] + rnorm(n, 0, 12e3), runif(n, 0, 320e3))
    cbind(x = pmin(pmax(x, 0), 500e3), y = pmin(pmax(y, 0), 320e3))
  }
  noise_layer <- function(seed, n = 900) {
    set.seed(seed)
    towns <- cbind(runif(8, 40e3, 460e3), runif(8, 30e3, 290e3))
    runif(40); runif(40); rnorm(60); rnorm(960); runif(40)
    xy <- sim_xy(n, towns)
    sim_xy(2000, towns)
    sf::st_as_sf(data.frame(x = 3e5 + xy[, 1], y = 4.5e6 + xy[, 2], z = rnorm(n)),
                 coords = c("x", "y"), crs = 32615)
  }
  for (seed in c(1010, 1018, 1011)) {
    p <- frf_quiet(resolution_profile(noise_layer(seed), "z", n_levels = 6, nstart = 2))
    expect_identical(min(p$levels), 2L)
    expect_false(attr(p, "bounds")$range_floor)
  }
})

test_that("select_resolution() reads reliability at or above the range floor", {
  pts <- frf_field()
  sac <- frf_sac(pts)
  p <- frf_quiet(resolution_profile(pts, "z", sac = sac, n_levels = 10, nstart = 3))
  b <- attr(p, "bounds")
  rel <- select_resolution(p, "reliability")
  # The curve's own maximum is below the floor; read from the floor up, the
  # pick is the floor, which is a level of the ladder, and the result says
  # the floor chose.
  free <- select_resolution(p, "reliability", range_floor = FALSE)
  expect_identical(free$best, p$levels[which.max(p$reliability)])
  expect_lt(free$best, b$floor)
  expect_true(is.na(free$range_floor))
  expect_identical(rel$best, b$floor)
  expect_identical(rel$edge, "the range floor (area / range^2)")
  expect_identical(rel$range_floor, b$floor)
  expect_true(rel$at_floor)
  expect_match(rel$edge, "range floor")
  expect_true(all(rel$flat >= b$floor))
  expect_output(print(rel), "range floor")
  # The other criteria are not held to it.
  expect_true(is.na(select_resolution(p, "cp")$range_floor))
  expect_error(select_resolution(p, "reliability", range_floor = NA), "TRUE or FALSE")

  # A ladder the caller set is read as given, as it always was.
  own <- frf_quiet(resolution_profile(pts, "z", sac = sac, levels = c(3, 6, 12, 24, 40),
                                      nstart = 3))
  expect_true(attr(own, "bounds")$levels_given)
  expect_identical(select_resolution(own, "reliability")$best,
                   own$levels[which.max(own$reliability)])
  # A ladder that starts at the floor gives the floor itself, as before.
  on <- frf_quiet(resolution_profile(pts, "z", sac = sac, n_levels = 10, nstart = 3,
                                     range_floor = TRUE))
  ron <- select_resolution(on, "reliability")
  expect_identical(ron$best, on$levels[which.max(on$reliability)])
  if (ron$best == attr(on, "bounds")$floor)
    expect_identical(ron$edge, "the range floor (area / range^2)")
})

test_that("the floor is a level of the ladder, so reliability returns the floor itself", {
  # With no level on the floor the pick was the next level up: 54 cells for
  # a floor of 45 at the default 20 levels, and a step of the ladder is a
  # factor of 2 to 4 on a short one.
  pts <- frf_field()
  sac <- frf_sac(pts)
  for (nl in c(5L, 8L, 20L)) {
    p <- frf_quiet(resolution_profile(pts, "z", sac = sac, n_levels = nl, nstart = 2))
    b <- attr(p, "bounds")
    expect_gt(b$floor, 2L)
    expect_lt(b$floor, b$ceiling)
    expect_true(b$floor %in% p$levels)
    expect_identical(range(p$levels), c(2L, b$ceiling))
    expect_lte(nrow(p), nl + 1L)
    rel <- select_resolution(p, "reliability")
    expect_identical(rel$best, b$floor)
    expect_identical(rel$edge, "the range floor (area / range^2)")
    expect_true(rel$at_floor)
  }
  # The same ladder without a range has no such level: the floor adds one.
  geo <- frf_quiet(resolution_profile(pts, n_levels = 5L, nstart = 2))
  p5 <- frf_quiet(resolution_profile(pts, "z", sac = sac, n_levels = 5L, nstart = 2))
  expect_identical(setdiff(p5$levels, geo$levels), attr(p5, "bounds")$floor)

  # No level is added where there is none to add: a floor of 2 is the first
  # level already, a floor the layer cannot support is not a level of the
  # ladder, and a ladder the caller set is scored as given.
  low <- frf_quiet(resolution_profile(pts, "z", sac = frf_sac(pts, range = 2000),
                                      n_levels = 5L, nstart = 2))
  expect_identical(attr(low, "bounds")$floor, 2L)
  expect_identical(low$levels, geo$levels)
  high <- frf_quiet(resolution_profile(pts, "z", sac = frf_sac(pts, range = 100),
                                       n_levels = 5L, nstart = 2))
  expect_false(attr(high, "bounds")$supported)
  expect_identical(high$levels, geo$levels)
  own <- frf_quiet(resolution_profile(pts, "z", sac = sac, levels = c(3, 6, 12, 24, 40),
                                      nstart = 2))
  expect_identical(own$levels, c(3L, 6L, 12L, 24L, 40L))
})
