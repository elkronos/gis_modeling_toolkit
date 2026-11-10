# cv_block_size_sweep()'s default ladder gives the same blocks in every unit.
#
# The ladder's bottom rung is a twenty-fifth of the shorter side, so that side
# holds exactly 25 blocks in exact arithmetic.  The rung was computed as
# exp(log(side / 25)), a hair above or below side / 25 by the rounding of the
# unit, and the blocks per side were floor(side / rung): 24 in one unit and 25
# in another.  On an 8 x 8 lattice at 250 m the bottom rung's blocks were
# 72.9 m wide in metres and 70 m in units of 1e-6 m.  The ladder's ends are
# now side / 25 and the top exactly, and a quotient within .TIE_TOL below a
# whole number is that number.

# An 8 x 8 lattice at 250 m and 60 points in a 4000 x 3000 m box, both in
# UTM zone 17N.
.sl_lattice <- function() {
  g <- expand.grid(i = 0:7, j = 0:7)
  list(xy = cbind(600000 + 250 * g$i, 3990000 + 250 * g$j),
       att = data.frame(z = 2 + sin(0.9 * g$i) + cos(0.7 * g$j),
                        a = ((3 * g$i + 5 * g$j) %% 7) / 7 - 0.5))
}
.sl_points <- function(n = 60L) {
  withr::with_seed(20261106L, {
    px <- round(8 * (600000 + stats::runif(n, 0, 4000))) / 8
    py <- round(8 * (3990000 + stats::runif(n, 0, 3000))) / 8
    a  <- stats::runif(n, -1, 1)
    z  <- 2 + sin((px - 600000) / 800) + cos((py - 3990000) / 700) + 0.5 * a +
      stats::rnorm(n, 0, 0.2)
  })
  list(xy = cbind(px, py), att = data.frame(z = z, a = a))
}

# The same layer in metres, in units of 1e-6 m and in US survey feet: the
# same projection, only the unit changes.
.sl_units <- list(m = 1, um = 1e-6, us_ft = 1200 / 3937)
.sl_crs <- function(u) switch(u,
  m = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs",
  um = "+proj=utm +zone=17 +datum=WGS84 +to_meter=1e-06 +no_defs",
  us_ft = "+proj=utm +zone=17 +datum=WGS84 +units=us-ft +no_defs")
.sl_sf <- function(L, u) {
  f <- .sl_units[[u]]
  sf::st_as_sf(cbind(L$att, X = L$xy[, 1] / f, Y = L$xy[, 2] / f),
               coords = c("X", "Y"), crs = .sl_crs(u))
}
.sl_fit <- function(train_sf) lm_spatial_fit(train_sf, "z", "a")
.sl_sweep <- function(d) suppressWarnings(cv_block_size_sweep(
  d, "z", "a", fit_fn = .sl_fit, k = 3, n_sizes = 3, include_random = FALSE,
  sac = NA, quiet = TRUE, seed = 1))

test_that("the default ladder gives the same blocks in metres, 1e-6 m and US survey feet", {
  for (nm in c("lattice", "points")) {
    L <- if (nm == "lattice") .sl_lattice() else .sl_points()
    ref <- .sl_sweep(.sl_sf(L, "m"))
    for (u in c("um", "us_ft")) {
      f <- .sl_units[[u]]
      got <- .sl_sweep(.sl_sf(L, u))
      lbl <- sprintf("%s in %s", nm, u)
      expect_identical(got$blocks_used, ref$blocks_used, info = lbl)
      # Lengths scaled back to metres agree to the rounding of the unit.
      for (col in c("block_size", "block_width", "block_height"))
        expect_equal(got[[col]] * f, ref[[col]], tolerance = 1e-12,
                     info = paste(lbl, col))
    }
  }
})

test_that("the bottom rung holds 25 blocks along the lattice's side in every unit", {
  L <- .sl_lattice()
  for (u in names(.sl_units)) {
    f <- .sl_units[[u]]
    sw <- .sl_sweep(.sl_sf(L, u))
    # 1750 m / 25 = 70 m.
    expect_equal(sw$block_width[1L] * f, 70, tolerance = 1e-12, info = u)
    expect_equal(sw$block_height[1L] * f, 70, tolerance = 1e-12, info = u)
  }
})

test_that("the ladder's ends are side / 25 and the top, to the bit", {
  for (nm in c("lattice", "points")) {
    L <- if (nm == "lattice") .sl_lattice() else .sl_points()
    for (u in names(.sl_units)) {
      d  <- .sl_sf(L, u)
      bb <- sf::st_bbox(d)
      side <- min(as.numeric(bb["xmax"] - bb["xmin"]), as.numeric(bb["ymax"] - bb["ymin"]))
      top  <- .sweep_largest_size_for_k(bb, 3, cap = side / 2 * (1 - 1e-9))
      sw <- .sl_sweep(d)
      lbl <- sprintf("%s in %s", nm, u)
      expect_identical(sw$block_size[1L], side / 25, info = lbl)
      expect_identical(sw$block_size[nrow(sw)], top, info = lbl)
    }
  }
})

test_that("a size within .TIE_TOL above side / i gives i blocks along that side", {
  bb <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1750, ymax = 1750))
  # 1750 / 25 rounded up by a few units in the last place.
  b <- 70 * (1 + 4 * .Machine$double.eps)
  expect_lt(1750 / b, 25)
  d <- .block_dims_from_size(bb, b)
  expect_identical(c(d$nx, d$ny), c(25, 25))
  # A size more than .TIE_TOL above it gives one block fewer.
  d <- .block_dims_from_size(bb, 70 * (1 + 1e-6))
  expect_identical(c(d$nx, d$ny), c(24, 24))
})

test_that("a size up to .TIE_TOL above side / i gives i blocks, and not beyond", {
  bb <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1750, ymax = 1750))
  # 1750 / b is 25 (1 - 9e-8) and 25 (1 - 1.1e-7): the first is within
  # .TIE_TOL below 25, the second is not.
  d <- .block_dims_from_size(bb, 70 * (1 + 9e-8))
  expect_identical(c(d$nx, d$ny), c(25, 25))
  d <- .block_dims_from_size(bb, 70 * (1 + 1.1e-7))
  expect_identical(c(d$nx, d$ny), c(24, 24))
})
