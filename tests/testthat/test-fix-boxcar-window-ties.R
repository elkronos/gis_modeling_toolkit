# The boxcar kernel kept a neighbour when d <= h exactly.  On a regular grid
# many neighbours tie with the bw-th one; computed from coordinates in
# another unit (US survey feet, a millionth of a metre) their distances and
# h round apart by a unit in the last place, so which of them were in the
# window, and the fitted values, moved with the unit.  A squared distance
# within .TIE_TOL of h^2 is now at the edge, and inside.

skip_if_not_installed("sp")

# An 8 x 8 lattice at 250 m spacing, so many distances tie exactly in
# metres, with a response and two predictors.
.bx_lattice <- function() {
  g <- expand.grid(i = 0:7, j = 0:7)
  withr::with_seed(5, {
    a <- rnorm(64); b <- rnorm(64)
    z <- 2 + sin(g$i / 2) + cos(g$j / 3) + 0.6 * a - 0.4 * b +
      rnorm(64, sd = 0.2)
  })
  list(xy = cbind(600000 + 250 * g$i, 3990000 + 250 * g$j),
       att = data.frame(a = a, b = b, z = z))
}

.bx_layer <- function(L, xy, crs, perm = seq_len(nrow(xy))) {
  sf::st_as_sf(cbind(L$att, X = xy[, 1], Y = xy[, 2])[perm, ],
               coords = c("X", "Y"), crs = crs)
}

.bx_fit <- function(d, ...) suppressWarnings(
  fit_gwr_model(d, "z", c("a", "b"), kernel = "boxcar", ...))

FOOT_US <- 1200 / 3937
.bx_crs <- function(to_meter) sf::st_crs(sprintf(
  "+proj=utm +zone=17 +datum=WGS84 +to_meter=%s +no_defs",
  format(to_meter, digits = 17)))

test_that("an adaptive boxcar on a lattice gives the same fitted values in any unit, at any origin and in any row order", {
  L <- .bx_lattice()
  ref_layer <- .bx_layer(L, L$xy, 32617)
  perm <- withr::with_seed(9, sample(64))
  for (bw in list(20, NULL)) {
    ref <- .bx_fit(ref_layer, bandwidth = bw)
    variants <- list(
      feet  = .bx_layer(L, L$xy / FOOT_US, .bx_crs(FOOT_US)),
      micro = .bx_layer(L, L$xy / 1e-6, .bx_crs(1e-6)),
      shift = .bx_layer(L, L$xy + 1e8, sf::NA_crs_),
      perm  = .bx_layer(L, L$xy, 32617, perm))
    for (nm in names(variants)) {
      fit <- .bx_fit(variants[[nm]], bandwidth = bw)
      expect_identical(fit$info$bandwidth, ref$info$bandwidth)
      fv <- fitted(fit)
      if (nm == "perm") fv <- fv[order(perm)]
      # The windows hold the same rows, so the fits agree to the rounding
      # of the coordinates' unit, which the local solve does not read.
      expect_lt(max(abs(fv - fitted(ref))), 1e-9)
      nw <- fit$info$local_collinearity$n_window
      if (nm == "perm") nw <- nw[order(perm)]
      expect_identical(nw, ref$info$local_collinearity$n_window)
    }
  }
})

test_that("the boxcar keeps a distance within .TIE_TOL of the edge and drops one beyond it", {
  h <- 250 * sqrt(5)
  d <- c(0, 100, h, h * (1 + 1e-12), h * (1 - 1e-12), h * (1 + 1e-6), 1e3)
  w <- .gw_kernel_weights(d, 3, "boxcar", adaptive = TRUE, h = h)
  expect_identical(w, c(1, 1, 1, 1, 1, 0, 0))
  win <- .gw_window(d, h, "boxcar")
  expect_identical(win$rows, c(1L, 2L, 3L, 4L, 5L))
  expect_identical(win$w, rep(1, 5))
  # The edge is d^2 <= h^2 (1 + .TIE_TOL), so d <= h (1 + 5e-8) to the
  # first order: 4e-8 beyond h is kept and 6e-8 beyond it is dropped, at
  # any size of h, where d^2 itself would overflow or vanish.
  for (hh in c(h, 1e200, 1e-200)) {
    de <- hh * c(1 + 4e-8, 1 + 6e-8)
    expect_identical(.gw_kernel_weights(de, 3, "boxcar", adaptive = TRUE,
                                        h = hh), c(1, 0))
    expect_identical(.gw_window(de, hh, "boxcar")$rows, 1L)
  }
  # A zero-width kernel keeps the rows at distance 0 only, as before.
  expect_identical(.gw_window(c(0, 0, 1e-300, 1), 0, "boxcar")$rows, 1:2)
  # Distances away from the edge give the rows !(d > h) gives.
  dd <- withr::with_seed(2, runif(500, 0, 10))
  hh <- sort(dd)[100]
  expect_identical(.gw_window(dd, hh, "boxcar")$rows, which(!(dd > hh)))
  # The other truncated kernels are unchanged at the edge.
  expect_identical(.gw_kernel_weights(h * (1 + 1e-12), 3, "bisquare",
                                      adaptive = TRUE, h = h), 0)
})
