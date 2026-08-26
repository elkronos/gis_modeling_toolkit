# ===========================================================================
# GWR bandwidth selection and its fallback.
#
# The fallback bandwidth has no relationship to the data's spatial
# structure, so the thing that matters is that a fit which used it says so
# ($info$bandwidth_is_fallback) rather than presenting an arbitrary number
# as a selected one.
# ===========================================================================

test_that(".fallback_bandwidth adaptive clamp respects [10, 0.9n]", {
  skip_if_not_installed("sp")
  fb <- spatialkit:::.fallback_bandwidth

  make_sp <- function(n) {
    sp::SpatialPointsDataFrame(
      cbind(seq_len(n), seq_len(n)),
      data = data.frame(x = seq_len(n))
    )
  }

  # Tiny n: lower clamp of 10 applies (old code returned floor(0.9n) < 10)
  expect_equal(fb(make_sp(3), adaptive = TRUE), 10L)
  # Mid n: 0.9n binds below the 50-neighbour target
  expect_equal(fb(make_sp(30), adaptive = TRUE), 27L)
  # Large n: the ~50-neighbour target binds
  expect_equal(fb(make_sp(200), adaptive = TRUE), 50L)
})


test_that(".fallback_bandwidth fixed branch is one third of the bbox diagonal", {
  skip_if_not_installed("sp")
  fb <- spatialkit:::.fallback_bandwidth

  # A 300 x 400 bounding box: diagonal 500 by construction, so the documented
  # rule gives exactly 500/3.  Asserting the value rather than `> 0` is what
  # separates this from any positive number a broken implementation returns.
  sp_dat <- sp::SpatialPointsDataFrame(
    cbind(c(100, 400, 250), c(50, 450, 200)),
    data = data.frame(x = 1:3)
  )
  expect_equal(fb(sp_dat, adaptive = FALSE), 500 / 3)

  # Degenerate extent: all points coincident gives a zero diagonal, which must
  # come back as the positive epsilon floor rather than 0 (a zero bandwidth
  # makes every GWmodel kernel weight undefined).
  coincident <- sp::SpatialPointsDataFrame(
    cbind(rep(7, 3), rep(9, 3)), data = data.frame(x = 1:3)
  )
  expect_equal(fb(coincident, adaptive = FALSE), .Machine$double.eps)

  # The two branches are genuinely different code paths, not one rounded:
  # the adaptive branch counts neighbours, the fixed branch measures distance.
  expect_false(isTRUE(all.equal(fb(sp_dat, adaptive = TRUE),
                                fb(sp_dat, adaptive = FALSE))))
})


test_that("a successful bandwidth selection is not labelled a fallback", {
  # The contract that matters downstream: compare_models() warns, and
  # $info$bandwidth_is_fallback drives that warning.  Assert it on a real fit
  # -- a mock `info` list literal can only re-assert its own contents.
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  set.seed(4)
  n <- 80
  px <- runif(n, 0, 5000); py <- runif(n, 0, 5000)
  dat <- sf::st_as_sf(
    data.frame(x1 = rnorm(n), px = px, py = py,
               y = 2 + 0.001 * px + 0.5 * rnorm(n)),
    coords = c("px", "py"), crs = 32632
  )

  # An explicit bandwidth is never a fallback, and is carried through verbatim.
  explicit <- fit_gwr_model(dat, "y", "x1", bandwidth = 30, adaptive = TRUE)
  expect_false(explicit$info$bandwidth_is_fallback)
  expect_equal(explicit$info$bandwidth, 30)

  # Automatic selection on well-behaved data succeeds, so the flag stays FALSE
  # and the chosen bandwidth is a real selection rather than the arbitrary
  # fallback constant.
  auto <- fit_gwr_model(dat, "y", "x1", adaptive = TRUE)
  expect_false(auto$info$bandwidth_is_fallback)
  expect_true(is.finite(auto$info$bandwidth) && auto$info$bandwidth > 0)
})
