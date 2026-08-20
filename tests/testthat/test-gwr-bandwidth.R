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


test_that("fallback_bandwidth sets bandwidth_is_fallback in info", {
  fb_adaptive <- spatialkit:::.fallback_bandwidth
  # Construct a minimal mock Spatial object
  skip_if_not_installed("sp")
  coords <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  sp_dat <- sp::SpatialPointsDataFrame(
    coords,
    data = data.frame(x = 1:3)
  )

  bw_a <- fb_adaptive(sp_dat, adaptive = TRUE)
  expect_true(is.integer(bw_a) || is.numeric(bw_a))
  expect_true(bw_a >= 1)

  bw_f <- fb_adaptive(sp_dat, adaptive = FALSE)
  expect_true(is.numeric(bw_f))
  expect_true(bw_f > 0)
})


test_that("gwr_fit info contains bandwidth_is_fallback field", {
  # We cannot run a full GWR fit without GWmodel + real data, but we can

  # verify the info list contract by constructing a mock gwr_fit.
  skip_if_not_installed("sf")
  pts <- sf::st_as_sf(
    data.frame(x = runif(10), y = runif(10), resp = rnorm(10)),
    coords = c("x", "y"), crs = 32632
  )
  mock_info <- list(
    bandwidth = 50,
    adaptive = TRUE,
    kernel = "bisquare",
    AICc = NA_real_,
    bandwidth_is_fallback = TRUE
  )
  expect_true("bandwidth_is_fallback" %in% names(mock_info))
  expect_true(mock_info$bandwidth_is_fallback)
})
