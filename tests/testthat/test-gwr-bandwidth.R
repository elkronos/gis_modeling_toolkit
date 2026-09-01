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


# ---------------------------------------------------------------------------
# "Exactly two distinct values" is not the same thing as "binary"
# ---------------------------------------------------------------------------

.gwr2v_points <- function(z, n = 60, seed = 2) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), z = rep(z, length.out = n)),
    coords = c("x", "y"), crs = 32632)
}

test_that("fit_gwr_model fits a two-valued NON-INTEGER response with a warning", {
  # A left-censored measurement -- everything either at the detection limit
  # (0.0031) or at a single higher reading -- has exactly two distinct values
  # and is perfectly continuous.  Gaussian GWR on it is a well-defined
  # least-squares problem, and the old error's advice to switch to
  # family = "binomial" was nonsense for such values.  The same guard runs once
  # per fold inside cv_gwr(), where a small training fold can legitimately hold
  # only two distinct values, so a hard stop there aborted whole CV runs.
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  censored <- .gwr2v_points(c(0.0031, 12.7401))
  expect_equal(length(unique(censored$z)), 2L)
  expect_false(all(censored$z == round(censored$z)))

  expect_warning(fit <- fit_gwr_model(censored, "z", "a", bandwidth = 300),
                 "has only 2 distinct finite values")
  expect_s3_class(fit, "gwr_fit")
  # Warned, not refused: it is the design that is degenerate, not the model.
  expect_warning(fit_gwr_model(censored, "z", "a", bandwidth = 300),
                 "genuinely continuous \\(e.g. censored at a detection limit\\)")

  # An integer-valued pair is still a hard error, with the binomial advice.
  binary <- .gwr2v_points(c(0, 1))
  expect_true(all(binary$z == round(binary$z)))
  expect_error(fit_gwr_model(binary, "z", "a", bandwidth = 300),
               "is binary \\(2 distinct values")
  expect_error(fit_gwr_model(binary, "z", "a", bandwidth = 300),
               "family = 'binomial'")
})

test_that("the two-valued response guard sits behind the GWmodel requirement", {
  # Dependency-free companion to the skipped test above, and the reason it has
  # to skip: fit_gwr_model() checks for its backend BEFORE it looks at the
  # response, so nothing about the response can be observed without GWmodel.
  # Both response shapes therefore produce the same missing-backend error, and
  # neither the warning nor the binary stop is reachable on this install.
  skip_if(requireNamespace("GWmodel", quietly = TRUE),
          "GWmodel is installed, so the backend check passes")

  for (z in list(c(0.0031, 12.7401), c(0, 1), c(1.5, 2.5, 3.5))) {
    expect_error(fit_gwr_model(.gwr2v_points(z), "z", "a", bandwidth = 300),
                 "package 'GWmodel' is required")
  }
})
