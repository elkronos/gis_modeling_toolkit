# tests/testthat/test-lscale-prior.R
# ---------------------------------------------------------------------------
# .lscale_invgamma(): calibration of the GP length-scale prior.
#
# The previous prior was normal(0, sd) on a positive parameter -- a half-normal
# with its mode at zero, so most of its mass sat at length-scales shorter than
# the basis can resolve.  That is where the Hilbert-space GP funnels and the
# sampler divergences.  An inverse-gamma with both tails pinned to the
# estimated bounds is the standard remedy.
#
# All of this is testable without Stan: it is a numerical calibration.
# ---------------------------------------------------------------------------

ig <- function(...) spatialkit:::.lscale_invgamma(...)

# P(X < q) for X ~ InvGamma(shape, scale), via 1/X ~ Gamma(shape, rate = scale)
p_invgamma <- function(q, shape, scale) 1 - stats::pgamma(1 / q, shape = shape, rate = scale)


test_that("both tails land on the requested probability", {
  # These are the bounds gp_lengthscale_bounds() actually produces on
  # standardised coordinates.
  for (b in list(c(0.4654, 1.8450), c(0.4657, 1.9560), c(0.4549, 1.9300))) {
    s <- ig(b[1], b[2])
    expect_true(s$ok)
    expect_equal(p_invgamma(b[1], s$shape, s$scale), 0.01, tolerance = 0.005)
    expect_equal(1 - p_invgamma(b[2], s$shape, s$scale), 0.01, tolerance = 0.005)
  }
})

test_that("the prior concentrates between the bounds, not at zero", {
  # The whole point: mass must be pulled away from the short length-scales
  # where the approximation degenerates.
  s <- ig(0.4654, 1.8450)
  expect_true(s$ok)

  # Under a half-normal with the old sd, a large share of the mass sits below
  # the lower bound; under the calibrated inverse-gamma it is 1% by design.
  half_normal_below <- 2 * (stats::pnorm(0.4654, 0, 1.845) - 0.5)
  expect_gt(half_normal_below, 0.15)
  expect_lt(p_invgamma(0.4654, s$shape, s$scale), 0.02)

  # Median should sit inside the bounds.
  med <- 1 / stats::qgamma(0.5, shape = s$shape, rate = s$scale)
  expect_gt(med, 0.4654)
  expect_lt(med, 1.8450)
})

test_that("it declines rather than returning a sloppy fit", {
  # A badly calibrated inverse-gamma would be worse than the half-normal it
  # replaces, so the caller needs an honest ok = FALSE to fall back on.
  expect_false(ig(1, 1)$ok)          # degenerate: lower == upper
  expect_false(ig(2, 1)$ok)          # inverted
  expect_false(ig(-1, 5)$ok)         # non-positive
  expect_false(ig(0, 5)$ok)
  expect_false(ig(NA_real_, 5)$ok)
  expect_false(ig(0.1, Inf)$ok)
})

test_that("shape and scale are finite and positive when ok", {
  s <- ig(0.2, 1.4)
  expect_true(s$ok)
  expect_true(is.finite(s$shape) && s$shape > 0)
  expect_true(is.finite(s$scale) && s$scale > 0)
})

test_that("a tighter tail target pushes more mass inside the bounds", {
  loose <- ig(0.4, 2.0, tail = 0.05)
  tight <- ig(0.4, 2.0, tail = 0.001)
  expect_true(loose$ok); expect_true(tight$ok)
  expect_lt(p_invgamma(0.4, tight$shape, tight$scale),
            p_invgamma(0.4, loose$shape, loose$scale))
})

test_that("the calibrated prior is a valid Stan inv_gamma string", {
  s <- ig(0.4654, 1.8450)
  spec <- sprintf("inv_gamma(%.6f, %.6f)", s$shape, s$scale)
  expect_match(spec, "^inv_gamma\\([0-9.]+, [0-9.]+\\)$")
  # parses as a call with two positive numeric arguments
  parsed <- parse(text = spec)[[1]]
  expect_equal(as.character(parsed[[1]]), "inv_gamma")
  expect_true(all(vapply(as.list(parsed)[-1], function(z) is.numeric(z) && z > 0,
                         logical(1))))
})
