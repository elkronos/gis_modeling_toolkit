# tests/testthat/test-gp-basis.R
# ---------------------------------------------------------------------------
# Specification for .gp_basis_spec(), which replaced the old
# min(n/3, max(15, sqrt(n))) heuristic.
#
# brms carries k^D basis functions for gp(..x, ..y, k = k), so D = 2 and the
# model cost is k^2.  The old rule made k ~ sqrt(n), hence k^2 == n: the
# approximation stopped approximating.  These tests pin the properties that
# make the replacement correct.
# ---------------------------------------------------------------------------

scaled_uniform <- function(n, seed = 1) {
  set.seed(seed)
  scale(matrix(runif(2 * n), ncol = 2))
}

test_that(".gp_basis_spec does not scale with n", {
  # The whole point of the change.  Resolution should track spatial structure,
  # not row count.
  #
  # Not asserted as exact equality: gp_lengthscale_bounds() subsamples through
  # .safe_dist(max_n = 1000), so above that size the 25th-percentile distance
  # is estimated rather than exact, and ceiling() can turn a tiny difference
  # into one basis function.  That subsample is itself seeded, so this is
  # deterministic -- just not identical across sizes.
  ns <- c(200, 1000, 10000)
  ks <- vapply(ns, function(n) {
    xy <- scaled_uniform(n)
    spatialkit:::.gp_basis_spec(xy, gp_lengthscale_bounds(xy))$k
  }, integer(1))

  expect_lte(diff(range(ks)), 1L)          # flat to within rounding
  expect_lt(max(ks) / min(ks), 1.1)        # ... and flat proportionally

  # The old rule gave 15 / 31 / 100 over these sizes -- a 6.7x spread, with
  # gp_k^2 == n by construction.  Guard against regressing to that shape.
  expect_lt(ks[3], 40L)
  expect_lt(ks[3]^2, ns[3] / 4)
})

test_that(".gp_basis_spec is stable across point-pattern geometry", {
  set.seed(2)
  unif <- scale(matrix(runif(2 * 2000), ncol = 2))
  elon <- scale(cbind(runif(2000, 0, 10), runif(2000, 0, 0.5)))
  ctr  <- matrix(runif(2 * 20), ncol = 2)
  lab  <- sample(20, 2000, replace = TRUE)
  clus <- scale(ctr[lab, ] + matrix(rnorm(2 * 2000, 0, 0.02), ncol = 2))

  ks <- vapply(list(unif, elon, clus), function(xy)
    spatialkit:::.gp_basis_spec(xy, gp_lengthscale_bounds(xy))$k, integer(1))

  # Derived k sits around 21-25 for realistic patterns; the band is loose
  # enough to survive RNG jitter but tight enough to catch a regression to an
  # n-dependent rule.
  expect_true(all(ks >= 15L & ks <= 35L))
})

test_that(".gp_basis_spec keeps the total basis count modest", {
  # The failure being guarded against: 10,000 basis functions at n = 10,000.
  xy <- scaled_uniform(10000)
  s  <- spatialkit:::.gp_basis_spec(xy, gp_lengthscale_bounds(xy))
  expect_lt(s$k^2, 2500L)
})

test_that(".gp_basis_spec respects its floor and its cap", {
  xy <- scaled_uniform(2000, seed = 3)
  b  <- gp_lengthscale_bounds(xy)

  # Defaults do not bind for realistic data, so exercise the clamps directly.
  expect_gte(spatialkit:::.gp_basis_spec(xy, b, k_min = 40L)$k, 40L)

  tight <- spatialkit:::.gp_basis_spec(xy, b, max_basis = 100L)
  expect_lte(tight$k^2, 100L)
  expect_true(tight$capped)

  loose <- spatialkit:::.gp_basis_spec(xy, b)
  expect_false(loose$capped)
})

test_that(".gp_basis_spec sets c from the UPPER length-scale bound", {
  # The boundary must contain the longest plausible correlation range.  The
  # previous hard-coded gp_c = 1.5 fails this whenever ell/S > 0.47.
  xy <- scaled_uniform(2000, seed = 4)
  b  <- gp_lengthscale_bounds(xy)
  s  <- spatialkit:::.gp_basis_spec(xy, b)

  # Compute in the same operation order as .gp_basis_spec() -- 3.2 * (upper/S),
  # not (3.2*upper)/S -- and allow a rounding tolerance: when this constraint
  # binds, c_val IS this product, so a strict >= is a floating-point coin flip.
  expect_gte(s$c, 3.2 * (b[["upper"]] / s$S) - 1e-9)
  expect_gte(s$c, 1.2)
  expect_gt(s$c, 1.5)          # would have failed under the old default
})

test_that(".gp_basis_spec returns the domain half-range it used", {
  xy <- scaled_uniform(500, seed = 5)
  s  <- spatialkit:::.gp_basis_spec(xy, gp_lengthscale_bounds(xy))
  expect_equal(s$S, max(apply(xy, 2, function(z) diff(range(z)) / 2)))
})

test_that(".gp_basis_spec survives degenerate coordinates", {
  xy <- cbind(rep(0, 10), rep(0, 10))          # zero extent -> S guard
  s  <- spatialkit:::.gp_basis_spec(xy, c(lower = 0.1, upper = 1))
  expect_true(is.finite(s$k) && s$k >= 10L)
  expect_true(is.finite(s$c))
  expect_equal(s$S, 1)
})

test_that("gp_lengthscale_bounds returns a sane, separated interval", {
  xy <- scaled_uniform(500, seed = 6)
  b  <- gp_lengthscale_bounds(xy)

  expect_named(b, c("lower", "upper"))
  expect_true(all(is.finite(b)))
  expect_gt(b[["lower"]], 0)
  expect_gte(b[["upper"]], b[["lower"]] * 1.2)   # the separation guard
})


test_that("the GP term disables brms's own covariate scaling", {
  # brms::gp() defaults to scale = TRUE, which normalises covariates so the max
  # pairwise distance is 1 and reports lscale in that space.  This package
  # already standardises coordinates, and gp_lengthscale_bounds(), gp_c and
  # gp_ell_min are all in those units -- so a second normalisation silently puts
  # the length-scale prior and the adequacy check in the wrong units (off by the
  # max pairwise distance, ~4.9 for standardised 2D coords).  That made the
  # adequacy warning fire on 100% of draws for every fit.
  term <- spatialkit:::.gp_formula_term(24, 3.5)
  expect_match(term, "scale = FALSE", fixed = TRUE)
})

test_that("the GP fits one length-scale per axis by default", {
  # Coordinates are standardised per axis, so a single shared length-scale
  # would make the kernel anisotropic in the original CRS by the ratio
  # sd(X)/sd(Y) -- a property of the sampling layout, not of the process.
  expect_match(spatialkit:::.gp_formula_term(24, 3.5), "iso = FALSE", fixed = TRUE)
  expect_match(spatialkit:::.gp_formula_term(24, 3.5, gp_iso = TRUE),
               "iso = TRUE", fixed = TRUE)
})

test_that("the GP term carries k and c through verbatim", {
  term <- spatialkit:::.gp_formula_term(17, 2.25)
  expect_match(term, "k = 17", fixed = TRUE)
  expect_match(term, "c = 2.25", fixed = TRUE)
})

test_that("the GP term is syntactically valid R", {
  # A malformed term would only surface at brm() time, after Stan compilation.
  for (iso in c(TRUE, FALSE)) {
    term <- spatialkit:::.gp_formula_term(24, 3.5, gp_iso = iso)
    f <- eval(parse(text = sprintf("y ~ x + %s", term))[[1]])
    expect_s3_class(f, "formula")
    expect_true(any(grepl("^gp\\(", attr(stats::terms(f), "term.labels"))))
  }
})
