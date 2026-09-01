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

test_that(".gp_basis_spec follows the Riutort-Mayol inequalities exactly", {
  # The documented specification, not an observed band:
  #     c >= 3.2 * (ell/S),  c >= 1.2
  #     m >= 1.75 * c / (ell/S)
  # with the boundary factor set from the UPPER bound (it must contain the
  # longest plausible range) and the basis count from the LOWER one (it must
  # resolve the shortest).  m is an integer, so the second inequality is met by
  # ceiling(): the slack is in [0, 1), never more.
  set.seed(2)
  unif <- scale(matrix(runif(2 * 2000), ncol = 2))
  elon <- scale(cbind(runif(2000, 0, 10), runif(2000, 0, 0.5)))
  ctr  <- matrix(runif(2 * 20), ncol = 2)
  lab  <- sample(20, 2000, replace = TRUE)
  clus <- scale(ctr[lab, ] + matrix(rnorm(2 * 2000, 0, 0.02), ncol = 2))

  for (xy in list(unif, elon, clus)) {
    b <- gp_lengthscale_bounds(xy)
    s <- spatialkit:::.gp_basis_spec(xy, b)

    r_lo <- b[["lower"]] / s$S          # ell/S that must be RESOLVED
    r_hi <- b[["upper"]] / s$S          # ell/S that must be CONTAINED

    expect_gte(s$c, 3.2 * r_hi - 1e-9)
    expect_gte(s$c, 1.2)
    # Neither clamp binds for these patterns, so k is exactly the ceiling of
    # the second inequality -- a stronger statement than ">=".
    expect_false(s$capped)
    expect_identical(s$k, as.integer(ceiling(1.75 * s$c / r_lo)))
  }
})

test_that(".gp_basis_spec tracks the length-scale ratio, not the coordinate scale", {
  # Composing the two inequalities (in the regime where c = 3.2 * r_hi, i.e.
  # whenever 3.2 * upper > 1.2 * S) gives
  #     k = ceiling(1.75 * 3.2 * upper / lower) = ceiling(5.6 * upper / lower)
  # -- S cancels.  k therefore depends ONLY on the shape of the pairwise
  # distance distribution, and is invariant to multiplying every coordinate by
  # a constant.  An implementation that reached for n, for S, or for an
  # absolute distance anywhere would break this.
  set.seed(2)
  xy <- scale(matrix(runif(2 * 2000), ncol = 2))
  b  <- gp_lengthscale_bounds(xy)
  k0 <- spatialkit:::.gp_basis_spec(xy, b)$k

  expect_identical(k0, as.integer(ceiling(5.6 * b[["upper"]] / b[["lower"]])))

  # S is brms's own domain measure: the pooled range of the COLUMN-CENTRED
  # coordinates, which is exactly what brms:::choose_L() multiplies by c.
  # Verified against brms itself -- the boundary recovered from
  # make_standata()'s slambda equals c * this quantity to the last digit at
  # every c, and equals TWICE the per-axis half-range the earlier
  # implementation used.
  brms_S <- function(m) {
    Xc <- sweep(as.matrix(m), 2L, colMeans(as.matrix(m)))
    max(1, max(Xc) - min(Xc))                      # the max(1, .) is brms's
  }
  # Scale invariance of k holds wherever brms's own max(1, .) floor does not
  # bite -- which is every realistic input, since this package hands gp() the
  # per-axis standardised coordinates (pooled centred range ~3.5 here).
  for (mult in c(1, 1e3, 1e6)) {
    xy_s <- xy * mult
    s_s  <- spatialkit:::.gp_basis_spec(xy_s, gp_lengthscale_bounds(xy_s))
    expect_gt(brms_S(xy_s), 1)                     # floor not engaged
    expect_identical(s_s$k, k0)                    # k unchanged ...
    expect_equal(s_s$S, brms_S(xy_s))              # ... while S scales with it
  }
  # Regression guard: the half-range convention is exactly half of this, and
  # using it built a GP boundary twice as wide as gp_k was sized for.
  expect_equal(spatialkit:::.gp_basis_spec(xy, b)$S, brms_S(xy))
  expect_gt(brms_S(xy), 1.5 * max(apply(xy, 2, function(z) diff(range(z)) / 2)))
})

test_that(".gp_basis_spec reproduces brms's max(1, range) boundary floor", {
  # choose_L() floors the domain range at 1, so below that extent brms builds
  # the boundary from 1 rather than from the data.  This function must do the
  # same, or gp_k and gp_ell_min would describe a basis brms does not build.
  # It is why k is NOT scale-invariant all the way down.
  set.seed(2)
  xy    <- scale(matrix(runif(2 * 2000), ncol = 2))
  tiny  <- xy * 1e-3                                # pooled centred range ~0.0035
  s_t   <- spatialkit:::.gp_basis_spec(tiny, gp_lengthscale_bounds(tiny))
  expect_equal(s_t$S, 1)                            # floored, not 0.0035
  expect_equal(s_t$c, 1.25)                         # so c falls to its floor
  # and k rises, because the basis must now resolve a much smaller ell/S.
  expect_gt(s_t$k, spatialkit:::.gp_basis_spec(xy, gp_lengthscale_bounds(xy))$k)
})

test_that(".gp_basis_spec demands more basis functions for finer structure", {
  # Three tight clusters far apart: a third of all pairs are within-cluster, so
  # the 25th-percentile distance is a within-cluster one and upper/lower jumps
  # from ~4 to ~90.  The required k rises with it and hits the cap -- which is
  # the behaviour a length-scale-driven rule must have and an n-driven one
  # cannot, since n is the same order as the uniform case above.
  set.seed(9)
  centres <- matrix(c(0, 0, 10, 0, 5, 9), ncol = 2, byrow = TRUE)
  tight   <- centres[sample(3, 1500, replace = TRUE), ] +
    matrix(rnorm(2 * 1500, 0, 0.05), ncol = 2)

  b <- gp_lengthscale_bounds(tight)
  s <- spatialkit:::.gp_basis_spec(tight, b)

  expect_gt(b[["upper"]] / b[["lower"]], 20)       # genuinely finer structure
  expect_gt(1.75 * s$c / (b[["lower"]] / s$S), 50) # more than the cap allows
  expect_true(s$capped)
  expect_identical(s$k, as.integer(floor(sqrt(2500L))))
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
  expect_gte(s$c, 1.25)        # brms's own default boundary factor, c = 5/4
  expect_gt(s$c, 1.5)          # would have failed under the old default
})

test_that(".gp_basis_spec returns brms's own domain measure, not the half-range", {
  # brms builds the GP boundary as choose_L(x, c) = c * max(1, max(x) - min(x))
  # over the column-centred, POOLED covariate matrix.  Deriving c against the
  # per-axis half-range and handing it to brms::gp() therefore produced a
  # boundary twice as wide as intended, under-resolving the GP by a factor of
  # two and making $info$gp_ell_min twice too lenient to catch it.
  xy <- scaled_uniform(500, seed = 5)
  s  <- spatialkit:::.gp_basis_spec(xy, gp_lengthscale_bounds(xy))
  Xc <- sweep(xy, 2L, colMeans(xy))
  expect_equal(s$S, max(1, max(Xc) - min(Xc)))
  # ... and that is twice the old half-range, which is the size of the bug.
  expect_equal(s$S / max(apply(xy, 2, function(z) diff(range(z)) / 2)), 2,
               tolerance = 0.05)
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

test_that("gp_lengthscale_bounds calibrates to 5 % SE correlation at the bound", {
  # Finiteness, positivity and a 1.2x gap are satisfied by an interval that is
  # 2.45x too wide at both ends, so none of them pins the calibration.  The
  # documented rule does: the squared-exponential kernel has
  #
  #     corr(d) = exp(-d^2 / (2 l^2)),
  #
  # so requiring corr = 0.05 at distance d gives d^2 / (2 l^2) = log 20, i.e.
  # l = d / sqrt(2 log 20).  The lower bound applies that to the q_small
  # quantile of the pairwise distances and the upper bound to the maximum, so
  # BOTH bounds must reproduce a correlation of exactly 0.05 at their own
  # distance.  Nothing here re-uses the implementation's constant.
  xy <- scaled_uniform(500, seed = 6)
  b  <- gp_lengthscale_bounds(xy)

  d    <- as.numeric(stats::dist(xy))          # n = 500 <= max_n, so no sample
  d    <- d[d > 0]
  dq   <- stats::quantile(d, 0.25, names = FALSE, type = 7)   # q_small default
  dmax <- max(d)
  # The separation guard must not be the thing setting `upper`, or the upper
  # assertion below would be testing max() rather than the calibration.
  expect_gt(dmax, 1.2 * dq)

  se_corr <- function(dist, l) exp(-dist^2 / (2 * l^2))
  expect_equal(se_corr(dq,   b[["lower"]]), 0.05, tolerance = 1e-12)
  expect_equal(se_corr(dmax, b[["upper"]]), 0.05, tolerance = 1e-12)

  # Equivalently, and more directly: the bounds ARE those distances divided by
  # sqrt(2 log 20) ~ 2.448.  A factor of 1 -- no conversion at all -- would put
  # both bounds 2.45x too high and make the correlation above exp(-1/2) = 0.61.
  expect_equal(unname(b[["lower"]]), dq   / sqrt(2 * log(20)), tolerance = 1e-12)
  expect_equal(unname(b[["upper"]]), dmax / sqrt(2 * log(20)), tolerance = 1e-12)

  # The separation guard still binds when the distance distribution is narrow
  # enough that dmax < 1.2 * dq -- points on one tight ring, say.
  ang  <- seq(0, 2 * pi, length.out = 60)[-60]
  ring <- cbind(cos(ang), sin(ang))
  br   <- gp_lengthscale_bounds(ring, q_small = 0.9)
  dr   <- as.numeric(stats::dist(ring)); dr <- dr[dr > 0]
  expect_lt(max(dr), 1.2 * stats::quantile(dr, 0.9, names = FALSE))
  expect_equal(unname(br[["upper"]]), unname(br[["lower"]]) * 1.2,
               tolerance = 1e-12)
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
