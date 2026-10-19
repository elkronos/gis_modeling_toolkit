# tests/testthat/test-fix-variogram-reml.R
# ---------------------------------------------------------------------------
# detrend = "reml" keeps the best of its starts, by REML log-likelihood.
#
# .reml_trend() returned the first start that gave finite parameters and never
# compared likelihoods.  Its starts were a tenth, a thirtieth and a third of
# the extent, used as the range parameter; when the true range is short
# against the extent they lie in a valley of the likelihood, and which side
# the optimiser climbs is arbitrary.  On the layer below (240 clustered sites
# over 600 km, effective range 10 km, a linear trend) the first start ends at
# an effective range of 502 km with a nugget proportion of 0.96, at a
# log-likelihood 1.18 below the fit the second start reaches (4.4 km), and
# estimate_sac_range() refused it as past the largest lag.
# ---------------------------------------------------------------------------

# Clustered sites over a 600 km extent, a 10 km exponential field with a
# nugget of 0.4 on a sill of 1, and a linear trend in x.
fvr_draw <- function(i, n_cl = 12L, per = 15L, bg = 60L, w = 600e3) {
  set.seed(500 + i)
  cen <- cbind(runif(n_cl, 0.1 * w, 0.9 * w), runif(n_cl, 0.1 * w, 0.9 * w))
  id <- rep(seq_len(n_cl), each = per)
  xy <- rbind(cen[id, ] + matrix(rnorm(2 * length(id), sd = w / 60), ncol = 2),
              cbind(runif(bg, 0, w), runif(bg, 0, w)))
  n <- nrow(xy)
  set.seed(900 + i)
  D <- as.matrix(stats::dist(xy))
  z <- drop(crossprod(chol(0.6 * exp(-3 * D / 10e3) + diag(1e-8, n)), rnorm(n))) +
    rnorm(n, 0, sqrt(0.4)) + 2 * xy[, 1] / w
  d <- data.frame(x = xy[, 1], y = xy[, 2], z = z, X1 = xy[, 1] / w)
  list(d = d, pts = sf::st_as_sf(d, coords = c("x", "y"), crs = 32615, remove = FALSE),
       extent = max(diff(range(xy[, 1])), diff(range(xy[, 2]))))
}
# The REML log-likelihood gls() reaches from one starting range, as the old
# code started it.
fvr_ll_from <- function(dr, r0) {
  d <- dr$d; d$.sac_x <- d$x; d$.sac_y <- d$y
  f <- suppressWarnings(nlme::gls(
    z ~ X1, data = d, method = "REML",
    correlation = nlme::corExp(value = c(r0, 0.1), form = ~ .sac_x + .sac_y, nugget = TRUE)))
  cs <- stats::coef(f$modelStruct$corStruct, unconstrained = FALSE)
  list(logLik = as.numeric(stats::logLik(f)), eff_range = 3 * cs[["range"]],
       nugget = cs[["nugget"]])
}

test_that("the REML fit kept is the best of the starts, not the first that converges", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  dr <- fvr_draw(7)
  first  <- fvr_ll_from(dr, dr$extent / 10)     # the start the old code kept
  second <- fvr_ll_from(dr, dr$extent / 30)
  # The premise: the first start ends lower, at a range beyond the extent's
  # scale, and the second finds the short range.
  expect_gt(second$logLik, first$logLik + 0.5)
  expect_gt(first$eff_range, 20 * second$eff_range)

  # structure_alpha = 1: the field's structure is weak on these sites (the
  # correlation adds 1.5 to the log-likelihood), the test for spatial
  # structure does not find it, and this test is about which REML fit is
  # kept when one is made.
  r <- suppressWarnings(estimate_sac_range(dr$pts, "z", "X1", detrend = "reml",
                                           structure_alpha = 1))
  expect_identical(attr(r, "detrend_method"), "reml")
  info <- attr(r, "reml")
  # At least as good as either start, and clearly better than the first.
  expect_gte(info$logLik, max(first$logLik, second$logLik) - 0.01)
  expect_gt(info$logLik, first$logLik + 0.5)
  # The range is the short one (the field's is 10 km), not hundreds of km.
  fitted <- if (is.finite(r)) as.numeric(r) else attr(r, "rejected_range")
  expect_lt(fitted, 50e3)
  expect_gt(fitted, 1e3)
  expect_false(identical(attr(r, "rejected_reason"), "fitted range exceeds the largest lag fitted"))
  # The model attached is that fit's.
  vm <- attr(r, "variogram_model")
  expect_equal(3 * vm$range[2], fitted)
  expect_equal(sum(vm$psill), info$sigma2, tolerance = 1e-8)
})

test_that("the REML record says how good the fit is and how many starts it took", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  dr <- fvr_draw(7)
  r <- suppressWarnings(estimate_sac_range(dr$pts, "z", "X1", detrend = "reml"))
  # The REML record rides on the result whatever the verdict on the range.
  info <- attr(r, "reml")
  expect_true(all(c("n_used", "subsampled", "nugget_prop", "sigma2", "logLik",
                    "logLik_gain", "n_fits") %in% names(info)))
  expect_identical(info$n_used, nrow(dr$pts))
  expect_true(is.finite(info$logLik))
  # The gain over the same trend with independent errors.  The first start,
  # which used to be kept, gains 0.35; the fit kept gains 1.5.
  ll0 <- as.numeric(stats::logLik(nlme::gls(z ~ X1, data = dr$d, method = "REML")))
  expect_equal(info$logLik_gain, info$logLik - ll0, tolerance = 1e-6)
  expect_gt(info$logLik_gain, 1)
  expect_true(info$n_fits >= 1L && info$n_fits <= 4L)
})

test_that(".reml_trend() starts from the likelihood and does not discard a nugget near zero", {
  skip_if_not_installed("nlme")
  dr <- fvr_draw(7)
  mf <- dr$d[, c("z", "X1")]; mf$.sac_x <- dr$d$x; mf$.sac_y <- dr$d$y
  fit <- spatialkit:::.reml_trend(mf, z ~ X1, extent = dr$extent)
  expect_type(fit, "list")
  expect_true(all(c("beta", "range", "nugget_prop", "sigma2", "n_used", "subsampled",
                    "pair_floor", "logLik", "logLik_null", "n_fits") %in% names(fit)))
  # One fit, from the best point of the likelihood grid, reaches the optimum.
  expect_identical(fit$n_fits, 1L)
  first <- fvr_ll_from(dr, dr$extent / 10)
  expect_gt(fit$logLik, first$logLik + 0.5)
  expect_lt(3 * fit$range, 50e3)
  # The optimum has almost no nugget here.  A fit is judged by its
  # likelihood, not thrown out for where its nugget lies.
  expect_lt(fit$nugget_prop, 0.2)
  # The same rows in another order give the same fit.
  set.seed(1); perm <- sample(nrow(mf))
  fit2 <- spatialkit:::.reml_trend(mf[perm, ], z ~ X1, extent = dr$extent)
  expect_equal(fit2$logLik, fit$logLik, tolerance = 1e-6)
  expect_equal(fit2$range, fit$range, tolerance = 1e-3)
})
