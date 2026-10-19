# The post-fit check of the GP basis against the posterior.
#
# It was one number, the share of length-scale draws below the smallest
# length-scale the basis resolves, logged when above 10 percent and shown
# nowhere else.  It could not fire when the data did not inform the
# length-scale: the posterior then equals the prior, which puts 1 percent
# below that floor by construction.  A Gaussian process that had collapsed
# onto noise over structure finer than the basis therefore passed in silence.
#
# Whether the Gaussian process is identified is decided on the share of the
# variance it carries, where the family has a residual SD to form one, and on
# the length-scales only where it has not.  The width of the length-scale
# posterior against its prior does not decide by itself: a surface smoother
# than the prior expects (a plane, a range longer than the extent) moves the
# length-scale to the long side of its prior, where the posterior is wider
# than the prior.  Three real fits whose GP carried 0.57 to 0.94 of the
# variance at the 2.5 percent quantile were reported as "no spatial structure
# found" on that ground.
#
# The posterior answers neither question about what the surface left behind.
# A Gaussian process that carries a share of the variance at a long range is
# what a regional trend gives, and what a regional trend beside structure
# finer than the basis gives.  The residuals are therefore read for every fit
# whose posterior was, and "ok" is not returned over residuals that are
# autocorrelated at short range.  They used to be read only for a Gaussian
# process that was not identified: a weak plane plus a field with an 8 km
# range under a 33 km floor came back "ok", with no warning of any class,
# over residuals with a Moran's I of 0.143 (p = 5.6e-06).
#
# No Stan.  The engine is a stand-in that carries posterior draws, the brms
# accessors are replaced, and the fitter's own code reads them.

# A brmsfit stand-in holding the draws the adequacy check reads.
.fba_engine <- function(ls, sdgp, sigma) {
  structure(list(ls = ls, sdgp = sdgp, sigma = sigma),
            class = c("fba_draws", "brmsfit"))
}
registerS3method("as.matrix", "fba_draws", function(x, variable = NULL,
                                                    regex = FALSE, ...) {
  if (identical(variable, "^lscale")) return(x$ls)
  if (identical(variable, "^sdgp"))
    return(matrix(x$sdgp, ncol = 1L, dimnames = list(NULL, "sdgp_gp..x..y")))
  if (identical(variable, "sigma"))
    return(matrix(x$sigma, ncol = 1L, dimnames = list(NULL, "sigma")))
  stop("no such variable")
})

.fba_pts <- function(n = 300, seed = 7, range = NULL) {
  set.seed(seed)
  xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  z <- if (is.null(range)) stats::rnorm(n) else {
    D <- as.matrix(stats::dist(xy))
    as.numeric(t(chol(exp(-D / (range / 3)) + diag(1e-8, n))) %*% stats::rnorm(n)) +
      stats::rnorm(n, 0, 0.3)
  }
  sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
               coords = c("x", "y"), crs = 32617)
}

# What the fitter will compute for this layer before sampling: the smallest
# length-scale the derived basis resolves and the prior it sets.
.fba_design <- function(d) {
  xy <- sf::st_coordinates(d)
  zz <- cbind((xy[, 1] - mean(xy[, 1])) / stats::sd(xy[, 1]),
              (xy[, 2] - mean(xy[, 2])) / stats::sd(xy[, 2]))
  b  <- gp_lengthscale_bounds(zz)
  sp <- spatialkit:::.gp_basis_spec(zz, b)
  ig <- spatialkit:::.lscale_invgamma(b[["lower"]], b[["upper"]])
  list(ell_min = 1.75 * sp$c * sp$S / sp$k, shape = ig$shape, scale = ig$scale)
}

# Fit with the stand-in engine and every sampler diagnostic clean.  The fitted
# values are zero, so the residuals the Moran's I step sees are the response.
# `epred` replaces that stand-in for the posterior draw; `sac` is the range
# the variogram before fitting is taken to have returned.
.fba_fit <- function(d, engine,
                     epred = function(object, newdata, ...)
                       matrix(0, nrow = 4L, ncol = nrow(newdata)),
                     sac = 9e4) {
  local_mocked_bindings(
    brm = function(...) engine,
    nuts_params = function(...)
      data.frame(Chain = 1L, Iteration = 1L, Parameter = "divergent__", Value = 0),
    rhat = function(...) c(b_Intercept = 1.001, sigma = 1.002),
    neff_ratio = function(...) c(b_Intercept = 0.9, sigma = 0.8),
    posterior_epred = epred,
    .package = "brms")
  local_mocked_bindings(estimate_sac_range = function(...) sac,
                        .package = "spatialkit")
  w <- list()
  fit <- withCallingHandlers(
    suppressMessages(fit_bayesian_spatial_model(d, "z", character(0),
                                                compute_loo = FALSE)),
    warning = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning") })
  list(fit = fit, warnings = w)
}

.fba_classes <- function(w) vapply(w, function(c) class(c)[1L], character(1))


test_that(".gp_adequacy() separates not identified, below resolution and ok", {
  set.seed(1)
  N <- 2000
  ig_a <- 11.5; ig_b <- 9.7                       # prior 1% point about 0.46
  spec <- sprintf("inv_gamma(%s, %s)", ig_a, ig_b)
  prior_draws <- function() 1 / stats::rgamma(N, shape = ig_a, rate = ig_b)
  ell_min <- 0.46

  # The data did not move the length-scale and the GP carries no variance.
  ni <- spatialkit:::.gp_adequacy(
    cbind(lx = prior_draws(), ly = prior_draws()),
    sdgp_draws = abs(stats::rnorm(N, 0, 0.05)), sigma_draws = rep(1.2, N),
    gp_ell_min = ell_min, prior_spec = spec, xy_scale = c(28000, 29000))
  expect_identical(ni$status, "not_identified")
  expect_false(ni$identified)
  expect_identical(ni$identified_by, "spatial_share")
  expect_true(all(abs(ni$contraction) < 0.25))
  expect_true(all(ni$prior_gap < 0.10))
  expect_lt(ni$spatial_share[["q2.5"]], 0.01)
  # The share below the floor says nothing here: it restates the prior's 1%.
  expect_lt(ni$frac_below, 0.05)

  # Identified, and the posterior wants a shorter length-scale than the basis
  # resolves.
  br <- spatialkit:::.gp_adequacy(
    cbind(lx = stats::rnorm(N, 0.42, 0.03), ly = stats::rnorm(N, 0.50, 0.03)),
    sdgp_draws = stats::rnorm(N, 1, 0.05), sigma_draws = stats::rnorm(N, 0.5, 0.02),
    gp_ell_min = ell_min, prior_spec = spec, xy_scale = c(28000, 29000))
  expect_identical(br$status, "below_resolution")
  expect_true(br$identified)
  expect_gt(br$frac_below, 0.10)
  expect_true(all(br$contraction > 0.9))

  # Identified and clear of the floor.
  ok <- spatialkit:::.gp_adequacy(
    cbind(lx = stats::rnorm(N, 0.9, 0.05), ly = stats::rnorm(N, 1.0, 0.05)),
    sdgp_draws = stats::rnorm(N, 1, 0.05), sigma_draws = stats::rnorm(N, 0.5, 0.02),
    gp_ell_min = ell_min, prior_spec = spec, xy_scale = c(28000, 29000))
  expect_identical(ok$status, "ok")
  expect_identical(dim(ok$effective_range), c(2L, 3L))
  # The effective range per axis in map units: 2.45 length-scales times the
  # axis's own SD.
  expect_equal(unname(ok$effective_range[, "q50"]),
               sqrt(2 * log(20)) * c(28000, 29000) * c(0.9, 1.0), tolerance = 0.01)

  # Without a prior this function set, and without sigma, identification is
  # unknown and the share alone decides.
  un <- spatialkit:::.gp_adequacy(
    cbind(lx = stats::rnorm(N, 0.9, 0.05)), sdgp_draws = stats::rnorm(N, 1, 0.05),
    gp_ell_min = ell_min, prior_spec = "student_t(3, 0, 1)")
  expect_true(is.na(un$identified))
  expect_true(is.na(un$identified_by))
  expect_true(all(is.na(un$prior_gap)))
  expect_identical(un$status, "ok")
  expect_null(un$spatial_share)
})

test_that(".gp_prior_gap() is the Kolmogorov-Smirnov distance to the prior", {
  set.seed(11)
  a <- 12.08121194; b <- 10.05144272
  cdf <- function(x) 1 - stats::pgamma(1 / x, shape = a, rate = b)
  x <- stats::rlnorm(500, log(1.4), 0.3)
  expect_equal(spatialkit:::.gp_prior_gap(x, cdf),
               unname(suppressWarnings(stats::ks.test(x, cdf))$statistic))
  # Draws from the prior itself: a small gap.  Draws far from it: near 1.
  expect_lt(spatialkit:::.gp_prior_gap(1 / stats::rgamma(2000, a, rate = b), cdf), 0.06)
  expect_gt(spatialkit:::.gp_prior_gap(stats::rnorm(500, 4, 0.1), cdf), 0.99)
  expect_true(is.na(spatialkit:::.gp_prior_gap(c(1, NA), cdf)))
  expect_identical(spatialkit:::.gp_share_floor, 0.01)
  expect_identical(spatialkit:::.gp_contraction_floor, 0.25)
  expect_identical(spatialkit:::.gp_prior_gap_floor, 0.30)
})

test_that("a length-scale posterior wider than its prior, on the long side, is not 'not identified'", {
  # Draws shaped like a real fit of a plane plus noise (n = 200, 100 km
  # square, default basis): length-scale mean/SD 1.48/0.44 and 2.06/0.52
  # against a prior with mean 0.91 and SD 0.29, GP share of the variance 0.71
  # to 0.95, R2 of the fitted surface 0.85.  The contraction is negative on
  # both axes, and on that alone the status was "not_identified".
  set.seed(1)
  N  <- 1000
  lg <- function(m, s) { v <- log(1 + s^2 / m^2); stats::rlnorm(N, log(m) - v / 2, sqrt(v)) }
  spec <- "inv_gamma(12.08121194, 10.05144272)"
  ls <- cbind(lx = lg(1.48, 0.44), ly = lg(2.06, 0.52))
  ad <- spatialkit:::.gp_adequacy(ls, sdgp_draws = lg(1.2, 0.3), sigma_draws = lg(0.5, 0.03),
                                  gp_ell_min = 0.46, prior_spec = spec,
                                  xy_scale = c(28000, 29000))
  expect_true(all(ad$contraction < 0))           # wider than the prior
  expect_true(all(ad$prior_gap > 0.5))           # and moved almost entirely
  expect_gt(ad$spatial_share[["q2.5"]], 0.5)
  expect_true(ad$identified)
  expect_identical(ad$identified_by, "spatial_share")
  expect_identical(ad$status, "ok")
  expect_equal(ad$frac_below, 0)

  # A family with no sigma has no share.  The same length-scales are then
  # identified because they moved, although they did not narrow.
  ns <- spatialkit:::.gp_adequacy(ls, sdgp_draws = lg(1.2, 0.3), gp_ell_min = 0.46,
                                  prior_spec = spec)
  expect_null(ns$spatial_share)
  expect_true(ns$identified)
  expect_identical(ns$identified_by, "length_scale")
  expect_identical(ns$status, "ok")
})

test_that("without a share the length-scales decide: narrowed or moved is identified, neither is not", {
  set.seed(2)
  N <- 2000
  a <- 12.08121194; b <- 10.05144272
  spec <- sprintf("inv_gamma(%.10g, %.10g)", a, b)
  prior_draws <- function() 1 / stats::rgamma(N, shape = a, rate = b)
  med <- 1 / stats::qgamma(0.5, shape = a, rate = b)
  psd <- b / ((a - 1) * sqrt(a - 2))
  sdgp <- stats::rnorm(N, 1, 0.05)

  # Equal to the prior on both axes: nothing narrowed, nothing moved.
  eq <- spatialkit:::.gp_adequacy(cbind(lx = prior_draws(), ly = prior_draws()),
                                  sdgp_draws = sdgp, gp_ell_min = 0.46, prior_spec = spec)
  expect_false(eq$identified)
  expect_identical(eq$identified_by, "length_scale")
  expect_identical(eq$status, "not_identified")
  expect_true(all(eq$contraction < 0.25) && all(eq$prior_gap < 0.30))

  # Narrowed around the prior's median (0.6 of its SD): the distribution
  # functions stay close, the variance does not.  A fit whose true
  # length-scale sat at the prior median looked like this (contraction 0.63
  # and 0.81, gaps 0.15 and 0.25).
  nar <- spatialkit:::.gp_adequacy(
    cbind(lx = stats::rnorm(N, med, 0.6 * psd), ly = prior_draws()),
    sdgp_draws = sdgp, gp_ell_min = 0.46, prior_spec = spec)
  expect_gt(nar$contraction[["lx"]], 0.25)
  expect_lt(nar$prior_gap[["lx"]], 0.30)
  expect_true(nar$identified)
  expect_identical(nar$status, "ok")

  # Moved to the long side on one axis and wider there: identified.
  mov <- spatialkit:::.gp_adequacy(
    cbind(lx = prior_draws(), ly = stats::rlnorm(N, log(1.9), 0.25)),
    sdgp_draws = sdgp, gp_ell_min = 0.46, prior_spec = spec)
  expect_lt(mov$contraction[["ly"]], 0)
  expect_gt(mov$prior_gap[["ly"]], 0.30)
  expect_true(mov$identified)
  expect_identical(mov$status, "ok")
})

test_that("where there is a share it decides, whatever the length-scales did", {
  set.seed(3)
  N <- 2000
  a <- 12.08121194; b <- 10.05144272
  spec <- sprintf("inv_gamma(%.10g, %.10g)", a, b)
  prior_draws <- function() 1 / stats::rgamma(N, shape = a, rate = b)
  # Length-scales the data did not move, under a GP that carries most of the
  # variance: identified.
  hi <- spatialkit:::.gp_adequacy(
    cbind(lx = prior_draws(), ly = prior_draws()),
    sdgp_draws = stats::rnorm(N, 1, 0.05), sigma_draws = stats::rnorm(N, 0.5, 0.02),
    gp_ell_min = 0.46, prior_spec = spec)
  expect_true(hi$identified)
  expect_identical(hi$status, "ok")
  # Length-scales far from their prior, under a GP that carries nothing: not
  # identified.
  lo <- spatialkit:::.gp_adequacy(
    cbind(lx = stats::rlnorm(N, log(1.9), 0.25), ly = stats::rlnorm(N, log(1.9), 0.25)),
    sdgp_draws = abs(stats::rnorm(N, 0, 0.03)), sigma_draws = rep(1.2, N),
    gp_ell_min = 0.46, prior_spec = spec)
  expect_false(lo$identified)
  expect_identical(lo$identified_by, "spatial_share")
  expect_identical(lo$status, "not_identified")
})

test_that(".lscale_prior_parse() reads the two priors the fitter sets", {
  p <- spatialkit:::.lscale_prior_parse("inv_gamma(11.48820413, 9.596223034)")
  # 1 - pgamma(1 / 1.7596, shape = 11.48820413, rate = 9.596223034)
  expect_equal(p$cdf(1.7596), 0.984, tolerance = 1e-3)
  expect_equal(p$var, 9.596223034^2 / (10.48820413^2 * 9.48820413))
  h <- spatialkit:::.lscale_prior_parse("normal(0, 2)")
  expect_equal(h$cdf(2), 2 * stats::pnorm(1) - 1)
  expect_equal(h$var, 4 * (1 - 2 / pi))
  expect_null(spatialkit:::.lscale_prior_parse("student_t(3, 0, 1)"))
  expect_null(spatialkit:::.lscale_prior_parse("inv_gamma(1, 2); inv_gamma(3, 4)"))
  expect_null(spatialkit:::.lscale_prior_parse(NULL))
})

test_that("a GP collapsed onto noise over short-range structure is reported", {
  skip_if_not_installed("brms")
  # An exponential field with an 8 km range; the default basis resolves 33 km.
  d   <- .fba_pts(range = 8000)
  des <- .fba_design(d)
  set.seed(2)
  N <- 800
  eng <- .fba_engine(
    ls = cbind(`lscale_gp..x..y..x` = 1 / stats::rgamma(N, des$shape, rate = des$scale),
               `lscale_gp..x..y..y` = 1 / stats::rgamma(N, des$shape, rate = des$scale)),
    sdgp = abs(stats::rnorm(N, 0, 0.05)), sigma = stats::rnorm(N, 1.05, 0.03))
  r  <- .fba_fit(d, eng)
  ad <- r$fit$info$gp_adequacy
  # The old check: about 1% of draws below the floor, so nothing was said.
  expect_lt(ad$frac_below, 0.10)
  expect_identical(ad$status, "unresolved_short_range")
  expect_lt(ad$residual_moran$p, 0.05)
  expect_identical(ad$residual_moran$k, 4L)
  expect_true("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  msg <- conditionMessage(r$warnings[[which(.fba_classes(r$warnings) ==
                                              "spatialkit_gp_adequacy")]])
  expect_match(msg, "the Gaussian process carries no variance that can be told from zero (its share of the variance has a 95% interval of 0.000 to",
               fixed = TRUE)
  expect_match(msg, "residuals are autocorrelated at short range", fixed = TRUE)
  expect_match(msg, "metre along x", fixed = TRUE)
  # It is a statement about the approximation, not the sampler.
  expect_true(isTRUE(r$fit$info$convergence_ok))
  expect_false("spatialkit_convergence" %in% .fba_classes(r$warnings))
  expect_identical(r$fit$info$convergence_diagnostics$gp_adequacy,
                   "unresolved_short_range")
  # print() and summary() both show it.
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "** GP adequacy: unresolved_short_range -- the GP's share of the variance is 0.000 to",
               fixed = TRUE)
  expect_false(grepl("Convergence warnings present", txt, fixed = TRUE))
  stx <- paste(utils::capture.output(print(summary(r$fit))), collapse = "\n")
  expect_match(stx, "** GP adequacy: unresolved_short_range", fixed = TRUE)
})

test_that("the same posterior over white noise is not called a failure", {
  skip_if_not_installed("brms")
  d   <- .fba_pts(seed = 21)
  des <- .fba_design(d)
  set.seed(3)
  N <- 800
  eng <- .fba_engine(
    ls = cbind(`lscale_gp..x..y..x` = 1 / stats::rgamma(N, des$shape, rate = des$scale),
               `lscale_gp..x..y..y` = 1 / stats::rgamma(N, des$shape, rate = des$scale)),
    sdgp = abs(stats::rnorm(N, 0, 0.05)), sigma = stats::rnorm(N, 1, 0.03))
  logged <- capture_spatialkit_log(r <- .fba_fit(d, eng), level = logger::WARN)
  ad <- r$fit$info$gp_adequacy
  expect_identical(ad$status, "no_signal_above_floor")
  expect_gt(ad$residual_moran$p, 0.05)
  # Logged and printed, not raised.
  expect_false("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  expect_true(log_has(logged, "carries no variance that can be told from zero"))
  expect_true(log_has(logged, "the residuals show no autocorrelation at short range"))
  expect_true(log_has(logged, "A response with no spatial structure gives this"))
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "GP adequacy: no_signal_above_floor -- the GP's share of the variance is 0.000 to",
               fixed = TRUE)
  expect_match(txt, "and the residuals show no autocorrelation at short range", fixed = TRUE)
  expect_false(grepl("** GP adequacy", txt, fixed = TRUE))
})

# The posterior of a real fit of a plane: both length-scales on the long side
# of the prior and wider than it, the GP carrying most of the variance.
.fba_plane_engine <- function(N = 800, seed = 6) {
  set.seed(seed)
  lg <- function(m, s) { v <- log(1 + s^2 / m^2); stats::rlnorm(N, log(m) - v / 2, sqrt(v)) }
  .fba_engine(
    ls = cbind(`lscale_gp..x..y..x` = lg(1.48, 0.44), `lscale_gp..x..y..y` = lg(2.06, 0.52)),
    sdgp = lg(1.2, 0.3), sigma = lg(0.5, 0.03))
}

test_that("a GP that carries most of the variance at a long range is 'ok' when its residuals show nothing left over", {
  skip_if_not_installed("brms")
  # The residuals are white noise (the stand-in's fitted values are zero).
  # The length-scale posterior is wider than its prior, which was read as
  # "not identified" and reported as "found no structure it can represent".
  logged <- capture_spatialkit_log(r <- .fba_fit(.fba_pts(seed = 21), .fba_plane_engine()),
                                   level = logger::WARN)
  ad <- r$fit$info$gp_adequacy
  expect_true(all(ad$contraction < 0))
  expect_gt(ad$spatial_share[["q2.5"]], 0.5)
  expect_true(ad$identified)
  expect_identical(ad$status, "ok")
  # "ok" now rests on the residuals as well, and the test is on the object.
  expect_named(ad$residual_moran, c("I", "z", "p", "k"))
  expect_gt(ad$residual_moran$p, 0.05)
  expect_identical(ad$residual_moran$k, 4L)
  expect_null(ad$residual_untested)
  expect_false("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  expect_false(log_has(logged, "carries no variance"))
  expect_false(log_has(logged, "could not be tested"))
  expect_identical(r$fit$info$convergence_diagnostics$gp_adequacy, "ok")
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "GP range: ", fixed = TRUE)
  expect_false(grepl("GP adequacy", txt, fixed = TRUE))
})

test_that("the same posterior over residuals with an 8 km range is not 'ok'", {
  skip_if_not_installed("brms")
  # An exponential field with an 8 km range under a 33 km floor, left entirely
  # in the residuals: what a regional surface beside short-range structure
  # leaves.  The posterior is the one above, so every figure the status used
  # to be read from says "ok": a share of the variance of 0.7 to 0.95, and
  # under 1 percent of the length-scale draws below the floor.  This returned
  # "ok" with no warning.
  d <- .fba_pts(range = 8000)
  r <- .fba_fit(d, .fba_plane_engine())
  ad <- r$fit$info$gp_adequacy
  expect_gt(ad$spatial_share[["q2.5"]], 0.5)
  expect_true(ad$identified)
  expect_lt(ad$frac_below, 0.01)
  expect_identical(ad$status, "unresolved_short_range")
  expect_lt(ad$residual_moran$p, 0.001)
  expect_gt(ad$residual_moran$I, 0.1)
  # The statistic is the one residual_morans_i() gives for these residuals.
  ref <- suppressWarnings(suppressMessages(
    residual_morans_i(r$fit, alternative = "greater", k = 4L)))
  expect_equal(ad$residual_moran$I, ref$observed)
  expect_equal(ad$residual_moran$p, ref$p_value)
  expect_identical(r$fit$info$convergence_diagnostics$gp_adequacy,
                   "unresolved_short_range")
  # Warned, in words for a Gaussian process that did find a surface.
  expect_true("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  msg <- conditionMessage(r$warnings[[which(.fba_classes(r$warnings) ==
                                              "spatialkit_gp_adequacy")]])
  expect_match(msg, sprintf("the Gaussian process carries %.2f to %.2f of the variance (95%% interval) at an effective range of ",
                            ad$spatial_share[["q2.5"]], ad$spatial_share[["q97.5"]]),
               fixed = TRUE)
  expect_match(msg, sprintf("and the residuals are still autocorrelated at short range (Moran's I = %.3f over 4 nearest neighbours, one-sided p = %.2g)",
                            ad$residual_moran$I, ad$residual_moran$p), fixed = TRUE)
  expect_match(msg, "The surface left structure behind", fixed = TRUE)
  expect_match(msg, "cannot be held by this basis (gp_k = ", fixed = TRUE)
  expect_match(msg, "A `gp_resolution` below", fixed = TRUE)
  expect_false(grepl("carries no variance", msg, fixed = TRUE))
  # About the approximation, not the sampler.
  expect_true(isTRUE(r$fit$info$convergence_ok))
  # print() and summary() flag it, with the share the GP does carry and the
  # test, and print() does not offer the fitted range as the range.
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, sprintf("** GP adequacy: unresolved_short_range -- the GP's share of the variance is %.3f to %.3f (95%% interval), and the residuals are autocorrelated at short range (Moran's I = %.3f, one-sided p = %.2g)",
                            ad$spatial_share[["q2.5"]], ad$spatial_share[["q97.5"]],
                            ad$residual_moran$I, ad$residual_moran$p), fixed = TRUE)
  expect_false(grepl("GP range: ", txt, fixed = TRUE))
  stx <- paste(utils::capture.output(print(summary(r$fit))), collapse = "\n")
  expect_match(stx, "** GP adequacy: unresolved_short_range", fixed = TRUE)
})

test_that("a variogram range below the floor stays flagged in print() unless the residuals answer it", {
  skip_if_not_installed("brms")
  # The 8 km field again, with the variogram before fitting taken to have
  # returned 9 km.  Both of its supporting tests agree, so the fit warns
  # before sampling.
  d <- .fba_pts(range = 8000)
  r <- .fba_fit(d, .fba_plane_engine(), sac = 9000)
  expect_identical(r$fit$info$gp_prefit$status, "below_floor")
  expect_true("spatialkit_gp_resolution" %in% .fba_classes(r$warnings))
  expect_identical(r$fit$info$gp_adequacy$status, "unresolved_short_range")
  txt <- utils::capture.output(print(r$fit))
  # This fit printed "the fitted model does not bear that out", unflagged.
  expect_true(any(grepl("** Range before fitting: a variogram put it at 9,000 metre, below what the basis resolves (an estimate; see GP adequacy)",
                        txt, fixed = TRUE)))
  expect_false(any(grepl("does not bear that out", txt, fixed = TRUE)))
  # The warning's advice carries the variogram's range, which the residuals
  # agree with.
  msg <- conditionMessage(r$warnings[[which(.fba_classes(r$warnings) ==
                                              "spatialkit_gp_adequacy")]])
  expect_match(msg, "The range a variogram gave before fitting (9,000 metre) is below that",
               fixed = TRUE)

  # The posterior draw fails, so there are no residuals to test.  The status
  # is "ok" on the posterior alone, which is said, and print() keeps the flag.
  logged <- capture_spatialkit_log(
    r2 <- .fba_fit(d, .fba_plane_engine(), sac = 9000,
                   epred = function(...) stop("no draws in a stand-in")),
    level = logger::WARN)
  ad2 <- r2$fit$info$gp_adequacy
  expect_identical(ad2$status, "ok")
  expect_null(ad2$residual_moran)
  expect_match(ad2$residual_untested, "no draws in a stand-in", fixed = TRUE)
  expect_false("spatialkit_gp_adequacy" %in% .fba_classes(r2$warnings))
  expect_true(log_has(logged, "the residuals could not be tested for autocorrelation at short range"))
  expect_true(log_has(logged, "rests on the posterior alone"))
  txt2 <- utils::capture.output(print(r2$fit))
  expect_true(any(grepl("** Range before fitting: a variogram put it at 9,000 metre, below what the basis resolves (an estimate; the residuals were not tested after fitting)",
                        txt2, fixed = TRUE)))
  expect_true(any(grepl("^  GP range: ", txt2)))
})

test_that(".gp_fit_resid_moran() returns the test, or the reason there is none", {
  skip_if_not_installed("brms")
  r <- .fba_fit(.fba_pts(seed = 21), .fba_plane_engine())
  local_mocked_bindings(
    posterior_epred = function(object, newdata, ...)
      matrix(0, nrow = 4L, ncol = nrow(newdata)),
    .package = "brms")
  got <- spatialkit:::.gp_fit_resid_moran(r$fit, "gaussian")
  expect_named(got, c("test", "reason"))
  expect_true(is.na(got$reason))
  expect_identical(got$test$k, 4L)
  ref <- suppressWarnings(suppressMessages(
    residual_morans_i(r$fit, alternative = "greater", k = 4L)))
  expect_equal(got$test[c("I", "z", "p")],
               list(I = ref$observed, z = ref$z, p = ref$p_value))
  # A family with a probability per category has no residual per row.
  cat_fam <- spatialkit:::.gp_fit_resid_moran(r$fit, "cumulative")
  expect_null(cat_fam$test)
  expect_match(cat_fam$reason, "no residual per row", fixed = TRUE)
  # A posterior that cannot be drawn: no test, the cause, and no warning or
  # log line of its own.
  clear_fitted_cache(r$fit)
  local_mocked_bindings(posterior_epred = function(...) stop("boom"), .package = "brms")
  logged <- capture_spatialkit_log(
    expect_no_warning(bad <- spatialkit:::.gp_fit_resid_moran(r$fit, "gaussian")),
    level = logger::WARN)
  expect_null(bad$test)
  expect_match(bad$reason, "residuals() failed: ", fixed = TRUE)
  expect_match(bad$reason, "boom", fixed = TRUE)
  expect_length(logged, 0L)
  expect_identical(spatialkit:::.gp_resid_k, 4L)
  expect_identical(spatialkit:::.gp_resid_level, 0.05)
})

test_that("an identified GP below the floor warns, and one clear of it prints its range", {
  skip_if_not_installed("brms")
  d   <- .fba_pts(range = 8000)
  des <- .fba_design(d)
  set.seed(4)
  N <- 800
  short <- .fba_engine(
    ls = cbind(`lscale_gp..x..y..x` = stats::rnorm(N, 0.9 * des$ell_min, 0.02),
               `lscale_gp..x..y..y` = stats::rnorm(N, 1.05 * des$ell_min, 0.02)),
    sdgp = stats::rnorm(N, 1, 0.05), sigma = stats::rnorm(N, 0.5, 0.02))
  r <- .fba_fit(d, short)
  expect_identical(r$fit$info$gp_adequacy$status, "below_resolution")
  expect_gt(r$fit$info$gp_adequacy$frac_below, 0.10)
  expect_equal(r$fit$info$convergence_diagnostics$gp_lscale_below_resolution,
               r$fit$info$gp_adequacy$frac_below)
  expect_true("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  msg <- conditionMessage(r$warnings[[which(.fba_classes(r$warnings) ==
                                              "spatialkit_gp_adequacy")]])
  expect_match(msg, "of the posterior length-scale draws fall below", fixed = TRUE)
  expect_match(msg, "A `gp_resolution` below", fixed = TRUE)
  # The residuals are read here too (the 8 km field, the stand-in's fitted
  # values being zero), and the warning adds what they show.
  expect_lt(r$fit$info$gp_adequacy$residual_moran$p, 0.05)
  expect_match(msg, "The residuals are autocorrelated at short range (Moran's I = ",
               fixed = TRUE)
  # This was a log line, and print() showed "GP basis: 23 per dim" and no more.
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "** GP adequacy: below_resolution", fixed = TRUE)
  expect_true(isTRUE(r$fit$info$convergence_ok))

  # Clear of the floor, over residuals with nothing in them (white noise:
  # over the 8 km field this posterior is "unresolved_short_range", as the
  # test above has it, and it used to be asserted "ok" here).
  d2   <- .fba_pts(seed = 21)
  des2 <- .fba_design(d2)
  long <- .fba_engine(
    ls = cbind(`lscale_gp..x..y..x` = stats::rnorm(N, 2 * des2$ell_min, 0.03),
               `lscale_gp..x..y..y` = stats::rnorm(N, 2.2 * des2$ell_min, 0.03)),
    sdgp = stats::rnorm(N, 1, 0.05), sigma = stats::rnorm(N, 0.5, 0.02))
  r2 <- .fba_fit(d2, long)
  ad <- r2$fit$info$gp_adequacy
  expect_identical(ad$status, "ok")
  expect_gt(ad$residual_moran$p, 0.05)
  expect_length(r2$warnings, 0L)
  # The posterior effective range, in map units, twice the floor on x.
  expect_equal(ad$effective_range["x", "q50"], 2 * r2$fit$info$gp_resolution[["x"]],
               tolerance = 0.01)
  txt2 <- paste(utils::capture.output(print(r2$fit)), collapse = "\n")
  expect_match(txt2, "GP range: ", fixed = TRUE)
  expect_false(grepl("GP adequacy", txt2, fixed = TRUE))
  # summary() adds nothing for a fit that is fine.  (fitted() was not needed
  # for this fit, so the stand-in for the posterior draw is set again here.)
  local_mocked_bindings(
    posterior_epred = function(object, newdata, ...)
      matrix(0, nrow = 4L, ncol = nrow(newdata)),
    .package = "brms")
  expect_false(grepl("GP ", paste(utils::capture.output(print(summary(r2$fit))),
                                  collapse = "\n"), fixed = TRUE))
})

test_that("an unchecked fit says so and an older fit prints as it did", {
  skip_if_not_installed("brms")
  d <- .fba_pts(n = 120)
  local_mocked_bindings(brm = function(...) structure(list(), class = "fba_stub"),
                        .package = "brms")
  local_mocked_bindings(estimate_sac_range = function(...) 9e4,
                        .package = "spatialkit")
  fit <- suppressWarnings(suppressMessages(fit_bayesian_spatial_model(
    d, "z", character(0), compute_loo = FALSE, check_convergence = FALSE)))
  expect_identical(fit$info$gp_adequacy$status, "not_checked")
  expect_length(fit$info$convergence_diagnostics, 0L)
  expect_false(grepl("GP adequacy", paste(utils::capture.output(print(fit)),
                                          collapse = "\n"), fixed = TRUE))
  # A fit saved before these fields existed prints without them.
  old <- fit
  old$info[c("gp_adequacy", "gp_resolution", "gp_prefit", "b_prior",
             "b_prior_source", "predictor_levels")] <- NULL
  txt <- paste(utils::capture.output(print(old)), collapse = "\n")
  expect_match(txt, "GP basis: ", fixed = TRUE)
  expect_false(grepl("GP resolves", txt, fixed = TRUE))
  expect_no_error(utils::capture.output(print(old)))
})


# ---------------------------------------------------------------------------
# Real sampler (opt-in).  Six default fits of two chains, three of 600
# iterations and three of 1,000, one to two minutes each.
# ---------------------------------------------------------------------------

.fba_stan_fit <- function(pts, iter = 600, seed = 123) {
  w <- list(); fit <- NULL
  utils::capture.output(
    fit <- withCallingHandlers(
      suppressMessages(fit_bayesian_spatial_model(
        pts, "z", character(0), chains = 2, iter = iter, cores = 1, seed = seed,
        compute_loo = FALSE, backend = "rstan")),
      warning = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning") }),
    type = "output")
  list(fit = fit, warnings = w)
}

test_that("a field the basis resolves is 'ok' and its range is recovered in map units", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  # A squared-exponential field on a 100 km square whose length-scale is the
  # median of the prior this layout gets: an effective range near 0.6 of the
  # side, well above the 32 km floor.
  n <- 400
  set.seed(1)
  xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  des <- .fba_design(sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
                                  coords = c("x", "y"), crs = 32617))
  med <- 1 / stats::qgamma(0.5, des$shape, rate = des$scale)
  ell <- med * mean(c(stats::sd(xy[, 1]), stats::sd(xy[, 2])))
  e <- eigen(exp(-as.matrix(stats::dist(xy))^2 / (2 * ell^2)), symmetric = TRUE)
  set.seed(5)
  f <- as.numeric(e$vectors %*% (sqrt(pmax(e$values, 0)) * stats::rnorm(n)))
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2],
                                 z = f + stats::rnorm(n, 0, sqrt(0.3))),
                      coords = c("x", "y"), crs = 32617)
  r  <- .fba_stan_fit(pts)
  ad <- r$fit$info$gp_adequacy
  expect_identical(ad$status, "ok")
  expect_true(ad$identified)
  expect_lt(ad$frac_below, 0.10)
  # And nothing is left in the residuals at short range.
  expect_gt(ad$residual_moran$p, 0.05)
  expect_false("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  # The posterior effective range, in metres, brackets the truth on each axis.
  truth <- sqrt(2 * log(20)) * ell
  expect_true(all(ad$effective_range[, "q2.5"] < truth))
  expect_true(all(ad$effective_range[, "q97.5"] > truth))
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "GP range: ", fixed = TRUE)
})

test_that("white noise is 'no_signal_above_floor', logged and not raised", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  n <- 400
  set.seed(1)
  xy <- cbind(runif(n, 0, 8e5), runif(n, 0, 3e5))
  set.seed(77)
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2],
                                 z = stats::rnorm(n, 0, sqrt(1.3))),
                      coords = c("x", "y"), crs = 32617)
  r  <- .fba_stan_fit(pts)
  ad <- r$fit$info$gp_adequacy
  # The share of draws below the floor (0.0075 in the fit this reproduces)
  # said nothing; the status does.
  expect_lt(ad$frac_below, 0.10)
  expect_false(ad$identified)
  expect_identical(ad$status, "no_signal_above_floor")
  expect_false("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "GP adequacy: no_signal_above_floor", fixed = TRUE)
})

test_that("a plane plus noise is 'ok': the GP carries the variance at a range longer than the extent", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  # The fit that was reported as "no spatial structure found": R2 of the
  # fitted surface 0.85, the GP's share of the variance 0.71 to 0.95, and a
  # length-scale posterior on the long side of its prior and wider than it
  # (contraction -1.35 and -2.27).
  n <- 200
  set.seed(42)
  xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  z  <- 3 * xy[, 1] / 1e5 + 2 * xy[, 2] / 1e5 + stats::rnorm(n, 0, 0.5)
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
                      coords = c("x", "y"), crs = 32617)
  r  <- .fba_stan_fit(pts)
  ad <- r$fit$info$gp_adequacy
  expect_true(all(ad$contraction < 0))
  expect_true(all(ad$prior_gap > 0.30))
  expect_gt(ad$spatial_share[["q2.5"]], 0.3)
  expect_true(ad$identified)
  expect_identical(ad$status, "ok")
  expect_lt(ad$frac_below, 0.10)
  # The residuals of a plane the surface took whole show nothing at short
  # range (p was 0.90 in the fit this reproduces).
  expect_gt(ad$residual_moran$p, 0.05)
  expect_false("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  # The posterior effective range is longer than the 32 km floor on both axes.
  expect_true(all(ad$effective_range[, "q2.5"] > max(r$fit$info$gp_resolution)))
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "GP range: ", fixed = TRUE)
  expect_false(grepl("no_signal_above_floor", txt, fixed = TRUE))
  # The surface it fitted explains most of the response.
  expect_gt(1 - stats::var(stats::residuals(r$fit)) / stats::var(z), 0.6)
})

test_that("a field with a 200 km range and a bowl stay 'ok': nothing is left in their residuals", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  # Two more surfaces that are smooth over the whole extent, with the plane
  # above the three fits whose length-scale posterior is wider than its
  # prior.  Reading the residuals of every fit must not cost them their
  # status: p was 0.90 and 0.97 in the fits these reproduce.
  n <- 200
  set.seed(42)
  xy  <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  ell <- 2e5 / sqrt(2 * log(20))
  e   <- eigen(exp(-as.matrix(stats::dist(xy))^2 / (2 * ell^2)), symmetric = TRUE)
  z_se200 <- 2 * as.numeric(e$vectors %*% (sqrt(pmax(e$values, 0)) * stats::rnorm(n))) +
    stats::rnorm(n, 0, 0.5)
  set.seed(42)
  xy_b   <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  z_bowl <- 4 * ((xy_b[, 1] / 1e5 - 0.5)^2 + (xy_b[, 2] / 1e5 - 0.5)^2) / 0.5 +
    stats::rnorm(n, 0, 0.3)
  for (z in list(z_se200, z_bowl)) {
    pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
                        coords = c("x", "y"), crs = 32617)
    r  <- .fba_stan_fit(pts, iter = 1000, seed = 3)
    ad <- r$fit$info$gp_adequacy
    expect_true(ad$identified)
    expect_identical(ad$status, "ok")
    expect_gt(ad$residual_moran$p, 0.05)
    expect_lt(ad$residual_moran$I, 0.05)
    expect_false("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
    txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
    expect_match(txt, "GP range: ", fixed = TRUE)
    expect_false(grepl("GP adequacy", txt, fixed = TRUE))
  }
})

test_that("a weak plane plus a field with an 8 km range is not 'ok': the residuals hold what the surface left", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  # The fit that came back "ok" with no warning of any class and
  # "Convergence: passed": a plane rising 1 along x and 0.7 along y, an
  # exponential field with an 8 km effective range (about four fifths of the
  # variance) and noise, on a 100 km square whose basis resolves 33 km.  The
  # Gaussian process took the plane (a share of the variance of 0.09 to 0.54
  # at 55 and 72 km, 2 percent of length-scale draws below the floor), its
  # surface explained 19 percent of the response, and the residuals had a
  # Moran's I of 0.143 (p = 5.6e-06) that nothing read.
  n <- 400
  set.seed(1001)
  xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  e  <- eigen(exp(-as.matrix(stats::dist(xy)) / (8000 / 3)), symmetric = TRUE)
  f  <- as.numeric(e$vectors %*% (sqrt(pmax(e$values, 0)) * stats::rnorm(n)))
  z  <- xy[, 1] / 1e5 + 0.7 * xy[, 2] / 1e5 + f + stats::rnorm(n, 0, sqrt(0.1))
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
                      coords = c("x", "y"), crs = 32617)
  r  <- .fba_stan_fit(pts, iter = 1000, seed = 3)
  ad <- r$fit$info$gp_adequacy
  # Everything the posterior says is what an adequate fit says.
  expect_true(ad$identified)
  expect_gt(ad$spatial_share[["q2.5"]], 0.01)
  expect_lt(ad$frac_below, 0.10)
  # The residuals do not.
  expect_false(identical(ad$status, "ok"))
  expect_identical(ad$status, "unresolved_short_range")
  expect_lt(ad$residual_moran$p, 0.001)
  expect_gt(ad$residual_moran$I, 0.05)
  expect_true("spatialkit_gp_adequacy" %in% .fba_classes(r$warnings))
  msg <- conditionMessage(r$warnings[[which(.fba_classes(r$warnings) ==
                                              "spatialkit_gp_adequacy")]])
  expect_match(msg, "and the residuals are still autocorrelated at short range",
               fixed = TRUE)
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "** GP adequacy: unresolved_short_range", fixed = TRUE)
  expect_false(grepl("GP range: ", txt, fixed = TRUE))
  # Most of the response is still in the residuals.
  expect_lt(1 - stats::var(stats::residuals(r$fit)) / stats::var(z), 0.5)
})
