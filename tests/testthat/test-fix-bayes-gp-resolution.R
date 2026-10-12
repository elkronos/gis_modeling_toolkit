# What the GP basis of fit_bayesian_spatial_model() resolves, in map units,
# and what is said about it before anything is compiled.
#
# The basis has a floor: an effective range below which it cannot represent
# structure.  With the defaults that floor is 14 to 33 percent of the extent
# along each axis, and nothing reported it in map units, compared it with the
# data, or let a caller ask for a finer one.
#
# No Stan: brms::brm() is replaced by a recorder.

.fbr_capture <- function(...) {
  cap <- new.env()
  local_mocked_bindings(
    brm = function(...) { cap$args <- list(...); structure(list(), class = "fbr_stub") },
    .package = "brms")
  w <- list()
  fit <- withCallingHandlers(
    suppressMessages(fit_bayesian_spatial_model(
      ..., compute_loo = FALSE, check_convergence = FALSE)),
    warning = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning") })
  list(fit = fit, args = cap$args, warnings = w)
}

.fbr_classes <- function(w) {
  vapply(w, function(c) class(c)[1L], character(1))
}

# 500 points uniform on a 100 km square; the response is white noise unless
# `range` gives an exponential field's effective range.
.fbr_pts <- function(n = 500, seed = 20261011, side = 1e5, range = NULL,
                     nugget = 0.09) {
  set.seed(seed)
  xy <- cbind(runif(n, 0, side), runif(n, 0, side))
  z <- if (is.null(range)) stats::rnorm(n) else {
    D <- as.matrix(stats::dist(xy))
    as.numeric(t(chol(exp(-D / (range / 3)) + diag(1e-8, n))) %*% stats::rnorm(n)) +
      stats::rnorm(n, 0, sqrt(nugget))
  }
  sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z, a = stats::rnorm(n)),
               coords = c("x", "y"), crs = 32617)
}

# A stand-in for estimate_sac_range() that returns what it is told.
.fbr_sac <- function(value, reason = NULL) {
  function(...) if (is.null(reason)) value else
    structure(NA_real_, rejected_reason = reason)
}


test_that("the fit stores the basis floor as an effective range in map units", {
  skip_if_not_installed("brms")
  d <- .fbr_pts()
  local_mocked_bindings(estimate_sac_range = .fbr_sac(9e4), .package = "spatialkit")
  r <- .fbr_capture(d, "z", character(0))
  inf <- r$fit$info
  expect_named(inf$gp_resolution, c("x", "y"))
  # By definition: sqrt(2 log 20) length-scales, each axis on its own SD.
  cs <- inf$coord_scaling
  expect_equal(unname(inf$gp_resolution),
               sqrt(2 * log(20)) * inf$gp_ell_min * c(cs$x_scale, cs$y_scale))
  # Recomputed from the rule itself, outside the package: the lower bound is
  # the 25th percentile of the pairwise distances as a length-scale, the
  # boundary factor follows the upper bound, and k resolves the lower one.
  xy <- sf::st_coordinates(d)
  z  <- cbind((xy[, 1] - mean(xy[, 1])) / stats::sd(xy[, 1]),
              (xy[, 2] - mean(xy[, 2])) / stats::sd(xy[, 2]))
  dd <- as.numeric(stats::dist(z)); dd <- dd[dd > 0]
  Fse <- sqrt(2 * log(20))
  lower <- unname(stats::quantile(dd, 0.25)) / Fse
  upper <- max(dd) / Fse
  Xc <- sweep(z, 2, colMeans(z)); S <- max(1, max(Xc) - min(Xc))
  cc <- max(3.2 * upper / S, 1.25)
  k  <- min(max(ceiling(1.75 * cc / (lower / S)), 10), 50)
  expect_identical(inf$gp_k, as.integer(k))
  expect_equal(unname(inf$gp_resolution),
               Fse * (1.75 * cc * S / k) * c(stats::sd(xy[, 1]), stats::sd(xy[, 2])))
  # The value the review recorded for this layout: 32.3 km and 32.5 km, a
  # third of the side.
  expect_equal(unname(inf$gp_resolution) / 1000, c(32.3, 32.5), tolerance = 5e-3)
  expect_null(inf$gp_resolution_requested)

  # print() gives it, with the CRS's unit.
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "GP resolves: effective ranges from 32,300 (x) and 32,500 (y) metre",
               fixed = TRUE)
})

test_that("a response whose range is below the floor is warned about before compiling", {
  skip_if_not_installed("brms")
  skip_if_not_installed("gstat")
  # An exponential field with a 6 km effective range on a 100 km square whose
  # default basis resolves 33 km.  Seeded: the range estimate is 7.3 km.
  d <- .fbr_pts(n = 400, seed = 5, range = 6000)
  r <- .fbr_capture(d, "z", character(0))
  expect_true("spatialkit_gp_resolution" %in% .fbr_classes(r$warnings))
  w   <- r$warnings[[which(.fbr_classes(r$warnings) == "spatialkit_gp_resolution")]]
  msg <- conditionMessage(w)
  pf  <- r$fit$info$gp_prefit
  expect_identical(pf$status, "below_floor")
  expect_lt(pf$range, pf$floor)
  expect_equal(pf$floor, max(r$fit$info$gp_resolution))
  # Both numbers, in map units, and the basis size that would be needed.
  fmt <- function(v) format(signif(v, 3), big.mark = ",", scientific = FALSE, trim = TRUE)
  expect_match(msg, sprintf("at %s metre, shorter than the %s metre", fmt(pf$range),
                            fmt(pf$floor)), fixed = TRUE)
  expect_match(msg, sprintf("gp_k = %s per dimension", format(pf$k_needed)), fixed = TRUE)
  # It is said as what it is: a variogram estimate made before fitting, which
  # the check after sampling tests on the residuals.  It used to state the
  # range as fact and say "this model cannot represent it"; and then to say
  # the later check "supersedes" it, when that check read the residuals only
  # for a Gaussian process that had found nothing.
  expect_match(msg, "a variogram fitted before sampling puts the autocorrelation range of the response at",
               fixed = TRUE)
  expect_match(msg, "This is an estimate and can be wrong", fixed = TRUE)
  expect_match(msg, "the check after sampling ($info$gp_adequacy, shown by print()) reads the fitted model and tests its residuals for what the surface left behind",
               fixed = TRUE)
  expect_false(grepl("supersedes", msg, fixed = TRUE))
  expect_match(msg, "If the range is right, structure at that range is counted as noise",
               fixed = TRUE)
  expect_false(grepl("cannot represent it", msg, fixed = TRUE))
  # Both tests that have to agree with the range do: autocorrelated at short
  # range, and no longer at the floor's distance.
  expect_lt(pf$moran$p, 0.05)
  expect_false(pf$moran_floor$I > 0.05 && pf$moran_floor$p < 0.05)
  expect_equal(pf$moran_floor$band, c(0.75, 1.25) * pf$floor)
  expect_match(msg, "the autocorrelation at the floor's distance agrees with it (Moran's I = ",
               fixed = TRUE)
  # 104 per dimension is past the cap of 50: the message says so and points
  # elsewhere, on the condition that the fitted model bears the range out.
  expect_gt(pf$k_needed, 50)
  expect_match(msg, "above the cap of 50", fixed = TRUE)
  expect_match(msg, "if the fitted model bears the range out, use kriging", fixed = TRUE)
  # The warning came before the model was built.  Nothing was checked after
  # it here, so print() flags the line.
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, sprintf("** Range before fitting: a variogram put it at %s metre",
                            fmt(pf$range)), fixed = TRUE)

  # A range the cap can reach is answered with the argument to pass.
  local_mocked_bindings(estimate_sac_range = .fbr_sac(2.4e4), .package = "spatialkit")
  r2 <- .fbr_capture(d, "z", character(0))
  m2 <- conditionMessage(r2$warnings[[which(.fbr_classes(r2$warnings) ==
                                              "spatialkit_gp_resolution")]])
  expect_match(m2, "gp_resolution = 24000 fits that basis", fixed = TRUE)
})

test_that("a short variogram range on a response still autocorrelated at the floor's distance is not warned about", {
  skip_if_not_installed("brms")
  # Smooth over long distances: a plane, and a sine with a 100 km wavelength,
  # each plus noise.  A variogram returned ranges of 5 to 7 km for surfaces
  # like these (true ranges 60 to 90 km) in 1 to 9 draws of 60, and the
  # short-lag Moran's I is significant for any real structure, so the range
  # was warned about as fact.  Here the variogram is handed that range.
  n <- 400
  set.seed(31)
  xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  mk <- function(z) sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z),
                                 coords = c("x", "y"), crs = 32617)
  local_mocked_bindings(estimate_sac_range = .fbr_sac(6000), .package = "spatialkit")
  for (z in list(3 * xy[, 1] / 1e5 + 2 * xy[, 2] / 1e5 + stats::rnorm(n, 0, 0.5),
                 sin(2 * pi * xy[, 1] / 1e5) + stats::rnorm(n))) {
    logged <- capture_spatialkit_log(r <- .fbr_capture(mk(z), "z", character(0)),
                                     level = logger::WARN)
    pf <- r$fit$info$gp_prefit
    expect_identical(pf$status, "below_floor_long_range")
    expect_lt(pf$moran$p, 0.05)                  # the gate there was
    expect_gt(pf$moran_floor$I, 0.05)            # the one that was missing
    expect_lt(pf$moran_floor$p, 0.05)
    expect_false("spatialkit_gp_resolution" %in% .fbr_classes(r$warnings))
    # Logged, with both numbers, and printed without a flag.
    expect_true(log_has(logged, "puts the autocorrelation range of the response at 6,000 metre"))
    expect_true(log_has(logged, "is still autocorrelated at that distance"))
    expect_true(log_has(logged, "no warning is raised on it"))
    txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
    expect_match(txt, "  Range before fitting: a variogram put it at 6,000 metre, but the response is still autocorrelated at the floor's distance",
                 fixed = TRUE)
    expect_false(grepl("** Range before fitting", txt, fixed = TRUE))
  }
})

test_that(".gp_floor_moran() is Moran's I over the pairs about one floor apart", {
  n <- 300
  set.seed(8)
  xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  fl <- 3.2e4
  # A plane: strongly correlated at any distance short of the extent.
  pl <- spatialkit:::.gp_floor_moran(xy, xy[, 1] / 1e5 + stats::rnorm(n, 0, 0.1), floor = fl)
  expect_gt(pl$I, 0.2)
  expect_lt(pl$p, 1e-6)
  expect_equal(pl$band, c(0.75, 1.25) * fl)
  # Recomputed by hand: binary weights over the band, residuals of the mean.
  z <- xy[, 1] / 1e5 + stats::rnorm(n, 0, 0.1)
  got <- spatialkit:::.gp_floor_moran(xy, z, floor = fl)
  D <- as.matrix(stats::dist(xy)); W <- (D >= 0.75 * fl & D < 1.25 * fl) * 1
  e <- z - mean(z)
  expect_equal(got$I, (n / sum(W)) * sum(e * (W %*% e)) / sum(e^2))
  expect_equal(got$n_pairs, sum(W) / 2)
  # White noise: nothing there.
  wn <- spatialkit:::.gp_floor_moran(xy, stats::rnorm(n), floor = fl)
  expect_lt(abs(wn$I), 0.05)
  # With predictors the residuals are the ones tested: a response that is the
  # predictor has none.
  a <- xy[, 1] / 1e5
  ap <- spatialkit:::.gp_floor_moran(xy, 2 * a + stats::rnorm(n, 0, 0.1), X = cbind(a = a),
                                     floor = fl)
  expect_lt(abs(ap$I), 0.05)
  # More rows than max_n are thinned at even steps, with no random draw.
  set.seed(1); before <- .Random.seed
  th <- spatialkit:::.gp_floor_moran(xy, z, floor = fl, max_n = 100L)
  expect_identical(.Random.seed, before)
  expect_lt(th$n_pairs, got$n_pairs / 4)
  expect_gt(th$I, 0.2)
  # No pair in the band, a floor that is not a distance, too few rows: NULL.
  expect_null(spatialkit:::.gp_floor_moran(xy, z, floor = 1e7))
  expect_null(spatialkit:::.gp_floor_moran(xy, z, floor = NA_real_))
  expect_null(spatialkit:::.gp_floor_moran(xy[1:5, ], z[1:5], floor = fl))
  expect_identical(spatialkit:::.gp_range_correlation, 0.05)
})

test_that("print() flags the pre-fit range unless the residuals after sampling answer it", {
  lines <- spatialkit:::.bayes_gp_print_lines
  eff <- matrix(c(5e4, 6e4, 7e4, 5e4, 6e4, 7e4), nrow = 2, byrow = TRUE,
                dimnames = list(c("x", "y"), c("q2.5", "q50", "q97.5")))
  info <- function(status, moran = NULL, prefit = "below_floor") list(
    gp_resolution = c(x = 32300, y = 32500),
    gp_prefit = list(status = prefit, range = 6000, unit = "metre"),
    gp_adequacy = list(status = status, unit = "metre", effective_range = eff,
                       frac_below = 0.4, spatial_share = c(0, 0.01, 0.1),
                       residual_moran = moran))
  clear <- list(I = -0.012, z = -0.3, p = 0.62, k = 4L)
  left  <- list(I = 0.143, z = 4.4, p = 5.6e-06, k = 4L)
  # "ok" on a residual test that found nothing: the variogram's range stays
  # on the page as the estimate it was, unflagged, beside the test that
  # answered it.  It used to be flagged whatever the fit found, and then to
  # read "the fitted model does not bear that out" whenever the posterior
  # alone looked adequate, which was printed over a field with a 10 km range.
  ok <- lines(info("ok", clear))
  expect_true(any(ok == "  Range before fitting: a variogram put it at 6,000 metre, below what the basis resolves; after fitting, the residuals show no autocorrelation at short range (Moran's I = -0.012, one-sided p = 0.62)"))
  expect_false(any(grepl("** Range before fitting", ok, fixed = TRUE)))
  expect_false(any(grepl("does not bear that out", ok, fixed = TRUE)))
  expect_true(any(grepl("^  GP range: ", ok)))
  # "ok" with no residual test behind it keeps the flag and says why.
  nt <- lines(info("ok"))
  expect_true(any(nt == "  ** Range before fitting: a variogram put it at 6,000 metre, below what the basis resolves (an estimate; the residuals were not tested after fitting)"))
  expect_true(any(grepl("^  GP range: ", nt)))
  # Not checked, or found short: flagged.
  for (st in c("not_checked", "below_resolution", "unresolved_short_range",
               "no_signal_above_floor")) {
    got <- lines(info(st, if (st == "unresolved_short_range") left))
    expect_true(any(grepl("** Range before fitting: a variogram put it at 6,000 metre, below what the basis resolves (an estimate; see GP adequacy)",
                          got, fixed = TRUE)), info = st)
  }
  # The adequacy line carries the test it rests on.
  us <- lines(info("unresolved_short_range", left))
  expect_true(any(us == "  ** GP adequacy: unresolved_short_range -- the GP's share of the variance is 0.000 to 0.100 (95% interval), and the residuals are autocorrelated at short range (Moran's I = 0.143, one-sided p = 5.6e-06)"))
  expect_false(any(grepl("^  GP range: ", us)))

  # A short range beside structure that reaches the floor is not read as the
  # range, and is not a flag, until the residuals show the short range was
  # there as well.
  lr_ok <- lines(info("ok", clear, prefit = "below_floor_long_range"))
  expect_true(any(lr_ok == "  Range before fitting: a variogram put it at 6,000 metre, but the response is still autocorrelated at the floor's distance, so that is not read as its range"))
  for (st in c("unresolved_short_range", "below_resolution")) {
    lr_us <- lines(info(st, left, prefit = "below_floor_long_range"))
    expect_true(any(lr_us == "  ** Range before fitting: a variogram put it at 6,000 metre, below what the basis resolves; the response is also autocorrelated at the floor's distance (see GP adequacy)"),
                info = st)
  }
  # A fit with no adequacy record (saved by an earlier version, or a stand-in)
  # prints the unflagged line and does not stop.
  bare <- list(gp_resolution = c(x = 32300, y = 32500),
               gp_prefit = list(status = "below_floor_long_range", range = 6000,
                                unit = "metre"))
  expect_true(any(grepl("^  Range before fitting: a variogram put it at 6,000 metre, but the response is still autocorrelated",
                        lines(bare))))
})

test_that("a range that cannot be estimated is said, not passed over", {
  skip_if_not_installed("brms")
  d <- .fbr_pts(n = 400, seed = 5, range = 6000)
  local_mocked_bindings(
    estimate_sac_range = .fbr_sac(NA, reason = "variogram model did not converge"),
    .package = "spatialkit")
  logged <- capture_spatialkit_log(r <- .fbr_capture(d, "z", character(0)),
                                   level = logger::WARN)
  pf <- r$fit$info$gp_prefit
  expect_identical(pf$status, "not_estimated")
  expect_identical(pf$reason, "variogram model did not converge")
  expect_true(is.na(pf$range))
  # Not an R warning (the range is unknown, not short) ...
  expect_false("spatialkit_gp_resolution" %in% .fbr_classes(r$warnings))
  # ... but logged with the reason, and shown by print().
  expect_true(log_has(logged, "could not be estimated before fitting \\(variogram model did not converge\\)"))
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "Range before fitting: not estimated (variogram model did not converge)",
               fixed = TRUE)
  # This field's range is 6 km: autocorrelated at short range and no longer
  # at the floor's distance, so a range below the floor is not ruled out, and
  # the log says so with both statistics.
  expect_lt(pf$moran$p, 0.05)
  expect_false(pf$moran_floor$I > 0.05 && pf$moran_floor$p < 0.05)
  expect_true(log_has(logged, "It is autocorrelated at short range \\(Moran's I = "))
  expect_true(log_has(logged, "and not at the floor's distance"))
  expect_true(log_has(logged, "so a range below that floor is not ruled out"))

  # A plane is autocorrelated at short range too, as any real structure is.
  # It used to get the same "not ruled out"; its structure reaches past the
  # floor, and the log now says that.
  n <- 400
  set.seed(31)
  xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5))
  pl <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2],
                                z = 3 * xy[, 1] / 1e5 + 2 * xy[, 2] / 1e5 +
                                  stats::rnorm(n, 0, 0.5)),
                     coords = c("x", "y"), crs = 32617)
  logged2 <- capture_spatialkit_log(r2 <- .fbr_capture(pl, "z", character(0)),
                                    level = logger::WARN)
  pf2 <- r2$fit$info$gp_prefit
  expect_identical(pf2$status, "not_estimated")
  expect_lt(pf2$moran$p, 0.05)
  expect_gt(pf2$moran_floor$I, 0.05)
  expect_true(log_has(logged2, "It is still autocorrelated at the floor's distance"))
  expect_true(log_has(logged2, "so its structure reaches at least as far as the basis resolves"))
  expect_false(log_has(logged2, "not ruled out"))
})

test_that("a short variogram range on a response with no autocorrelation is not read as structure", {
  skip_if_not_installed("brms")
  # White noise.  A variogram fitted to it returns a finite range below the
  # floor in about one draw in three; here it is handed one.
  d <- .fbr_pts(n = 400, seed = 12)
  local_mocked_bindings(estimate_sac_range = .fbr_sac(7000), .package = "spatialkit")
  r <- .fbr_capture(d, "z", character(0))
  expect_identical(r$fit$info$gp_prefit$status, "below_floor_no_autocorrelation")
  expect_gt(r$fit$info$gp_prefit$moran$p, 0.05)
  expect_false("spatialkit_gp_resolution" %in% .fbr_classes(r$warnings))
  # A range above the floor is simply fine.
  local_mocked_bindings(estimate_sac_range = .fbr_sac(6e4), .package = "spatialkit")
  r2 <- .fbr_capture(d, "z", character(0))
  expect_identical(r2$fit$info$gp_prefit$status, "ok")
  expect_length(r2$warnings, 0L)
})

test_that("gp_resolution derives a finer basis and moves the prior with it", {
  skip_if_not_installed("brms")
  d <- .fbr_pts()
  local_mocked_bindings(estimate_sac_range = .fbr_sac(9e4), .package = "spatialkit")
  def <- .fbr_capture(d, "z", character(0))$fit$info
  fin <- .fbr_capture(d, "z", character(0), gp_resolution = 20000)$fit$info

  expect_identical(fin$gp_resolution_requested, 20000)
  expect_gt(fin$gp_k, def$gp_k)
  # The basis resolves what was asked for, on both axes.
  expect_true(all(fin$gp_resolution <= 20000))
  # And not much more: one basis function fewer would not.
  cs <- fin$coord_scaling
  one_fewer <- sqrt(2 * log(20)) * (1.75 * fin$gp_c * fin$gp_S / (fin$gp_k - 1)) *
    max(cs$x_scale, cs$y_scale)
  expect_gt(one_fewer, 20000)
  # The prior's lower 1% tail sits at the length-scale that basis resolves.
  expect_equal(fin$gp_lengthscale_bounds[["lower"]], fin$gp_ell_min)
  expect_lt(fin$gp_lengthscale_bounds[["lower"]], def$gp_lengthscale_bounds[["lower"]])
  expect_false(identical(fin$gp_lscale_prior, def$gp_lscale_prior))
  ig <- spatialkit:::.lscale_invgamma(fin$gp_ell_min, fin$gp_lengthscale_bounds[["upper"]])
  expect_identical(fin$gp_lscale_prior,
                   sprintf("inv_gamma(%.10g, %.10g)", ig$shape, ig$scale))
  expect_lt(fin$gp_adequacy$prior_mass_below, 0.02)

  # A request the default basis already meets changes nothing.
  same <- .fbr_capture(d, "z", character(0), gp_resolution = 60000)$fit$info
  expect_identical(same$gp_k, def$gp_k)
  expect_identical(same$gp_lscale_prior, def$gp_lscale_prior)
})

test_that("gp_resolution past the cap, or beside gp_k, is an error that says what is possible", {
  skip_if_not_installed("brms")
  d <- .fbr_pts()
  calls <- 0L
  local_mocked_bindings(
    brm = function(...) { calls <<- calls + 1L; structure(list(), class = "fbr_stub") },
    .package = "brms")
  # 5 km on a 100 km square needs about 150 functions per dimension.
  err <- tryCatch(
    suppressMessages(fit_bayesian_spatial_model(d, "z", character(0),
                                                gp_resolution = 5000)),
    error = conditionMessage)
  expect_match(err, "`gp_resolution` = 5,000 metre needs gp_k = ", fixed = TRUE)
  expect_match(err, "above the cap of 50 (2500)", fixed = TRUE)
  expect_match(err, "At the cap this layout resolves effective ranges from", fixed = TRUE)
  expect_match(err, "kriging", fixed = TRUE)
  expect_error(suppressMessages(fit_bayesian_spatial_model(
    d, "z", character(0), gp_resolution = 20000, gp_k = 30)),
    "pass `gp_resolution` or `gp_k`, not both")
  for (bad in list(-1, 0, c(1, 2), NA_real_, "20 km"))
    expect_error(suppressMessages(fit_bayesian_spatial_model(
      d, "z", character(0), gp_resolution = bad)),
      "`gp_resolution` must be a single positive number")
  # Nothing was compiled for any of them.
  expect_identical(calls, 0L)
})

test_that("an explicit gp_k finer than the derived one moves the prior's lower tail; the default does not", {
  skip_if_not_installed("brms")
  d <- .fbr_pts()
  local_mocked_bindings(estimate_sac_range = .fbr_sac(9e4), .package = "spatialkit")
  def <- .fbr_capture(d, "z", character(0))$fit$info
  k50 <- .fbr_capture(d, "z", character(0), gp_k = 50)$fit$info
  # The percentile bounds, as the default prior has always been calibrated.
  xy <- sf::st_coordinates(d)
  zz <- cbind((xy[, 1] - mean(xy[, 1])) / stats::sd(xy[, 1]),
              (xy[, 2] - mean(xy[, 2])) / stats::sd(xy[, 2]))
  b  <- gp_lengthscale_bounds(zz)
  ig <- spatialkit:::.lscale_invgamma(b[["lower"]], b[["upper"]])
  expect_identical(def$gp_lscale_prior,
                   sprintf("inv_gamma(%.10g, %.10g)", ig$shape, ig$scale))
  expect_equal(def$gp_lengthscale_bounds, b)
  # At gp_k = 50 the prior kept those bounds, so the finer basis was held to
  # the coarser one's length-scales.  Now its lower tail is what k = 50
  # resolves.
  expect_equal(k50$gp_lengthscale_bounds[["lower"]], k50$gp_ell_min)
  expect_lt(k50$gp_ell_min, b[["lower"]])
  expect_equal(k50$gp_lengthscale_bounds[["upper"]], b[["upper"]])
  expect_false(identical(k50$gp_lscale_prior, def$gp_lscale_prior))
})

test_that("the prior mass below the floor is known before sampling", {
  skip_if_not_installed("brms")
  set.seed(11)
  n <- 150
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1e5), y = runif(n, 0, 1e5),
                               z = stats::rnorm(n)),
                    coords = c("x", "y"), crs = 32617)
  local_mocked_bindings(estimate_sac_range = .fbr_sac(9e4), .package = "spatialkit")
  def <- .fbr_capture(d, "z", character(0))$fit
  logged <- capture_spatialkit_log(k6 <- .fbr_capture(d, "z", character(0), gp_k = 6)$fit,
                                   level = logger::WARN)
  # A derived basis sits at the prior's 1% point; gp_k = 6 leaves 98% of the
  # prior below what it resolves, and the post-fit share only restated that.
  expect_lt(def$info$gp_adequacy$prior_mass_below, 0.02)
  expect_gt(k6$info$gp_adequacy$prior_mass_below, 0.9)
  expect_equal(k6$info$gp_adequacy$prior_mass_below, 0.984, tolerance = 2e-3)
  # A coarser basis does not move the prior up to meet it.
  expect_identical(k6$info$gp_lscale_prior, def$info$gp_lscale_prior)
  expect_true(log_has(logged, "98% of the length-scale prior .* lies below the smallest length-scale gp_k = 6 resolves"))
  expect_true(log_has(logged, "the derived value is 23"))
  txt <- paste(utils::capture.output(print(k6)), collapse = "\n")
  expect_match(txt, "98% of the length-scale prior lies below what this basis resolves",
               fixed = TRUE)
  expect_false(grepl("of the length-scale prior lies below",
                     paste(utils::capture.output(print(def)), collapse = "\n"),
                     fixed = TRUE))
})

test_that(".gp_k_needed() inverts the resolution rule", {
  # A basis sized by .gp_k_needed() for a range resolves that range, and one
  # function fewer does not.
  Fse <- sqrt(2 * log(20))
  for (rng in c(8000, 20000, 45000)) {
    k <- spatialkit:::.gp_k_needed(rng, gp_c = 1.6, S = 3.5, xy_scale = c(28000, 29000))
    floor_at <- function(kk) Fse * (1.75 * 1.6 * 3.5 / kk) * 29000
    expect_lte(floor_at(k), rng)
    expect_gt(floor_at(k - 1), rng)
  }
})


# ---------------------------------------------------------------------------
# Real sampler (opt-in)
# ---------------------------------------------------------------------------

test_that("a requested resolution compiles and samples with the prior it moved", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  d <- .fbr_pts(n = 150, seed = 8, range = 25000)
  fit <- NULL
  utils::capture.output(suppressWarnings(suppressMessages(
    fit <- fit_bayesian_spatial_model(
      d, "z", character(0), gp_resolution = 25000, chains = 1, iter = 300,
      cores = 1, seed = 1, compute_loo = FALSE, check_convergence = FALSE,
      backend = "rstan"))), type = "output")
  expect_s3_class(fit, "bayesian_fit")
  expect_true(all(fit$info$gp_resolution <= 25000))
  # The prior Stan was given is the moved one.
  ps <- brms::prior_summary(fit$engine)
  expect_true(all(ps$prior[ps$class == "lscale" & nzchar(ps$coef)] ==
                    fit$info$gp_lscale_prior))
  expect_equal(fit$info$gp_lengthscale_bounds[["lower"]], fit$info$gp_ell_min)
  expect_true(all(is.finite(suppressMessages(predict(fit, newdata = d[1:5, ])))))
})
