# fit_bayesian_spatial_model() hands brms the response less its centre where
# it lies far from zero, and divided by a power of two where its standard
# deviation lies outside 2^-20 to 2^20, and puts every returned quantity back
# in the response's units.
#
# It handed brms the response as given.  Stan starts every parameter between
# -2 and 2 and brms floors its default priors at a scale of 2.5 in the
# response's units, so at an offset of 1e8 the posterior mean fitted values
# missed the response by the offset, in units of 1e-20 they were 1e15
# response SDs off after divergent transitions, near 1e160 the fit had no
# draws, and a response of 1e15 and the next double was sampled where a
# constant one is refused.
#
# The first tests replace brms::brm() by a recorder and read what brms is
# handed; the helpers are checked on their own.  The tests that sample are
# opt-in via SPATIALKIT_TEST_BRMS, like the other sampler tests.

.frs_pts <- function(n = 60, seed = 20261109) {
  set.seed(seed)
  u <- runif(n, 0, 4000); v <- runif(n, 0, 3000)
  # Multiples of 2^-13, so y + 1e8 and y * 2^530 are exact.
  x1 <- round(rnorm(n) * 2^13) / 2^13
  y  <- 2 + 1.5 * sin(1.2 * u / 1000) + cos(1.4 * v / 1000) + 0.6 * x1 +
    rnorm(n, sd = 0.2)
  y  <- round(y * 2^13) / 2^13
  sf::st_as_sf(data.frame(y = y, x1 = x1, X = u + 6e5, Y = v + 3.99e6),
               coords = c("X", "Y"), crs = 32617)
}

.frs_with_y <- function(d, y) { d$y <- y; d }

# fit_bayesian_spatial_model() with brms::brm() replaced by a recorder.
.frs_capture <- function(d, ...) {
  cap <- new.env()
  local_mocked_bindings(
    brm = function(...) { cap$args <- list(...); structure(list(), class = "frs_stub") },
    .package = "brms")
  cls <- character(0)
  fit <- withCallingHandlers(
    suppressMessages(fit_bayesian_spatial_model(
      d, "y", "x1", gp_k = 5L, compute_loo = FALSE, check_convergence = FALSE,
      ...)),
    warning = function(w) {
      cls <<- c(cls, class(w))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, args = cap$args, classes = cls)
}


test_that("a response of ordinary size reaches brms as it is", {
  skip_if_not_installed("brms")
  d <- .frs_pts()
  r <- .frs_capture(d)
  # The same Stan data, so the same draws under the same seed.
  expect_identical(r$args$data$y, d$y)
  expect_null(r$fit$info$response_scaling)
  expect_false("spatialkit_response_scale" %in% r$classes)
})


test_that("a response far from zero reaches brms less its centre", {
  skip_if_not_installed("brms")
  d <- .frs_pts()
  r <- .frs_capture(.frs_with_y(d, d$y + 1e8))
  rs <- r$fit$info$response_scaling
  # It was handed 1e8 + y, and no response_scaling was recorded.
  expect_false(is.null(rs))
  expect_identical(rs$scale, 1)
  expect_lt(abs(rs$centre - 1e8), 10)
  # Exact: y + 1e8 less a centre within a factor of two of every value.
  expect_identical(r$args$data$y, (d$y + 1e8) - rs$centre)
  expect_identical(r$args$data$y + rs$centre, d$y + 1e8)
  expect_lt(abs(mean(r$args$data$y)), 10)
})


test_that("a response in tiny or huge units reaches brms divided by a power of two", {
  skip_if_not_installed("brms")
  d <- .frs_pts()
  for (a in c(1e-20, 1e-160, 2^530, 1e160)) {
    ya <- d$y * a
    r  <- .frs_capture(.frs_with_y(d, ya))
    rs <- r$fit$info$response_scaling
    expect_false(is.null(rs), info = format(a))
    expect_identical(rs$centre, 0, info = format(a))
    # A power of two, so the division is exact.
    expect_identical(log2(rs$scale), round(log2(rs$scale)), info = format(a))
    expect_identical(r$args$data$y * rs$scale, ya, info = format(a))
    # The power of two nearest the SD: what brms sees has an SD within a
    # factor of sqrt(2) of 1.  It saw the SD of y times a.
    s <- stats::sd(r$args$data$y)
    expect_gt(s, 2^-0.5 * 0.999)
    expect_lt(s, 2^0.5 * 1.001)
  }
})


test_that("a response that differs only by rounding is refused as a constant one", {
  skip_if_not_installed("brms")
  d  <- .frs_pts()
  ulp <- 2^(floor(log2(1e15)) - 52)
  yr <- 1e15 + rep(c(0, 1), length.out = nrow(d)) * ulp
  cap <- new.env()
  local_mocked_bindings(
    brm = function(...) { cap$called <- TRUE; structure(list(), class = "frs_stub") },
    .package = "brms")
  # It was sampled: a fit to the last bit of 1e15.
  expect_error(suppressWarnings(suppressMessages(fit_bayesian_spatial_model(
    .frs_with_y(d, yr), "y", "x1", gp_k = 5L, compute_loo = FALSE,
    check_convergence = FALSE))),
    "takes a single value up to rounding")
  expect_null(cap$called)
  # A response with spread at the same size is fitted, about its centre.
  r <- .frs_capture(.frs_with_y(d, 1e15 + (d$y - 2) * 64))
  expect_false(is.null(r$fit$info$response_scaling))
})


test_that("a prior on the response's scale keeps the response as given, with a warning", {
  skip_if_not_installed("brms")
  d  <- .frs_pts()
  y8 <- d$y + 1e8
  r <- .frs_capture(.frs_with_y(d, y8),
                    prior = brms::set_prior("normal(1e8, 10)", class = "Intercept"))
  expect_true("spatialkit_response_scale" %in% r$classes)
  expect_identical(r$args$data$y, y8)
  expect_null(r$fit$info$response_scaling)
  # A prior on the length-scales only does not read the response's scale.
  r2 <- .frs_capture(.frs_with_y(d, y8),
                     prior = brms::set_prior("inv_gamma(5, 5)", class = "lscale"))
  expect_false("spatialkit_response_scale" %in% r2$classes)
  expect_false(is.null(r2$fit$info$response_scaling))
})


test_that("a count family is never centred or scaled", {
  skip_if_not_installed("brms")
  d <- .frs_pts()
  r <- .frs_capture(.frs_with_y(d, round(abs(d$y) * 1e7)), family = stats::poisson())
  expect_null(r$fit$info$response_scaling)
  expect_identical(r$args$data$y, round(abs(d$y) * 1e7))
})


test_that("fixef and LOO are put back in the response's units", {
  fx <- rbind(Intercept = c(0.5, 0.1, 0.3, 0.7), x1 = c(-1, 0.2, -1.4, -0.6))
  colnames(fx) <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
  rs <- list(centre = 1e8, scale = 2^-66)
  got <- .bayes_unscale_fixef(fx, rs)
  expect_identical(got["x1", ], fx["x1", ] * 2^-66)
  expect_identical(got["Intercept", "Est.Error"], 0.1 * 2^-66)
  expect_identical(got["Intercept", c("Estimate", "Q2.5", "Q97.5")],
                   fx["Intercept", c("Estimate", "Q2.5", "Q97.5")] * 2^-66 + 1e8)
  expect_identical(.bayes_unscale_fixef(fx, NULL), fx)

  skip_if_not_installed("loo")
  # The log-likelihood of y = s * y_t under N(s mu, s sigma) is that of y_t
  # under N(mu, sigma) less log(s) per row, so its LOO is the scaled one's
  # put back.
  set.seed(5)
  S <- 400; n <- 30
  y_t   <- rnorm(n)
  mu    <- matrix(rnorm(S * n, 0, 0.3), S, n)
  sigma <- exp(rnorm(S, 0, 0.1))
  s     <- 2^-66
  ll_t <- stats::dnorm(matrix(y_t, S, n, byrow = TRUE), mu, sigma, log = TRUE)
  ll_y <- stats::dnorm(matrix(s * y_t, S, n, byrow = TRUE), s * mu, s * sigma, log = TRUE)
  lo_t <- suppressWarnings(loo::loo(ll_t, r_eff = rep(1, n)))
  lo_y <- suppressWarnings(loo::loo(ll_y, r_eff = rep(1, n)))
  back <- .bayes_unscale_loo(lo_t, list(centre = 0, scale = s))
  expect_s3_class(back, "loo")
  # Tolerance: the two log-likelihoods differ by rounding of log(s) = -45.7
  # in each entry, about 1e-14 relative to the elpd.
  expect_equal(back$estimates, lo_y$estimates, tolerance = 1e-10)
  expect_equal(back$pointwise[, "elpd_loo"], lo_y$pointwise[, "elpd_loo"],
               tolerance = 1e-10)
  expect_identical(.bayes_unscale_loo(lo_t, NULL), lo_t)
})


test_that("cv_bayes() takes the predictive SD of draws at any size", {
  set.seed(2)
  x <- rnorm(200)
  expect_identical(.sd_any_size(x), stats::sd(x))
  # sd() squares the deviations: Inf at 2^530 and near 1e160, and its bits
  # lost at 2^-530 and near 1e-160 (4.5e-6 and 3.9e-5 relative error).
  # Compared as a ratio: an absolute tolerance would pass any value this
  # small.  1e-12 allows the rounding of a sum of 200 squares.
  for (a in c(2^530, 1e160, 2^-530, 1e-160))
    expect_equal(.sd_any_size(x * a) / (stats::sd(x) * a), 1, tolerance = 1e-12,
                 label = paste("sd ratio at", format(a)))
  expect_identical(.sd_any_size(c(x[1:3], NA)), NA_real_)
})


test_that("the spatial share of the GP is the same at any size of the draws", {
  # fit_bayesian_spatial_model() hands .gp_adequacy() the sdgp and sigma
  # draws in the response's units.  It squared them, so near 1e160 the share
  # was Inf / Inf and below about 1e-162 it was 0 / 0: NA, and the
  # length-scale rule decided whether the GP is identified in its place,
  # without a word.
  set.seed(31)
  S  <- 400
  ls <- cbind(lx = exp(rnorm(S, log(0.8), 0.05)), ly = exp(rnorm(S, log(0.8), 0.05)))
  sdgp  <- exp(rnorm(S, log(1.2), 0.2))
  sigma <- exp(rnorm(S, log(0.5), 0.05))
  pr <- "inv_gamma(14.08, 11.48)"
  base <- .gp_adequacy(ls, sdgp, sigma, gp_ell_min = 0.46, prior_spec = pr)
  expect_identical(base$identified_by, "spatial_share")
  for (a in c(2^530, 2^-560, 1e160, 1e-170)) {
    ga <- .gp_adequacy(ls, sdgp * a, sigma * a, gp_ell_min = 0.46, prior_spec = pr)
    lab <- format(a)
    expect_identical(ga$identified_by, base$identified_by, info = lab)
    expect_identical(ga$identified, base$identified, info = lab)
    expect_identical(ga$status, base$status, info = lab)
    # Exact for a power of two; otherwise the rounding of x * a in each draw.
    expect_equal(ga$spatial_share, base$spatial_share, tolerance = 1e-12, info = lab)
    expect_equal(ga$sdgp / a, base$sdgp, tolerance = 1e-12, info = lab)
  }
  # At ordinary sizes the share is the same to the bit as the plain formula.
  sh <- sdgp^2 / (sdgp^2 + sigma^2)
  expect_identical(unname(base$spatial_share),
                   stats::quantile(sh, c(0.025, 0.5, 0.975), names = FALSE))
})


test_that("the automatic slope prior is taken on the response brms is given", {
  skip_if_not_installed("brms")
  d  <- .frs_pts()
  ya <- d$y * 1e-20
  r  <- .frs_capture(.frs_with_y(d, ya), standardize_predictors = TRUE)
  rs <- r$fit$info$response_scaling
  expect_false(is.null(rs))
  pr <- as.data.frame(r$args$prior)
  b  <- pr$prior[pr$class == "b"]
  expect_length(b, 1L)
  # 2.5 SD of what brms is given, near 2.5: taken on the raw response it was
  # normal(0, 2.4e-20), which shrinks every slope on a response of SD 1 to 0.
  sd_b <- as.numeric(sub("^normal\\(0, (.*)\\)$", "\\1", b))
  expect_equal(sd_b, 2.5 * stats::sd(r$args$data$y), tolerance = 1e-9)
  # The record gives it in the response's units as well.
  expect_match(r$fit$info$b_prior, "in the units of 'y'", fixed = TRUE)
  rec <- as.numeric(sub("^normal\\(0, ([^)]*)\\).*$", "\\1", r$fit$info$b_prior))
  expect_equal(rec, 2.5 * stats::sd(ya), tolerance = 1e-9)
})


test_that("a prior on the slopes keeps the response as given, with a warning", {
  skip_if_not_installed("brms")
  d  <- .frs_pts()
  y8 <- d$y + 1e8
  r  <- .frs_capture(.frs_with_y(d, y8),
                     prior = brms::set_prior("normal(0, 1)", class = "b"))
  # brms reads a prior on class b in the units of the response it is given.
  expect_true("spatialkit_response_scale" %in% r$classes)
  expect_identical(r$args$data$y, y8)
  expect_null(r$fit$info$response_scaling)
})


test_that("a response under a link other than identity is never centred", {
  skip_if_not_installed("brms")
  d  <- .frs_pts()
  y8 <- d$y + 1e8
  # On the log link y - c is a different model from y, not the same one
  # moved.
  r <- .frs_capture(.frs_with_y(d, y8), family = stats::gaussian(link = "log"))
  expect_null(r$fit$info$response_scaling)
  expect_identical(r$args$data$y, y8)
})


test_that("the response is divided only where its SD lies outside 2^-20 to 2^20", {
  skip_if_not_installed("brms")
  d <- .frs_pts()
  s <- stats::sd(d$y)
  # The centre is near twice the spread, so it is never taken off here.
  inside  <- .frs_capture(.frs_with_y(d, d$y * (2^19 / s)))
  expect_null(inside$fit$info$response_scaling)
  expect_identical(inside$args$data$y, d$y * (2^19 / s))
  outside <- .frs_capture(.frs_with_y(d, d$y * (2^21 / s)))
  expect_identical(outside$fit$info$response_scaling,
                   list(centre = 0, scale = 2^21))
  small <- .frs_capture(.frs_with_y(d, d$y * (2^-21 / s)))
  expect_identical(small$fit$info$response_scaling$scale, 2^-21)
  tiny_in <- .frs_capture(.frs_with_y(d, d$y * (2^-19 / s)))
  expect_null(tiny_in$fit$info$response_scaling)
})


test_that("predictions, fitted values, coefficients and LOO come back in the response's units", {
  skip_if_not_installed("brms")
  skip_if_not_installed("loo")
  d <- .frs_pts()
  n <- nrow(d)
  # Draws on the scale brms was given: four draws of 0.1, 0.2, ... per row.
  ep <- function(object, newdata, ...)
    matrix(rep(seq_len(nrow(newdata)) / 10, each = 4), 4)
  fx <- rbind(Intercept = c(0.5, 0.1, 0.3, 0.7), x1 = c(-1, 0.2, -1.4, -0.6))
  colnames(fx) <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
  set.seed(8)
  ll <- matrix(stats::dnorm(rnorm(200 * n), log = TRUE), 200, n)
  lo <- suppressWarnings(loo::loo(ll, r_eff = rep(1, n)))
  local_mocked_bindings(
    brm = function(...) structure(list(), class = "brmsfit"),
    posterior_epred = ep, posterior_predict = ep,
    fixef = function(x, ...) fx,
    loo = function(x, ...) lo,
    .package = "brms")
  for (cs in list(list(a = 2^530, b = 0), list(a = 1, b = 1e8))) {
    lab <- sprintf("y * %s + %s", format(cs$a), format(cs$b))
    fit <- suppressWarnings(suppressMessages(fit_bayesian_spatial_model(
      .frs_with_y(d, d$y * cs$a + cs$b), "y", "x1", gp_k = 5L,
      compute_loo = TRUE, check_convergence = FALSE)))
    rs <- fit$info$response_scaling
    expect_false(is.null(rs), info = lab)
    want <- rs$centre + rs$scale * seq_len(n) / 10
    expect_identical(fitted(fit), want, info = lab)
    expect_identical(predict(fit, newdata = d), want, info = lab)
    expect_identical(predict(fit, newdata = d, type = "predict"), want, info = lab)
    dr <- predict(fit, newdata = d, draws = TRUE)
    expect_identical(dr[1L, ], want, info = lab)
    cf <- coef(fit)
    expect_identical(cf["x1", ], fx["x1", ] * rs$scale, info = lab)
    expect_identical(cf["Intercept", "Estimate"], 0.5 * rs$scale + rs$centre,
                     info = lab)
    # looic moves by 2 n log(scale) and not at all for a centre.
    expect_equal(fit$info$looic,
                 lo$estimates["looic", "Estimate"] + 2 * n * log(rs$scale),
                 tolerance = 1e-12, info = lab)
  }
})


# ---- exgaussian() and asym_laplace() -----------------------------------------
#
# Both are location-scale families on the response's own scale under their
# identity link: a shift of the response moves mu only, and a scale
# multiplies mu, sigma and the exponential mean beta of exgaussian() and
# leaves the quantile of asym_laplace() as it is.  They were handed to brms
# as given, so they met the offsets and units the three families above were
# taken out of.

.frs_ls_families <- function() {
  list(exgaussian = brms::exgaussian(), asym_laplace = brms::asym_laplace())
}

test_that("exgaussian() and asym_laplace() hand brms the response centred and scaled", {
  skip_if_not_installed("brms")
  d <- .frs_pts()
  for (nm in names(.frs_ls_families())) {
    fam <- .frs_ls_families()[[nm]]
    # Far from zero: less its centre, exactly.
    r8 <- .frs_capture(.frs_with_y(d, d$y + 1e8), family = fam)
    rs <- r8$fit$info$response_scaling
    expect_false(is.null(rs), info = nm)
    expect_identical(rs$scale, 1, info = nm)
    expect_identical(r8$args$data$y, (d$y + 1e8) - rs$centre, info = nm)
    expect_lt(abs(mean(r8$args$data$y)), 10, label = paste("mean,", nm))
    # In tiny or huge units: divided by the power of two nearest the SD.
    for (a in c(1e-20, 2^530)) {
      ra <- .frs_capture(.frs_with_y(d, d$y * a), family = fam)
      rs <- ra$fit$info$response_scaling
      lab <- paste(nm, format(a))
      expect_false(is.null(rs), info = lab)
      expect_identical(rs$centre, 0, info = lab)
      expect_identical(log2(rs$scale), round(log2(rs$scale)), info = lab)
      expect_identical(ra$args$data$y * rs$scale, d$y * a, info = lab)
    }
    # A response of ordinary size is handed over as it is.
    r1 <- .frs_capture(d, family = fam)
    expect_identical(r1$args$data$y, d$y, info = nm)
    expect_null(r1$fit$info$response_scaling, info = nm)
  }
  # Named as a string, as brms takes it.
  rs <- .frs_capture(.frs_with_y(d, d$y + 1e8), family = "exgaussian")
  expect_false(is.null(rs$fit$info$response_scaling))
})


test_that("a prior on the quantile keeps the transform and a prior on beta does not", {
  skip_if_not_installed("brms")
  d  <- .frs_pts()
  y8 <- d$y + 1e8
  # The quantile does not change under a shift or a scale of the response,
  # so its prior reaches brms as it was given and the response is centred.
  pq <- brms::set_prior("beta(2, 2)", class = "quantile")
  r  <- .frs_capture(.frs_with_y(d, y8), family = brms::asym_laplace(), prior = pq)
  expect_false("spatialkit_response_scale" %in% r$classes)
  expect_false(is.null(r$fit$info$response_scaling))
  pr <- as.data.frame(r$args$prior)
  expect_identical(pr$prior[pr$class == "quantile"], "beta(2, 2)")
  # beta, the mean of exgaussian()'s exponential part, is in the response's
  # units, so brms reads its prior on the response it is given.
  pb <- brms::set_prior("gamma(1, 1e-8)", class = "beta")
  r2 <- .frs_capture(.frs_with_y(d, y8), family = brms::exgaussian(), prior = pb)
  expect_true("spatialkit_response_scale" %in% r2$classes)
  expect_identical(r2$args$data$y, y8)
  expect_null(r2$fit$info$response_scaling)
})


test_that("exgaussian() and asym_laplace() take a response off its centre and scale at moderate sizes", {
  skip_if_not_installed("brms")
  # brms's default prior on exgaussian()'s beta, gamma(1, 0.1), and Stan's
  # starting values do not follow the response.  Handed over as given, an
  # exgaussian fit of SD 136 came back with beta at a twentieth of its
  # value, one of SD 1,360 or 300 SDs from zero had no draws, and an
  # asym_laplace fit 300 SDs from zero or of SD 1e-4 had sigma at twice its
  # value or more.  So these two families are centred at more than 4
  # spreads from zero and scaled outside an SD of 2^-2 to 2^2.
  d <- .frs_pts()
  s <- stats::sd(d$y)                        # 1.28; spread 3.42, centre 2.09
  for (nm in names(.frs_ls_families())) {
    fam <- .frs_ls_families()[[nm]]
    for (a in c(8, 100, 1e3, 1e-4, 0.15)) {
      ya <- d$y * a
      r  <- .frs_capture(.frs_with_y(d, ya), family = fam)
      rs <- r$fit$info$response_scaling
      lab <- paste(nm, format(a))
      expect_false(is.null(rs), info = lab)
      expect_identical(rs$centre, 0, info = lab)
      expect_identical(log2(rs$scale), round(log2(rs$scale)), info = lab)
      expect_identical(r$args$data$y * rs$scale, ya, info = lab)
      # The power of two nearest the SD.
      expect_lt(abs(log2(stats::sd(r$args$data$y))), 0.5 + 1e-9, label = lab)
    }
    for (b in c(-300, 20)) {
      yb <- d$y + b
      r  <- .frs_capture(.frs_with_y(d, yb), family = fam)
      rs <- r$fit$info$response_scaling
      lab <- paste(nm, "+", format(b))
      expect_false(is.null(rs), info = lab)
      expect_identical(rs$scale, 1, info = lab)
      # Exact: every value lies within a factor of two of the centre.
      expect_identical(r$args$data$y + rs$centre, yb, info = lab)
      expect_lt(abs(stats::median(r$args$data$y)), 1e-9, label = lab)
    }
    # Within 4 spreads of zero and an SD of 2^-2 to 2^2: as given.
    for (yi in list(d$y * 2, d$y + 10, d$y * 0.5)) {
      r <- .frs_capture(.frs_with_y(d, yi), family = fam)
      expect_identical(r$args$data$y, yi, info = nm)
      expect_null(r$fit$info$response_scaling, info = nm)
    }
  }
  # The rule of the other families does not change: in band, as given.
  for (yi in list(d$y * 100, d$y + 300, d$y * 1e-4)) {
    r <- .frs_capture(.frs_with_y(d, yi))
    expect_identical(r$args$data$y, yi)
    expect_null(r$fit$info$response_scaling)
  }
  # The message names the rule that applied.
  rs <- .bayes_response_scaling(d$y * 1e3, brms::exgaussian())
  expect_match(.bayes_response_scaling_text(rs), "outside 2^-2 to 2^2", fixed = TRUE)
  rs <- .bayes_response_scaling(d$y - 300, brms::asym_laplace())
  expect_match(.bayes_response_scaling_text(rs), "more than 4 times its spread",
               fixed = TRUE)
  rs <- .bayes_response_scaling(d$y * 1e30, stats::gaussian())
  expect_match(.bayes_response_scaling_text(rs), "outside 2^-20 to 2^20", fixed = TRUE)
})


test_that("a prior on sigma keeps the response as given, with a warning", {
  skip_if_not_installed("brms")
  d  <- .frs_pts()
  y8 <- d$y + 1e8
  # sigma is in the response's units under every one of these families, so
  # brms reads its prior on the response it is given.
  ps <- brms::set_prior("student_t(3, 0, 2.5)", class = "sigma")
  fams <- c(list(gaussian = stats::gaussian()), .frs_ls_families())
  for (nm in names(fams)) {
    r <- .frs_capture(.frs_with_y(d, y8), family = fams[[nm]], prior = ps)
    expect_true("spatialkit_response_scale" %in% r$classes, info = nm)
    expect_identical(r$args$data$y, y8, info = nm)
    expect_null(r$fit$info$response_scaling, info = nm)
  }
})


test_that("lognormal() and zero_inflated_asym_laplace are never centred or scaled", {
  skip_if_not_installed("brms")
  d <- .frs_pts()
  # lognormal() reports an identity link, but its mu is on the log scale,
  # so y - c is another model; a shift moves the zero mass of the
  # zero-inflated family.
  fams <- list(lognormal = brms::lognormal(),
               zi_asym_laplace = brms::brmsfamily("zero_inflated_asym_laplace"))
  for (nm in names(fams)) {
    for (yi in list(d$y + 1e8, (d$y + 1) * 1e-20)) {
      r <- .frs_capture(.frs_with_y(d, yi), family = fams[[nm]])
      expect_identical(r$args$data$y, yi, info = nm)
      expect_null(r$fit$info$response_scaling, info = nm)
    }
  }
})


test_that("a Stan fit with no posterior draws is refused", {
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  # rstan returns an empty stanfit when it rejects every starting value.
  # The fit was returned with convergence_ok NA and no warning, and
  # fitted() stopped later.
  empty <- methods::new(methods::getClass("stanfit", where = asNamespace("rstan")))
  local_mocked_bindings(
    brm = function(...) structure(list(fit = empty), class = "brmsfit"),
    .package = "brms")
  d <- .frs_pts()
  for (chk in c(TRUE, FALSE)) {
    expect_error(suppressWarnings(suppressMessages(fit_bayesian_spatial_model(
      .frs_with_y(d, d$y * 1e3), "y", "x1", family = brms::exgaussian(),
      gp_k = 5L, compute_loo = TRUE, check_convergence = chk))),
      "no posterior draws.*less 0 and divided by 1024, with median 2\\.037.*'exgaussian'")
  }
  expect_error(suppressWarnings(suppressMessages(fit_bayesian_spatial_model(
    d, "y", "x1", gp_k = 5L, compute_loo = FALSE))),
    "response 'y' as it is, with median 2\\.086 and standard deviation 1\\.279")
})


test_that("exgaussian() and asym_laplace() fits come back in the response's units", {
  skip_if_not_installed("brms")
  skip_if_not_installed("loo")
  d <- .frs_pts()
  n <- nrow(d)
  ep <- function(object, newdata, ...)
    matrix(rep(seq_len(nrow(newdata)) / 10, each = 4), 4)
  fx <- rbind(Intercept = c(0.5, 0.1, 0.3, 0.7), x1 = c(-1, 0.2, -1.4, -0.6))
  colnames(fx) <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
  set.seed(8)
  ll <- matrix(stats::dnorm(rnorm(200 * n), log = TRUE), 200, n)
  lo <- suppressWarnings(loo::loo(ll, r_eff = rep(1, n)))
  local_mocked_bindings(
    brm = function(...) structure(list(), class = "brmsfit"),
    posterior_epred = ep, posterior_predict = ep,
    fixef = function(x, ...) fx,
    loo = function(x, ...) lo,
    .package = "brms")
  for (nm in names(.frs_ls_families())) {
    for (cs in list(list(a = 1e-20, b = 0), list(a = 1, b = 1e8))) {
      lab <- sprintf("%s, y * %s + %s", nm, format(cs$a), format(cs$b))
      yg  <- d$y * cs$a + cs$b
      fit <- suppressWarnings(suppressMessages(fit_bayesian_spatial_model(
        .frs_with_y(d, yg), "y", "x1", gp_k = 5L,
        family = .frs_ls_families()[[nm]],
        compute_loo = TRUE, check_convergence = FALSE)))
      rs <- fit$info$response_scaling
      expect_false(is.null(rs), info = lab)
      want <- rs$centre + rs$scale * seq_len(n) / 10
      expect_identical(fitted(fit), want, info = lab)
      expect_identical(predict(fit, newdata = d), want, info = lab)
      expect_identical(predict(fit, newdata = d, type = "predict"), want, info = lab)
      expect_identical(predict(fit, newdata = d, draws = TRUE)[1L, ], want, info = lab)
      expect_identical(unname(residuals(fit)), unname(yg - want), info = lab)
      cf <- coef(fit)
      expect_identical(cf["x1", ], fx["x1", ] * rs$scale, info = lab)
      expect_identical(cf["Intercept", "Estimate"], 0.5 * rs$scale + rs$centre,
                       info = lab)
      expect_equal(fit$info$looic,
                   lo$estimates["looic", "Estimate"] + 2 * n * log(rs$scale),
                   tolerance = 1e-12, info = lab)
    }
  }
})


# ---- Sampling (opt-in) -------------------------------------------------------

.frs_fit <- function(d) {
  fit <- NULL
  utils::capture.output(suppressWarnings(suppressMessages(
    fit <- fit_bayesian_spatial_model(
      d, "y", "x1", chains = 1, iter = 400, warmup = 200, gp_k = 5L,
      cores = 1, compute_loo = TRUE, check_convergence = TRUE,
      backend = "rstan", seed = 11))), type = "output")
  fit
}

# Fitted values on the response's own scale: no change of origin or unit of
# the response moves them.  Mean and SD taken on y / max|y|, so that a
# response near 1e160 does not overflow.
.frs_z <- function(fv, y) {
  a <- max(abs(y))
  (fv - mean(y / a) * a) / (stats::sd(y / a) * a)
}

test_that("the fit is the same from any origin and in any unit of the response", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  d  <- .frs_pts()
  f0 <- .frs_fit(d)
  z0 <- .frs_z(fitted(f0), d$y)
  c0 <- coef(f0)
  n0 <- f0$info$convergence_diagnostics$n_divergent
  cases <- list(list(a = 1, b = 1e8), list(a = 1e-20, b = 0), list(a = 2^530, b = 0))
  for (cs in cases) {
    lab <- sprintf("y * %s + %s", format(cs$a), format(cs$b))
    yg  <- d$y * cs$a + cs$b
    fg  <- .frs_fit(.frs_with_y(d, yg))
    fv  <- fitted(fg)
    expect_true(all(is.finite(fv)), info = lab)
    # Monte Carlo allowance, in response SDs: every row within 0.5 and the
    # rows on average within 0.15 of the base fit.  It was off by 7e7 at an
    # offset of 1e8 and 1e14 in units of 1e-20, and near 1e160 had no draws.
    dz <- abs(.frs_z(fv, yg) - z0)
    expect_lt(max(dz), 0.5, label = paste("max |dz|,", lab))
    expect_lt(mean(dz), 0.15, label = paste("mean |dz|,", lab))
    # predict() gives the same values on the training rows, in the units.
    pz <- .frs_z(predict(fg, newdata = .frs_with_y(d, yg)), yg)
    expect_lt(mean(abs(pz - z0)), 0.15, label = paste("predict,", lab))
    # The slope in the units of the response; the allowance is 0.15
    # response SDs per unit of x1, as for the fitted values.
    cg <- coef(fg)
    expect_lt(abs(cg["x1", "Estimate"] / cs$a - c0["x1", "Estimate"]),
              0.15 * stats::sd(d$y), label = paste("slope,", lab))
    # No divergent transition where the base fit has none (103 in units of
    # 1e-20).
    if (identical(n0, 0))
      expect_identical(fg$info$convergence_diagnostics$n_divergent, 0,
                       info = lab)
    # LOO on the response's scale: looic moves by 2 n log(a), here up to
    # 1e5, and is otherwise the base fit's to Monte Carlo error.
    expect_lt(abs(fg$info$looic - 2 * nrow(d) * log(cs$a) - f0$info$looic), 10,
              label = paste("looic,", lab))
    expect_equal(fg$info$gp_adequacy$sdgp / cs$a, f0$info$gp_adequacy$sdgp,
                 tolerance = 0.5, label = paste("sdgp,", lab))
    # The same rule decides whether the GP is identified: the share, which
    # near 1e160 was NA and handed the decision to the length-scales.
    expect_identical(fg$info$gp_adequacy$identified_by,
                     f0$info$gp_adequacy$identified_by, info = lab)
    expect_true(all(is.finite(fg$info$gp_adequacy$spatial_share)), info = lab)
  }
})


test_that("cv_bayes() scores the same in any unit of the response", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  d <- .frs_pts()
  folds <- rep(1:2, length.out = nrow(d))
  run <- function(y) {
    out <- NULL
    utils::capture.output(suppressWarnings(suppressMessages(
      out <- cv_bayes(.frs_with_y(d, y), "y", "x1", folds = folds, k = 2,
                      seed = 123,
                      fit_args = list(chains = 1, iter = 400, warmup = 200,
                                      gp_k = 5L, cores = 1, backend = "rstan",
                                      check_convergence = FALSE)))),
      type = "output")
    out
  }
  r0 <- run(d$y)
  for (a in c(1e-20, 2^530)) {
    ra <- run(d$y * a)
    # Every row predicted, as in the base run; near 1e160 none was.
    expect_identical(sum(is.finite(ra$predictions$yhat)),
                     sum(is.finite(r0$predictions$yhat)), info = format(a))
    # Pooled RMSE within 15% of the base run's after the scale (6.4e-5
    # against 2.7e-21 in units of 1e-20).
    expect_lt(abs(ra$overall$RMSE / (a * r0$overall$RMSE) - 1), 0.15,
              label = paste("RMSE ratio,", format(a)))
    cov <- unlist(ra$predictive_coverage)
    expect_true(all(is.finite(cov)), info = format(a))
    expect_true(all(is.finite(ra$predictions$yhat_sd)), info = format(a))
  }
})


# The data each family fits: an exponential part of mean 0.5 for
# exgaussian(), whose beta has nothing to hold on Gaussian noise (its looic
# then moved by 28 between seeds of the same fit), and Gaussian noise for
# asym_laplace().  Multiples of 2^-13, so y + 1e8 and y - 300 are exact.
.frs_ls_data <- function() {
  d0 <- .frs_pts()
  set.seed(20261110)
  e  <- round(stats::rexp(nrow(d0), rate = 2) * 2^13) / 2^13
  list(exgaussian = .frs_with_y(d0, d0$y + e), asym_laplace = d0)
}

.frs_ls_fit <- function(dd, fam, seed, prior = NULL) {
  fit <- NULL
  utils::capture.output(suppressWarnings(suppressMessages(
    fit <- fit_bayesian_spatial_model(
      dd, "y", "x1", family = fam, prior = prior, chains = 1, iter = 400,
      warmup = 200, gp_k = 5L, cores = 1, compute_loo = TRUE,
      check_convergence = TRUE, backend = "rstan", seed = seed))),
    type = "output")
  fit
}

# Posterior mean of a distributional parameter in the response's units.
# brms::as_draws_df() is posterior's, which brms exports: posterior is not
# among the package's Suggests, and R CMD check warns on a `::` to it.
.frs_dpar <- function(fit, nm) {
  dr <- brms::as_draws_df(fit$engine)
  s  <- fit$info$response_scaling$scale %||% 1
  mean(dr[[nm]]) * s
}

test_that("exgaussian() and asym_laplace() fits are the same from any origin and in any unit", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  dd <- .frs_ls_data()
  # One chain of 400 iterations.  Under seed 11 the exgaussian fit of y
  # itself settles 0.17 response SDs (mean) from the fits of seed 12 and of
  # the transformed responses under either seed, which agree within 0.01,
  # at 400 and at 1000 iterations alike.  So the fitted values are held to
  # 0.25 SDs of the base fit on average: a seed that moves the base chain
  # as far as seed 11 still passes, and the misses this test is there for
  # (1.6 SDs at y + 300, 4e5 to 2e15 SDs at y + 1e8 and y * 1e-20) do not.
  seeds <- c(exgaussian = 12, asym_laplace = 11)
  # Response sizes: far offsets and units (1e8, 1e-20), and moderate ones
  # inside 2^-20 to 2^20 where these two families failed as given.  Under
  # exgaussian, y * 1e3 and y - 300 had no draws and y * 100 gave beta at a
  # twentieth of its value; under asym_laplace, y - 300 and y * 1e-4 gave
  # sigma at twice its value or more and looic 40 to 150 higher.
  cases <- list(list(a = 1, b = 1e8), list(a = 1e-20, b = 0),
                list(a = 100, b = 0), list(a = 1e3, b = 0),
                list(a = 1e-4, b = 0), list(a = 1, b = -300))
  dpars <- list(exgaussian = c("sigma", "beta"), asym_laplace = "sigma")
  for (nm in names(.frs_ls_families())) {
    fam <- .frs_ls_families()[[nm]]
    d   <- dd[[nm]]
    f0  <- .frs_ls_fit(d, fam, seeds[[nm]])
    z0  <- .frs_z(fitted(f0), d$y)
    p0  <- vapply(dpars[[nm]], function(p) .frs_dpar(f0, p), numeric(1))
    for (cs in cases) {
      lab <- sprintf("%s, y * %s + %s", nm, format(cs$a), format(cs$b))
      yg  <- d$y * cs$a + cs$b
      fg  <- .frs_ls_fit(.frs_with_y(d, yg), fam, seeds[[nm]])
      fv  <- fitted(fg)
      expect_true(all(is.finite(fv)), info = lab)
      dz <- abs(.frs_z(fv, yg) - z0)
      expect_lt(max(dz), 0.75, label = paste("max |dz|,", lab))
      expect_lt(mean(dz), 0.25, label = paste("mean |dz|,", lab))
      expect_lt(abs(coef(fg)["x1", "Estimate"] / cs$a - coef(f0)["x1", "Estimate"]),
                0.25 * stats::sd(d$y), label = paste("slope,", lab))
      # sigma and beta in the response's units, within 30% of the base fit's:
      # over the transformed fits they ranged 0.16 to 0.18 and 0.45 to 0.46
      # for exgaussian, 0.11 to 0.12 for asym_laplace's sigma.
      pg <- vapply(dpars[[nm]], function(p) .frs_dpar(fg, p) / cs$a, numeric(1))
      expect_lt(max(abs(pg / p0 - 1)), 0.3,
                label = paste(paste(dpars[[nm]], collapse = " and "), lab))
      # looic moves by 2 n log(a).  Allowance 30: the exgaussian fit of one
      # response gave looic from 69 to 88 over three seeds of this sampler
      # (Pareto k up to 1.1 on 200 draws).  It was off by 5e4 at an offset
      # of 1e8 and by 4e3 in units of 1e-20.
      expect_lt(abs(fg$info$looic - 2 * nrow(d) * log(cs$a) - f0$info$looic), 30,
                label = paste("looic,", lab))
    }
  }
})


test_that("a sampler that rejects every starting value is an error", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  # A prior on beta keeps the response as given.  At y * 1e3 the
  # exgaussian density is zero at every starting value Stan tries, and the
  # empty fit was returned with convergence_ok NA and no warning.
  d  <- .frs_ls_data()$exgaussian
  pb <- brms::set_prior("gamma(1, 0.1)", class = "beta")
  expect_error(.frs_ls_fit(.frs_with_y(d, d$y * 1e3), brms::exgaussian(), 12,
                           prior = pb),
               "no posterior draws")
})
