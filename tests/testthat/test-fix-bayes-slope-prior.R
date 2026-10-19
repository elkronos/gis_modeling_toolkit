# The automatic slope prior of fit_bayesian_spatial_model(standardize_predictors
# = TRUE) follows the scale of the response.
#
# It was a fixed normal(0, 5) on every slope.  Under an identity link a slope
# per SD of a predictor is in the units of the response, so the same data gave
# a slope of 0.3 with income in dollars and 10.1 (thousand) with income in
# thousands, with convergence_ok TRUE both times.
#
# Nothing here compiles a Stan model except the last test, which is opt-in via
# SPATIALKIT_TEST_BRMS like the other sampler tests: brms::brm() is replaced by
# a recorder and brms::validate_prior() / make_stancode() run for real.

.fbs_pts <- function(n = 150, seed = 11) {
  set.seed(seed)
  x <- runif(n, 0, 1e5); y <- runif(n, 0, 1e5)
  edu <- rnorm(n, 30, 10)
  inc <- 60000 + 1500 * (edu - 30) + 4000 * sin(x / 3e4) + rnorm(n, 0, 8000)
  sf::st_as_sf(
    data.frame(x, y, edu, inc, inc_k = inc / 1000,
               cnt = stats::rpois(n, 5), bin = as.integer(edu > 30),
               pos = exp(inc / 3e4)),
    coords = c("x", "y"), crs = 32617)
}

.fbs_capture <- function(...) {
  cap <- new.env()
  local_mocked_bindings(
    brm = function(...) { cap$args <- list(...); structure(list(), class = "fbs_stub") },
    .package = "brms")
  fit <- suppressWarnings(suppressMessages(
    fit_bayesian_spatial_model(..., compute_loo = FALSE,
                               check_convergence = FALSE)))
  list(fit = fit, args = cap$args)
}

.fbs_b_prior <- function(args) {
  pr <- as.data.frame(args$prior)
  pr$prior[pr$class == "b"]
}

.fbs_scale <- function(spec) {
  as.numeric(sub("^normal\\(0, (.*)\\)$", "\\1", spec))
}


test_that("the automatic slope prior scales with the response under an identity link", {
  skip_if_not_installed("brms")
  d <- .fbs_pts()
  r_usd <- .fbs_capture(d, "inc",   "edu", standardize_predictors = TRUE, gp_k = 6)
  r_k   <- .fbs_capture(d, "inc_k", "edu", standardize_predictors = TRUE, gp_k = 6)

  b_usd <- .fbs_b_prior(r_usd$args)
  b_k   <- .fbs_b_prior(r_k$args)
  # 2.5 SDs of the response, in the response's own units ...
  expect_equal(.fbs_scale(b_usd), 2.5 * stats::sd(d$inc), tolerance = 1e-8)
  expect_equal(.fbs_scale(b_k),   2.5 * stats::sd(d$inc_k), tolerance = 1e-8)
  # ... so rescaling the response rescales the prior with it.  Both were
  # "normal(0, 5)".
  expect_equal(.fbs_scale(b_usd) / .fbs_scale(b_k), 1000, tolerance = 1e-8)

  # What reaches Stan carries that scale, not a fixed 5.
  code <- brms::make_stancode(r_usd$args$formula, data = r_usd$args$data,
                              family = r_usd$args$family,
                              prior = r_usd$args$prior)
  expect_true(grepl(sprintf("normal_lpdf(b | 0, %s)",
                            sub("^normal\\(0, (.*)\\)$", "\\1", b_usd)),
                    code, fixed = TRUE))
  expect_false(grepl("normal_lpdf(b | 0, 5)", code, fixed = TRUE))

  # The fit records the prior brms used, and print() shows it.
  expect_identical(r_usd$fit$info$b_prior, b_usd)
  expect_match(r_usd$fit$info$b_prior_source, "automatic")
  txt <- paste(utils::capture.output(print(r_usd$fit)), collapse = "\n")
  expect_match(txt, paste0("Slope prior: ", b_usd), fixed = TRUE)
})

test_that("the slope prior stays on a fixed unit scale for log and logit links", {
  skip_if_not_installed("brms")
  d <- .fbs_pts()
  expect_identical(.fbs_b_prior(.fbs_capture(
    d, "cnt", "edu", family = stats::poisson(), standardize_predictors = TRUE,
    gp_k = 6)$args), "normal(0, 5)")
  expect_identical(.fbs_b_prior(.fbs_capture(
    d, "bin", "edu", family = brms::bernoulli(), standardize_predictors = TRUE,
    gp_k = 6)$args), "normal(0, 5)")
  # An identity-link family fitted through a log link is on the link's scale.
  expect_identical(.fbs_b_prior(.fbs_capture(
    d, "pos", "edu", family = stats::gaussian(link = "log"),
    standardize_predictors = TRUE, gp_k = 6)$args), "normal(0, 5)")
})

test_that(".b_prior_auto() picks the scale the linear predictor is on", {
  skip_if_not_installed("brms")
  set.seed(4)
  y <- stats::rnorm(80, 6e4, 1.5e4)
  yp <- exp(stats::rnorm(80, 1, 0.6))
  auto <- spatialkit:::.b_prior_auto

  g <- auto(y, stats::gaussian())
  expect_identical(g$basis, "sd")
  expect_equal(.fbs_scale(g$spec), 2.5 * stats::sd(y), tolerance = 1e-8)
  expect_identical(auto(y, brms::student())$spec, g$spec)
  expect_identical(auto(y, brms::skew_normal())$spec, g$spec)
  # A family given as a bare string or as a function resolves the same way.
  expect_identical(auto(y, "gaussian")$spec, g$spec)
  expect_identical(auto(y, stats::gaussian)$spec, g$spec)
  # A mixture of identity-link components is on the response's scale too.
  mix <- suppressMessages(brms::mixture(stats::gaussian(), stats::gaussian()))
  expect_identical(auto(y, mix)$spec, g$spec)

  # lognormal reports an identity link, and its location is on the log scale:
  # sd(y) would be the wrong scale by a factor of sd(y) / sd(log y).
  ln <- auto(yp, brms::lognormal())
  expect_identical(ln$basis, "sd_log")
  expect_equal(.fbs_scale(ln$spec), 2.5 * stats::sd(log(yp)), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(.fbs_scale(ln$spec), 2.5 * stats::sd(yp))))

  # Non-identity links keep the fixed prior.
  for (fam in list(stats::poisson(), brms::bernoulli(), brms::negbinomial(),
                   stats::Gamma(link = "log"), brms::categorical(),
                   brms::cumulative(), stats::gaussian(link = "log")))
    expect_identical(auto(y, fam)$spec, "normal(0, 5)")

  # A constant response has no spread; normal(0, 0) is not a prior Stan takes.
  k <- auto(rep(3, 20), stats::gaussian())
  expect_gt(.fbs_scale(k$spec), 0)
  # A small scale is not rounded to a literal zero.
  expect_gt(.fbs_scale(auto(y * 1e-9, stats::gaussian())$spec), 0)
})

test_that("a user prior, or no standardisation, sets no automatic slope prior", {
  skip_if_not_installed("brms")
  d <- .fbs_pts()
  # Not standardised: brms's flat default, recorded as such.
  r0 <- .fbs_capture(d, "inc", "edu", gp_k = 6)
  expect_length(.fbs_b_prior(r0$args), 0L)
  expect_identical(r0$fit$info$b_prior, "(flat)")
  expect_identical(r0$fit$info$b_prior_source, "brms default")
  # A user's own prior is kept as given and recorded as theirs.
  up <- brms::set_prior("normal(0, 123)", class = "b")
  ru <- .fbs_capture(d, "inc", "edu", standardize_predictors = TRUE, gp_k = 6,
                     prior = up)
  expect_identical(.fbs_b_prior(ru$args), "normal(0, 123)")
  expect_identical(ru$fit$info$b_prior, "normal(0, 123)")
  expect_identical(ru$fit$info$b_prior_source, "from `prior`")
  # No slopes, no slope prior.
  ri <- .fbs_capture(d, "inc", character(0), gp_k = 6)
  expect_true(is.na(ri$fit$info$b_prior))
})


# ---------------------------------------------------------------------------
# Real sampler (opt-in)
# ---------------------------------------------------------------------------

test_that("the slope is recovered on a dollar-scale response and follows its units", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  d <- .fbs_pts()
  ols <- stats::coef(stats::lm(inc ~ scale(edu), data = sf::st_drop_geometry(d)))[[2L]]
  fit_one <- function(resp) {
    fit <- NULL
    utils::capture.output(suppressWarnings(suppressMessages(
      fit <- fit_bayesian_spatial_model(
        d, resp, "edu", standardize_predictors = TRUE, gp_k = 6,
        chains = 2, iter = 600, cores = 1, seed = 1, compute_loo = FALSE,
        backend = "rstan"))), type = "output")
    fit
  }
  f_usd <- fit_one("inc")
  f_k   <- fit_one("inc_k")
  b_usd <- coef(f_usd)["edu", ]
  b_k   <- coef(f_k)["edu", ]
  # The fixed prior gave 0.21 (95% interval -9.6 to 10.7) against an OLS slope
  # of about 14,000 per SD.
  expect_lt(b_usd[["Q2.5"]], ols)
  expect_gt(b_usd[["Q97.5"]], ols)
  expect_gt(b_usd[["Estimate"]], 0.8 * ols)
  # Dividing the response by 1000 divides the estimate by 1000, to Monte Carlo
  # error (it was 0.21 against 13.7).
  # Tolerance: half a posterior SD, five Monte Carlo standard errors at an
  # effective sample size of 100.
  expect_equal(b_usd[["Estimate"]] / 1000, b_k[["Estimate"]],
               tolerance = 0.5 * b_k[["Est.Error"]] / abs(b_k[["Estimate"]]))
})
