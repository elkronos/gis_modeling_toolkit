# Real-sampler counterparts of the mocked tests in test-review2-bayes-rf.R.
# Opt-in via SPATIALKIT_TEST_BRMS, like test-bayes-smoke.R, and for the same
# reasons: each test compiles a Stan model.  Kept to two tiny fits.

.r2bb_skip <- function() {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
}

.r2bb_pts <- function(n = 40, seed = 20240817) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); z <- rnorm(n)
  resp <- 0.004 * x + 1.5 * z + rnorm(n, sd = 0.5)
  d <- sf::st_as_sf(data.frame(x = x, y = y, z = z, resp = resp),
                    coords = c("x", "y"), crs = 32632)
  d$cat3 <- cut(d$resp, stats::quantile(d$resp, c(0, 1/3, 2/3, 1)),
                include.lowest = TRUE, labels = c("lo", "mid", "hi"))
  d
}

.r2bb_fit <- function(pts, ...) {
  fit <- NULL
  warns <- character(0)
  utils::capture.output(
    withCallingHandlers(
      suppressMessages(
        fit <- fit_bayesian_spatial_model(pts, predictor_vars = "z",
                                          chains = 1, iter = 200, warmup = 100,
                                          cores = 1, compute_loo = FALSE,
                                          seed = 1234, gp_k = 6, ...)),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    type = "output")
  list(fit = fit, warnings = warns)
}

test_that("a categorical fit samples, and predict() says what it can return", {
  .r2bb_skip()
  pts <- .r2bb_pts()
  # This fit used to fail before sampling: the automatic length-scale prior
  # dropped the dpar and brms refused the duplicated rows.
  r <- .r2bb_fit(pts, response_var = "cat3", family = brms::categorical(),
                 check_convergence = FALSE)
  fit <- r$fit
  expect_s3_class(fit, "bayesian_fit")
  expect_identical(fit$info$convergence_ok, NA)
  expect_match(fit$info$gp_lscale_prior, "inv_gamma|normal")

  nd <- .r2bb_pts(n = 5, seed = 99)
  expect_error(suppressMessages(predict(fit, newdata = nd)),
               "probability per response category")
  expect_error(fitted(fit), "probability per response category")
  d <- suppressMessages(predict(fit, newdata = nd, type = "predict", draws = TRUE))
  expect_true(is.matrix(d))
  expect_identical(ncol(d), 5L)
  expect_true(all(d %in% 1:3))
})

test_that("a checked gaussian fit is quiet about capped ESS and saves its engine once", {
  .r2bb_skip()
  pts <- .r2bb_pts()
  r <- .r2bb_fit(pts, response_var = "resp", check_convergence = TRUE)
  fit <- r$fit
  expect_false(any(grepl("ESS has been capped", r$warnings, fixed = TRUE)))
  expect_true(is.logical(fit$info$convergence_ok) &&
                !is.na(fit$info$convergence_ok))

  sz <- function(x) length(serialize(x, NULL))
  engine_sz <- sz(fit$engine)
  before <- sz(fit)
  # The formula's environment was the fitting frame, which holds the brmsfit
  # again, so a saved fit was about twice its engine.
  expect_lt(before, 1.25 * engine_sz)
  invisible(fitted(fit))
  # And the fitted-value cache held the engine itself: another copy.
  expect_lt(sz(fit) - before, 0.05 * engine_sz)
  # The cache still answers for this engine: an entry keyed to it is found.
  hit <- get(".fitted_values", envir = fit$info$.cache)
  expect_true(identical(hit$engine,
                        spatialkit:::.fitted_engine_token(fit$engine)))
  expect_true(is.environment(hit$engine))
})
