# Smoke test for the brms backend.
#
# Exactly one other test in the suite touches a brms code path, and the brms
# examples live in \dontrun{}, so R CMD check never runs them either.  Without
# this file the weekly check-brms workflow proves only that brms *installs* --
# it never samples a posterior.  This fits one deliberately tiny model so the
# Bayesian path is exercised end to end.
#
# Opt-in only, via SPATIALKIT_TEST_BRMS.  skip_on_cran() is NOT enough: the
# r-lib GitHub actions set NOT_CRAN=true, so these five Stan fits would run in
# every matrix job that happens to have brms available -- slow, and dependent
# on whether Stan sampled cleanly on that runner rather than on anything this
# package controls.  check-brms.yaml sets the variable; nothing else does.
#
# NOTE on argument names: fit_bayesian_spatial_model() has no `...`, so only
# its actual formals may be passed -- there is no `refresh`.  cv_bayes()
# likewise forwards fit arguments through `fit_args`, not as top-level
# arguments.  Stan's own sampling chatter goes to stdout via cat(), which
# suppressMessages() cannot catch, hence capture.output().

bayes_smoke_points <- function(n = 40, seed = 20240817) {
  set.seed(seed)
  x <- runif(n, 0, 1000)
  y <- runif(n, 0, 1000)
  # Response with a mild spatial trend plus a real predictor effect.
  z <- rnorm(n)
  resp <- 0.004 * x + 1.5 * z + rnorm(n, sd = 0.5)
  sf::st_as_sf(
    data.frame(x = x, y = y, z = z, resp = resp),
    coords = c("x", "y"), crs = 32632
  )
}

# Quietly fit the smallest useful model.  `...` here reaches
# fit_bayesian_spatial_model(), so callers may only pass its real formals.
fit_smoke_model <- function(pts, ...) {
  fit <- NULL
  utils::capture.output(
    suppressWarnings(suppressMessages(
      fit <- fit_bayesian_spatial_model(
        pts, response_var = "resp", predictor_vars = "z",
        chains = 1, iter = 200, warmup = 100,
        cores = 1, compute_loo = FALSE, check_convergence = FALSE,
        seed = 1234, ...
      )
    )),
    type = "output"
  )
  fit
}

test_that("fit_bayesian_spatial_model returns a usable bayesian_fit", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan smoke tests")
  skip_if_not_installed("brms")

  pts <- bayes_smoke_points()
  fit <- fit_smoke_model(pts)

  expect_s3_class(fit, "bayesian_fit")
  expect_s3_class(fit, "spatial_fit")
  expect_identical(fit$response_var, "resp")
  expect_identical(fit$predictor_vars, "z")
  expect_identical(fit$n, nrow(pts))

  # The GP basis metadata is the package's own contribution over a bare
  # brms::gp() call, so assert its internal consistency rather than a value.
  expect_true(is.numeric(fit$info$gp_k) && fit$info$gp_k >= 1)
  expect_identical(fit$info$gp_n_basis, as.numeric(fit$info$gp_k)^2)
  expect_true(is.character(fit$info$gp_lscale_prior))
  expect_match(fit$info$gp_lscale_prior, "inv_gamma|normal")
  # Whichever prior was chosen, its parameters must be strictly positive --
  # a regression here means sprintf() rounded a small scale to literal zero,
  # which Stan rejects with an opaque error.
  nums <- as.numeric(regmatches(
    fit$info$gp_lscale_prior,
    gregexpr("[0-9]+\\.?[0-9]*(e[-+]?[0-9]+)?", fit$info$gp_lscale_prior))[[1]])
  expect_true(length(nums) > 0)
  expect_true(all(nums > 0))
})

test_that("predict() and fitted() on a bayesian_fit are length-aligned and finite", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan smoke tests")
  skip_if_not_installed("brms")

  pts <- bayes_smoke_points()
  fit <- fit_smoke_model(pts)

  fv <- fitted(fit)
  expect_length(fv, nrow(pts))
  expect_true(all(is.finite(fv)))

  rs <- residuals(fit)
  expect_length(rs, nrow(pts))
  expect_equal(rs, sf::st_drop_geometry(pts)$resp - fv, tolerance = 1e-8)

  # Out-of-sample prediction on held-out rows.
  nd <- bayes_smoke_points(n = 10, seed = 99)
  p <- suppressMessages(predict(fit, newdata = nd))
  expect_length(p, nrow(nd))
  expect_true(all(is.finite(p)))
})

test_that("predict() honours summary and type when newdata is NULL", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan smoke tests")
  skip_if_not_installed("brms")

  # Regression: predict.bayesian_fit() used to short-circuit to the cached
  # fitted() values whenever newdata was NULL, silently ignoring both
  # arguments -- so summary = "median" returned means and type = "predict"
  # returned expected values with no observation noise.
  pts <- bayes_smoke_points()
  fit <- fit_smoke_model(pts)

  mean_epred <- suppressMessages(predict(fit))
  median_epred <- suppressMessages(predict(fit, summary = "median"))
  post_pred <- suppressMessages(predict(fit, type = "predict"))

  expect_length(median_epred, nrow(pts))
  expect_length(post_pred, nrow(pts))
  expect_false(isTRUE(all.equal(mean_epred, median_epred)))
  expect_false(isTRUE(all.equal(mean_epred, post_pred)))

  # draws = TRUE must return a matrix of draws x observations.
  d <- suppressMessages(predict(fit, draws = TRUE))
  expect_true(is.matrix(d))
  expect_identical(ncol(d), nrow(pts))
  expect_gt(nrow(d), 1L)
})

test_that("a supplied control list keeps the package's adapt_delta default", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan smoke tests")
  skip_if_not_installed("brms")

  # Regression: control = list(...) used to replace the defaults wholesale,
  # silently dropping adapt_delta = 0.9 back to Stan's 0.8 -- which is exactly
  # the setting the divergence warning then tells the user to raise.
  pts <- bayes_smoke_points()
  fit <- fit_smoke_model(pts, control = list(max_treedepth = 11))

  # Where the merged control ends up is a Stan-backend detail (rstan and
  # cmdstanr store it differently), so read it defensively and skip rather
  # than fail if this brms/backend combination does not expose it.
  ctl <- tryCatch(fit$engine$fit@stan_args[[1]]$control,
                  error = function(e) NULL)
  if (is.null(ctl$adapt_delta))
    skip("this Stan backend does not expose the merged control list")

  expect_equal(ctl$adapt_delta, 0.9)   # package default survived the merge
  expect_equal(ctl$max_treedepth, 11)  # user value applied
})

test_that("cv_bayes runs a small spatial cross-validation", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan smoke tests")
  skip_if_not_installed("brms")

  pts <- bayes_smoke_points(n = 60)
  folds <- make_folds(pts, k = 2, method = "random_kfold", seed = 7)

  cv <- NULL
  utils::capture.output(
    suppressWarnings(suppressMessages(
      # Fit arguments reach fit_bayesian_spatial_model() via `fit_args`;
      # cv_bayes() has no chains/iter/warmup arguments of its own.
      cv <- cv_bayes(pts, response_var = "resp", predictor_vars = "z",
                     folds = folds, seed = 7,
                     fit_args = list(chains = 1, iter = 200, warmup = 100,
                                     cores = 1, compute_loo = FALSE,
                                     check_convergence = FALSE))
    )),
    type = "output"
  )

  expect_true(is.list(cv))
  expect_identical(cv$n_folds_attempted, 2L)
  expect_gte(cv$n_folds_succeeded, 1L)
  expect_true(all(c("RMSE", "MAE", "R2") %in% names(cv$overall)))
  expect_true(is.finite(cv$overall$RMSE))
  expect_identical(nrow(cv$fold_metrics), cv$n_folds_succeeded)
  # Every observation is predicted exactly once across the folds.
  expect_identical(nrow(cv$predictions), nrow(pts))

  # predictive_coverage is weighted by per-fold n_pred; the coverage entries
  # must be genuine probabilities.
  if (!is.null(cv$predictive_coverage)) {
    cov_vals <- unlist(cv$predictive_coverage[
      grep("^coverage_", names(cv$predictive_coverage))])
    expect_true(all(cov_vals >= 0 & cov_vals <= 1))
  }
})
