# The GP adequacy status of a Bayesian fit is carried where fits are scored:
# compare_models() beside convergence_ok, and cv_bayes() fold_metrics.  A fit
# whose basis was found too coarse ("below_resolution",
# "unresolved_short_range") showed nothing in either table.  The Bayesian
# backend is mocked with helper-lmfit.R's lm fit, so no Stan toolchain is
# needed.

.fat_pts <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               w = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 5 + 2 * d$w + rnorm(n, 0, 0.5)
  d
}

.fat_fit <- function(d, status) {
  f <- lm_spatial_fit(d, "z", "w")
  f$info$convergence_ok <- TRUE
  if (!is.null(status)) f$info$gp_adequacy <- list(status = status)
  f
}

test_that("compare_models() carries the GP adequacy status beside convergence_ok", {
  d <- .fat_pts()
  as_bayes <- function(f) {
    class(f) <- c(class(f)[1L], "bayesian_fit", class(f)[-1L])
    f
  }
  fits <- list(fine    = as_bayes(.fat_fit(d, "ok")),
               coarse  = as_bayes(.fat_fit(d, "below_resolution")),
               missed  = as_bayes(.fat_fit(d, "unresolved_short_range")),
               older   = as_bayes(.fat_fit(d, NULL)),
               plain   = lm_spatial_fit(d, "z", "w"))
  cmp <- suppressWarnings(compare_models(fits))
  expect_true("gp_adequacy" %in% names(cmp))
  expect_type(cmp$gp_adequacy, "character")
  got <- stats::setNames(cmp$gp_adequacy, cmp$model)
  expect_identical(unname(got[c("fine", "coarse", "missed", "older", "plain")]),
                   c("ok", "below_resolution", "unresolved_short_range", NA, NA))
  # Beside convergence_ok, not in place of it.
  expect_true("convergence_ok" %in% names(cmp))
})

test_that("cv_bayes() fold_metrics carries each fold's GP adequacy status", {
  d <- .fat_pts()
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  statuses <- c("ok", "unresolved_short_range", "below_resolution")
  i <- 0L
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, response_var, predictor_vars,
                                          ..., seed = 123) {
      i <<- i + 1L
      .fat_fit(data_sf, statuses[i])
    },
    .package = "spatialkit")
  cv <- suppressWarnings(suppressMessages(cv_bayes(d, "z", "w", folds = f, seed = 1)))
  expect_true("gp_adequacy" %in% names(cv$fold_metrics))
  expect_identical(cv$fold_metrics$gp_adequacy, statuses)
  expect_identical(cv$fold_metrics$convergence_ok, rep(TRUE, 3L))
})

test_that("a fit without the status, and an all-failed run, keep the column", {
  d <- .fat_pts()
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, response_var, predictor_vars,
                                          ..., seed = 123) .fat_fit(data_sf, NULL),
    .package = "spatialkit")
  cv <- suppressWarnings(suppressMessages(cv_bayes(d, "z", "w", folds = f, seed = 1)))
  expect_identical(cv$fold_metrics$gp_adequacy, rep(NA_character_, 3L))

  local_mocked_bindings(
    fit_bayesian_spatial_model = function(...) stop("no fit here"),
    .package = "spatialkit")
  cv0 <- suppressWarnings(suppressMessages(cv_bayes(d, "z", "w", k = 3, seed = 1)))
  expect_true("gp_adequacy" %in% names(cv0$fold_metrics))
  expect_type(cv0$fold_metrics$gp_adequacy, "character")
})
