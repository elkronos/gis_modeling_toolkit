# tests/testthat/test-followups-convergence-crs.R
# ---------------------------------------------------------------------------
# A Bayesian fit whose sampler did not converge is flagged where fits are
# scored (cv_bayes(), compare_models()), not only in the log; and
# predict_surface() aligns a CRS-less `grid` or `covariates` with an R
# warning, as it does `boundary`.  The Bayesian backend is mocked with
# helper-lmfit.R's lm fit, so no Stan toolchain is needed.
# ---------------------------------------------------------------------------

.fu_pts <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               w = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 5 + 2 * d$w + rnorm(n, 0, 0.5)
  d
}

# The lm stand-in, carrying the convergence verdict fit_bayesian_spatial_model()
# records in $info$convergence_ok.
.fu_fit <- function(d, ok) {
  f <- lm_spatial_fit(d, "z", "w")
  f$info$convergence_ok <- ok
  f
}

test_that("cv_bayes() marks folds whose sampler did not converge and warns once", {
  d <- .fu_pts()
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  verdicts <- c(TRUE, FALSE, TRUE)
  i <- 0L
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, response_var, predictor_vars,
                                          ..., seed = 123) {
      i <<- i + 1L
      .fu_fit(data_sf, verdicts[i])
    },
    .package = "spatialkit")

  ws <- character(0)
  cv <- withCallingHandlers(
    suppressMessages(cv_bayes(d, "z", "w", folds = f, seed = 1)),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_true("convergence_ok" %in% names(cv$fold_metrics))
  expect_identical(cv$fold_metrics$convergence_ok, verdicts)
  conv <- grep("did not converge", ws, value = TRUE)
  expect_length(conv, 1L)
  expect_match(conv, "in 1 of 3 fold\\(s\\) \\(fold 2\\)")
  expect_match(conv, "fold_metrics$convergence_ok", fixed = TRUE)
})

test_that("cv_bayes() stays quiet when every fold converged or none was checked", {
  d <- .fu_pts()
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  for (ok in list(TRUE, NA)) {
    local_mocked_bindings(
      fit_bayesian_spatial_model = function(data_sf, response_var, predictor_vars,
                                            ..., seed = 123) .fu_fit(data_sf, ok),
      .package = "spatialkit")
    ws <- character(0)
    cv <- withCallingHandlers(
      suppressMessages(cv_bayes(d, "z", "w", folds = f, seed = 1)),
      warning = function(w) {
        ws <<- c(ws, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_identical(cv$fold_metrics$convergence_ok, rep(ok, 3L))
    expect_false(any(grepl("did not converge", ws)))
  }
})

test_that("an all-failed cv_bayes() run still has the convergence_ok column", {
  d <- .fu_pts()
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(...) stop("no fit here"),
    .package = "spatialkit")
  cv <- suppressWarnings(suppressMessages(cv_bayes(d, "z", "w", k = 3, seed = 1)))
  expect_true("convergence_ok" %in% names(cv$fold_metrics))
  expect_type(cv$fold_metrics$convergence_ok, "logical")
})

test_that("compare_models() carries convergence_ok and warns about a non-converged fit", {
  d <- .fu_pts()
  as_bayes <- function(f) {
    class(f) <- c(class(f)[1L], "bayesian_fit", class(f)[-1L])
    f
  }
  fits <- list(ok      = as_bayes(.fu_fit(d, TRUE)),
               bad     = as_bayes(.fu_fit(d, FALSE)),
               unknown = as_bayes(.fu_fit(d, NA)),
               plain   = lm_spatial_fit(d, "z", "w"))
  ws <- character(0)
  cmp <- withCallingHandlers(
    compare_models(fits),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  got <- stats::setNames(cmp$convergence_ok, cmp$model)
  expect_identical(unname(got[c("ok", "bad", "unknown", "plain")]),
                   c(TRUE, FALSE, NA, NA))
  conv <- grep("did not converge", ws, value = TRUE)
  expect_length(conv, 1L)
  expect_match(conv, "'bad'", fixed = TRUE)
})

test_that("predict_surface() warns when it stamps or reprojects a CRS-less grid or covariates", {
  d <- .fu_pts()
  fit <- lm_spatial_fit(d, "z", "w")

  # Projected-looking coordinates without a CRS: stamped with the fit's CRS.
  grid <- sf::st_as_sf(data.frame(x = 5e5 + c(100, 500, 900),
                                  y = 5e6 + c(100, 500, 900), w = 0),
                       coords = c("x", "y"))
  expect_warning(s <- predict_surface(fit, grid = grid),
                 "predict_surface\\(\\): `grid` has no CRS.*stamping")
  expect_equal(sf::st_crs(s), sf::st_crs(d))
  expect_equal(nrow(s), 3L)

  # The same for the covariates layer, named as such.
  g2  <- sf::st_as_sf(data.frame(x = 5e5 + c(200, 700), y = 5e6 + c(300, 800)),
                      coords = c("x", "y"), crs = 32632)
  cov <- sf::st_as_sf(data.frame(x = 5e5 + c(200, 700), y = 5e6 + c(300, 800),
                                 w = c(1, 2)),
                      coords = c("x", "y"))
  expect_warning(s2 <- predict_surface(fit, grid = g2, covariates = cov),
                 "predict_surface\\(\\): `covariates` has no CRS.*stamping")
  expect_equal(s2$w, c(1, 2))

  # Lon/lat-looking coordinates are reprojected from EPSG:4326, with a warning.
  ll <- sf::st_coordinates(sf::st_transform(g2, 4326))
  g3 <- sf::st_as_sf(data.frame(x = ll[, 1], y = ll[, 2], w = 0),
                     coords = c("x", "y"))
  expect_warning(s3 <- predict_surface(fit, grid = g3),
                 "predict_surface\\(\\): `grid` has no CRS; its coordinates look like lon/lat")
  expect_equal(unname(sf::st_coordinates(s3)), unname(sf::st_coordinates(g2)),
               tolerance = 1e-3)

  # A grid that has a CRS is used without any such warning.
  ws <- character(0)
  withCallingHandlers(predict_surface(fit, grid = sf::st_set_crs(grid, 32632)),
                      warning = function(w) {
                        ws <<- c(ws, conditionMessage(w))
                        invokeRestart("muffleWarning")
                      })
  expect_false(any(grepl("has no CRS", ws)))
})
