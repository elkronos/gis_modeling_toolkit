# ===========================================================================
# The spatial_fit S3 surface.
#
# DESCRIPTION sells "a common S3 class ... with consistent predict, fitted,
# residuals and plot methods", but print.spatial_fit(), summary.spatial_fit(),
# print.summary.spatial_fit(), model_metrics.spatial_fit(), new_spatial_fit()
# and clear_fitted_cache() had no coverage at all.  Every one of them is
# reachable through helper-lmfit.R's lm_spatial_fit(), so none of this needs
# an optional backend.
# ===========================================================================

.s3_printed <- function(x, ...) paste(utils::capture.output(print(x, ...)),
                                      collapse = "\n")


# ---------------------------------------------------------------------------
# new_spatial_fit()
# ---------------------------------------------------------------------------

test_that("new_spatial_fit builds the documented structure", {
  pts <- surf_test_points(40)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  expect_s3_class(fit, "spatial_fit")
  expect_equal(class(fit), c("lmsurf_fit", "spatial_fit"))
  expect_equal(fit$n, nrow(pts))
  expect_equal(fit$response_var, "z")
  expect_equal(fit$predictor_vars, "w")
  expect_s3_class(fit$formula, "formula")

  # A cache environment is always provided, with reference semantics so that
  # fitted.bayesian_fit()'s memo survives R's copy-on-modify.
  expect_true(is.environment(fit$info$.cache))
  copy <- fit
  assign(".probe", 42, envir = copy$info$.cache)
  expect_equal(get(".probe", envir = fit$info$.cache), 42)

  # A caller-supplied cache is not replaced.
  own <- new.env(parent = emptyenv())
  custom <- new_spatial_fit("x_fit", engine = NULL, formula = z ~ 1,
                            response_var = "z", predictor_vars = character(0),
                            data_sf = pts, info = list(.cache = own))
  expect_identical(custom$info$.cache, own)

  expect_error(new_spatial_fit(c("a", "b"), NULL, z ~ 1, "z", character(0), pts))
})


# ---------------------------------------------------------------------------
# print.spatial_fit()
# ---------------------------------------------------------------------------

test_that("print.spatial_fit reports the formula and n for any subclass", {
  pts <- surf_test_points(40)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  txt <- .s3_printed(fit)
  expect_match(txt, "<lmsurf_fit> spatial model fit", fixed = TRUE)
  expect_match(txt, "Formula : z ~ w", fixed = TRUE)
  expect_match(txt, "n       : 40", fixed = TRUE)
  # No backend-specific lines for an unknown subclass.
  expect_false(grepl("Bandwidth", txt))
  expect_false(grepl("GP basis", txt))

  # print() returns its argument invisibly.  expect_output() also keeps the
  # cat()ed text out of the suite log.
  expect_output(res <- withVisible(print(fit)))
  expect_false(res$visible)
  expect_identical(res$value, fit)
})


test_that("print.spatial_fit labels and details a gwr_fit", {
  pts <- surf_test_points(30)
  gwr <- new_spatial_fit(
    "gwr_fit", engine = list(), formula = z ~ w, response_var = "z",
    predictor_vars = "w", data_sf = pts,
    info = list(bandwidth = 37.5, adaptive = TRUE, kernel = "gaussian",
                AICc = 412.25, bandwidth_is_fallback = FALSE))

  txt <- .s3_printed(gwr)
  expect_match(txt, "<GWR (GWmodel)> spatial model fit", fixed = TRUE)
  expect_match(txt, "Bandwidth: 37.5 (adaptive, gaussian kernel)", fixed = TRUE)
  expect_match(txt, "AICc    : 412.25", fixed = TRUE)

  # A fixed bandwidth says "fixed", and the kernel defaults to bisquare.
  gwr$info$adaptive <- FALSE
  gwr$info$kernel   <- NULL
  gwr$info$AICc     <- NA_real_
  txt2 <- .s3_printed(gwr)
  expect_match(txt2, "(fixed, bisquare kernel)", fixed = TRUE)
  expect_false(grepl("AICc", txt2))       # a non-finite AICc is not printed
})


test_that("print.spatial_fit labels and details a bayesian_fit", {
  pts <- surf_test_points(30)
  bay <- new_spatial_fit(
    "bayesian_fit", engine = list(), formula = z ~ w, response_var = "z",
    predictor_vars = "w", data_sf = pts,
    info = list(gp_k = 24L, gp_n_basis = 576L, looic = 310.5,
                convergence_ok = TRUE))

  txt <- .s3_printed(bay)
  expect_match(txt, "<Bayesian Spatial GP (brms)> spatial model fit",
               fixed = TRUE)
  expect_match(txt, "GP basis: 24 per dim (576 total)", fixed = TRUE)
  expect_match(txt, "LOOIC   : 310.50", fixed = TRUE)
  expect_false(grepl("Convergence warnings", txt))

  # Convergence trouble is announced, and cannot be missed by omission.
  bay$info$convergence_ok <- FALSE
  expect_match(.s3_printed(bay), "Convergence warnings present", fixed = TRUE)
  bay$info$convergence_ok <- NULL
  expect_match(.s3_printed(bay), "Convergence warnings present", fixed = TRUE)
})


# ---------------------------------------------------------------------------
# summary.spatial_fit() / print.summary.spatial_fit()
# ---------------------------------------------------------------------------

test_that("summary.spatial_fit carries the metrics computed from fitted()", {
  pts <- surf_test_points(80)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")
  s   <- summary(fit)

  expect_s3_class(s, "summary.spatial_fit")
  expect_named(s, c("class", "formula", "n", "response_var", "predictor_vars",
                    "info", "in_sample"))
  expect_equal(s$class, "lmsurf_fit")
  expect_equal(s$n, 80L)
  expect_equal(s$response_var, "z")
  expect_equal(s$predictor_vars, "w")

  # The metrics are .compute_reg_metrics(y, fitted(fit)), hand-checked here.
  y   <- pts$z
  yh  <- fitted(fit)
  expect_equal(s$in_sample$RMSE, sqrt(sum((y - yh)^2) / length(y)))
  expect_equal(s$in_sample$MAE, mean(abs(y - yh)))
  expect_equal(s$in_sample$R2,
               1 - sum((y - yh)^2) / sum((y - mean(y))^2))
  # Adjusted R2 is deliberately suppressed (p = NULL) for every backend.
  expect_true(is.na(s$in_sample$Adj_R2))
  # ... and it agrees with model_metrics() on the same object.
  expect_equal(s$in_sample, model_metrics(fit))
})


test_that("print.summary.spatial_fit distinguishes in-sample from out-of-bag", {
  pts <- surf_test_points(60)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  txt <- .s3_printed(summary(fit))
  expect_match(txt, "Summary of <lmsurf_fit> fit (n = 60)", fixed = TRUE)
  expect_match(txt, "In-sample metrics:", fixed = TRUE)
  expect_match(txt, "RMSE  =")
  expect_match(txt, "MAE   =")
  # The R-squared label carries a superscript 2, which capture.output()
  # renders as a <U+00B2> escape in a C locale -- match around it, not on it.
  expect_match(txt, "\\n    R\\S*\\s*= 0\\.")
  expect_false(grepl("Out-of-bag", txt))

  # The mislabelling that was fixed: an rf_fit's fitted() values are OUT OF
  # BAG, so calling them "In-sample" understates a random hold-out as a
  # memorised fit -- the opposite of the truth, and of the caveat in
  # ?fit_rf_model.  $info$fitted_are_oob drives the label.
  oob <- fit
  oob$info$fitted_are_oob <- TRUE
  txt_oob <- .s3_printed(summary(oob))
  expect_match(txt_oob, "Out-of-bag metrics (NOT in-sample; see ?fit_rf_model)",
               fixed = TRUE)
  expect_false(grepl("In-sample metrics", txt_oob, fixed = TRUE))

  expect_output(res <- withVisible(print(summary(fit))))
  expect_false(res$visible)
})


test_that("print.summary.spatial_fit labels a real rf_fit's metrics out-of-bag", {
  skip_if_not_installed("ranger")
  set.seed(12)
  n <- 80
  dat <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), a = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  dat$z <- 2 * dat$a + rnorm(n, 0, 0.3)

  fit <- fit_rf_model(dat, "z", "a", num_trees = 100L)
  expect_true(isTRUE(fit$info$fitted_are_oob))

  txt <- .s3_printed(summary(fit))
  expect_match(txt, "Out-of-bag metrics (NOT in-sample", fixed = TRUE)
  expect_false(grepl("In-sample metrics", txt, fixed = TRUE))
})


test_that("print.summary.spatial_fit prints optional metrics only when finite", {
  pts <- surf_test_points(60)
  s   <- summary(lm_spatial_fit(pts, predictor_vars = "w"))

  # Adj R2 is NA by construction, so it must not be printed at all.
  expect_false(grepl("Adj R", .s3_printed(s), fixed = TRUE))
  s$in_sample$Adj_R2 <- 0.8123
  expect_match(.s3_printed(s), "Adj R")
  expect_match(.s3_printed(s), "= 0\\.8123")

  s$in_sample$SMAPE <- NA_real_
  expect_false(grepl("SMAPE", .s3_printed(s), fixed = TRUE))
})


# ---------------------------------------------------------------------------
# model_metrics.spatial_fit()
# ---------------------------------------------------------------------------

test_that("model_metrics scores fitted values or newdata predictions", {
  train <- surf_test_points(90, seed = 1)
  test  <- surf_test_points(40, seed = 77)
  fit   <- lm_spatial_fit(train, predictor_vars = "w")

  ins <- model_metrics(fit)
  expect_s3_class(ins, "data.frame")
  expect_named(ins, c("n", "RMSE", "MAE", "MAPE", "SMAPE", "R2", "Adj_R2"))
  expect_equal(ins$n, 90L)

  oos <- model_metrics(fit, newdata = test)
  expect_equal(oos$n, 40L)
  # The out-of-sample numbers really come from predict() on newdata.
  yh <- predict(fit, newdata = test)
  expect_equal(oos$RMSE, sqrt(sum((test$z - yh)^2) / nrow(test)))
  expect_equal(oos$R2, 1 - sum((test$z - yh)^2) /
                 sum((test$z - mean(test$z))^2))
  expect_false(isTRUE(all.equal(ins$RMSE, oos$RMSE)))

  # newdata without the response cannot be scored, and says why.
  no_y <- test
  no_y$z <- NULL
  expect_error(model_metrics(fit, newdata = no_y),
               "must contain the response variable 'z'")
  expect_error(model_metrics(fit, newdata = no_y), "predict\\(\\) can be used")
  # ... while predict() itself is still fine on it.
  expect_length(predict(fit, newdata = no_y), nrow(no_y))
})


# ---------------------------------------------------------------------------
# clear_fitted_cache()
# ---------------------------------------------------------------------------

test_that("clear_fitted_cache drops the memoised fitted values", {
  pts <- surf_test_points(30)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")
  cache <- fit$info$.cache

  assign(".fitted_values", rep(-1, nrow(pts)), envir = cache)
  expect_true(exists(".fitted_values", envir = cache, inherits = FALSE))

  out <- clear_fitted_cache(fit)
  expect_false(exists(".fitted_values", envir = cache, inherits = FALSE))
  expect_identical(out, fit)
  expect_invisible(clear_fitted_cache(fit))

  # Idempotent, and safe on a fit that carries no cache at all.
  expect_silent(clear_fitted_cache(fit))
  no_cache <- fit
  no_cache$info$.cache <- NULL
  expect_silent(clear_fitted_cache(no_cache))
})


# ---------------------------------------------------------------------------
# The coef() contract
# ---------------------------------------------------------------------------

test_that("coef() on a spatial_fit errors rather than returning NULL", {
  # ?new_spatial_fit documents this: a NULL return would be indistinguishable
  # from "this model genuinely has no fixed effects", so lapply(fits, coef)
  # would quietly come back shorter than the caller expected.  gwr_fit and
  # bayesian_fit used to return NULL; all three now stop().
  pts <- surf_test_points(30)

  gwr <- new_spatial_fit("gwr_fit", engine = list(), formula = z ~ w,
                         response_var = "z", predictor_vars = "w",
                         data_sf = pts, info = list())
  expect_error(coef(gwr), "no GWmodel `SDF` component")

  bay <- new_spatial_fit("bayesian_fit", engine = list(), formula = z ~ w,
                         response_var = "z", predictor_vars = "w",
                         data_sf = pts, info = list())
  expect_error(coef(bay), "coef\\.bayesian_fit\\(\\)")

  # A gwr_fit whose engine DOES carry an SDF returns the local coefficients.
  sdf <- sf::st_sf(Intercept = rnorm(30), w = rnorm(30),
                   geometry = sf::st_geometry(pts))
  gwr_ok <- gwr
  gwr_ok$engine <- list(SDF = sdf)
  co <- coef(gwr_ok)
  expect_s3_class(co, "data.frame")
  expect_equal(nrow(co), 30L)
  expect_false(inherits(co, "sf"))          # geometry dropped
  expect_named(co, c("Intercept", "w"))
})


test_that("coef() on an rf_fit points at importance instead", {
  skip_if_not_installed("ranger")
  set.seed(13)
  n <- 60
  dat <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  dat$z <- 2 * dat$a + rnorm(n, 0, 0.3)

  fit <- fit_rf_model(dat, "z", c("a", "b"), num_trees = 60L)
  expect_error(coef(fit), "no coefficients")
  expect_error(coef(fit), "importance")
})
