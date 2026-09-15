# tests/testthat/test-cv-metrics.R
# ---------------------------------------------------------------------------
# The `metrics` hook of the cv_*() functions: a user scoring function applied
# per fold and to the pooled predictions, whose names become columns beside
# the built-in Gaussian set.
# ---------------------------------------------------------------------------

cvm_site <- function(n = 120, seed = 1) {
  set.seed(seed)
  site <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  site$price <- 10 + 0.01 * sf::st_coordinates(site)[, 1] + 2 * site$elev + rnorm(n)
  site
}

cvm_lm_fit <- function(train_sf) {
  new_spatial_fit(
    subclass = "cvm_lm", engine = stats::lm(price ~ elev, sf::st_drop_geometry(train_sf)),
    formula = price ~ elev, response_var = "price", predictor_vars = "elev",
    data_sf = train_sf)
}
predict.cvm_lm <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data_sf
  as.numeric(stats::predict(object$engine, sf::st_drop_geometry(newdata)))
}
registerS3method("predict", "cvm_lm", predict.cvm_lm)

med_ae <- function(y, yhat) c(MedAE = stats::median(abs(y - yhat)), bias = mean(yhat - y))

run_cv <- function(..., metrics = med_ae, k = 3, seed = 1) {
  suppressMessages(cv_spatial(cvm_site(), "price", "elev", fit_fn = cvm_lm_fit,
                              k = k, seed = seed, metrics = metrics, ...))
}


test_that("a user metric is applied per fold and to the pooled predictions", {
  cv <- run_cv()
  expect_true(all(c("MedAE", "bias") %in% names(cv$overall)))
  expect_true(all(c("MedAE", "bias") %in% names(cv$fold_metrics)))
  # The built-in columns come first and are unchanged in position.
  expect_equal(names(cv$overall)[1:9],
               c("RMSE", "MAE", "MAPE", "SMAPE", "R2", "Adj_R2", "n_pred", "n_MAPE", "n_SMAPE"))

  p  <- cv$predictions
  ok <- is.finite(p$y) & is.finite(p$yhat)
  expect_equal(cv$overall$MedAE, stats::median(abs(p$y[ok] - p$yhat[ok])))
  expect_equal(cv$overall$bias, mean(p$yhat[ok] - p$y[ok]))
  for (f in cv$fold_metrics$fold) {
    q <- p[p$fold == f, ]
    expect_equal(cv$fold_metrics$MedAE[cv$fold_metrics$fold == f],
                 stats::median(abs(q$y - q$yhat)))
  }
  # Nothing else about the result moved.
  ref <- run_cv(metrics = NULL)
  expect_equal(cv$overall[, 1:9], ref$overall)
  expect_equal(cv$predictions, ref$predictions)
})

test_that("parallel and sequential runs agree on the user columns", {
  skip_on_os("windows")
  cv  <- run_cv()
  cvp <- run_cv(parallel = 2)
  expect_identical(cvp$overall, cv$overall)
  expect_identical(cvp$fold_metrics, cv$fold_metrics)
})

test_that("a named list and a one-row data frame are accepted as return values", {
  cv_l <- run_cv(metrics = function(y, yhat) list(q90 = unname(stats::quantile(abs(y - yhat), 0.9))))
  cv_d <- run_cv(metrics = function(y, yhat) data.frame(q90 = unname(stats::quantile(abs(y - yhat), 0.9))))
  expect_true("q90" %in% names(cv_l$overall))
  expect_identical(cv_l$overall$q90, cv_d$overall$q90)
  expect_identical(cv_l$fold_metrics$q90, cv_d$fold_metrics$q90)
})

test_that("the argument and the return value are validated", {
  site <- cvm_site()
  expect_error(cv_spatial(site, "price", "elev", fit_fn = cvm_lm_fit, k = 3, metrics = 3),
               "`metrics` must be a function of two arguments")
  expect_error(cv_spatial(site, "price", "elev", fit_fn = cvm_lm_fit, k = 3,
                          metrics = function(y) 1),
               "fewer than two arguments")
  # Reserved names, unnamed values, duplicated names and non-scalars are
  # errors at the first application, not silently reshaped.
  expect_error(run_cv(metrics = function(y, yhat) c(RMSE = 1)),
               "already columns of the metrics frames: RMSE")
  expect_error(run_cv(metrics = function(y, yhat) 1),
               "every element is named")
  expect_error(run_cv(metrics = function(y, yhat) c(a = 1, a = 2)),
               "duplicated names: a")
  expect_error(run_cv(metrics = function(y, yhat) list(a = 1, b = "x")),
               "element 'b' is not a single numeric value")
  expect_error(run_cv(metrics = function(y, yhat) data.frame(a = 1:2)),
               "data frame with 2 rows")
})

test_that("a metric that throws gives NA columns, not a dropped fold", {
  thrower <- function(y, yhat) if (length(y) > 30) stop("boom") else c(z = 1)
  lines <- capture_spatialkit_log(cv <- run_cv(metrics = thrower))
  expect_equal(cv$n_folds_succeeded, 3L)
  expect_true("z" %in% names(cv$fold_metrics))
  expect_true(all(is.na(cv$fold_metrics$z)))
  expect_true(is.na(cv$overall$z))
  expect_true(log_has(lines, "`metrics` failed on fold"))
  expect_true(log_has(lines, "`metrics` failed on the pooled predictions"))
  # Non-finite values come back as NA rather than Inf/NaN.
  cv2 <- run_cv(metrics = function(y, yhat) c(inf = Inf, nan = NaN))
  expect_true(is.na(cv2$overall$inf)); expect_true(is.na(cv2$overall$nan))
})

test_that("the empty frames of a run where every fold failed carry the metric columns", {
  bad_fit <- function(train_sf) stop("nope")
  cv <- suppressWarnings(suppressMessages(
    cv_spatial(cvm_site(), "price", "elev", fit_fn = bad_fit, k = 3, seed = 1,
               metrics = med_ae)))
  expect_equal(nrow(cv$fold_metrics), 0L)
  expect_true(all(c("MedAE", "bias") %in% names(cv$fold_metrics)))
  expect_type(cv$fold_metrics$MedAE, "double")
  expect_true(all(c("MedAE", "bias") %in% names(cv$overall)))
  expect_true(is.na(cv$overall$MedAE))
  expect_equal(cv$overall$n_pred, 0L)
})

test_that("cv_rf() and compare_models_cv() carry the metric through, and the latter protects it", {
  skip_if_not_installed("ranger")
  site <- cvm_site()
  r <- suppressMessages(cv_rf(site, "price", "elev", k = 3, seed = 1, num_trees = 50,
                              metrics = med_ae))
  expect_true(all(c("MedAE", "bias") %in% names(r$overall)))
  cmp <- suppressMessages(suppressWarnings(
    compare_models_cv(site, "price", "elev", models = "RF", k = 3,
                      rf_args = list(num_trees = 50), metrics = med_ae, quiet = TRUE)))
  expect_true(all(c("MedAE", "bias") %in% names(cmp$overall)))
  expect_equal(names(cmp$overall)[ncol(cmp$overall)], "model")
  expect_true(all(c("MedAE", "bias") %in% names(cmp$by_fold)))
  expect_equal(cmp$overall$MedAE, cmp$rf_cv$overall$MedAE)
  # A per-model metrics entry is ignored with a warning: one scoring
  # function for every row of `overall`.
  expect_warning(
    suppressMessages(compare_models_cv(site, "price", "elev", models = "RF", k = 3,
                                       rf_args = list(num_trees = 50, metrics = med_ae),
                                       quiet = TRUE)),
    "ignoring 'metrics' in `rf_args`")
  expect_error(compare_models_cv(site, "price", "elev", models = "RF", metrics = "x"),
               "compare_models_cv\\(\\): `metrics` must be a function")
})

test_that("cv_gwr() and cv_bayes() accept and validate the argument", {
  site <- cvm_site()
  expect_error(cv_gwr(site, "price", "elev", k = 3, metrics = "x"),
               "cv_gwr\\(\\): `metrics` must be a function")
  expect_error(cv_bayes(site, "price", "elev", k = 3, metrics = "x"),
               "cv_bayes\\(\\): `metrics` must be a function")
})
