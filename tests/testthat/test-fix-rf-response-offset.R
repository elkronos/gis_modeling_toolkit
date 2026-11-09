# fit_rf_model() handed ranger the response as given.  ranger scores a split
# by sums of squares it forms from the response, so at an offset of 1e8 or
# 1e12 those scores read the rounding of the offset and the trees split
# elsewhere: fitted values up to 1.3 from those at the origin.  At 1e160 the
# squares overflow and at 1e-160 the OOB error loses its bits.  The forest is
# now grown on the response taken about its centre where that centre is more
# than .PRED_CENTRE_RATIO times its spread, and divided by a power of two
# where it lies beyond .RF_Y_MAX or its spread below .RF_Y_MIN; the centre
# and scale are put back on everything the fit returns.

skip_if_not_installed("ranger")

rfo_layer <- function(n = 80, seed = 3) {
  set.seed(seed)
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), yy = runif(n, 0, 1000),
                               a = rnorm(n), b = rnorm(n)),
                    coords = c("x", "yy"), crs = 32632)
  d$z <- 2 * d$a - d$b + rnorm(n, 0, 0.3)
  d
}

rfo_fit <- function(d, ...)
  suppressWarnings(fit_rf_model(d, "z", c("a", "b"), num_trees = 50L,
                                num_threads = 1L, seed = 1L, ...))

# The split variables and values of every tree; the predictions of the
# terminal nodes are in the response's units and left out.
rfo_trees <- function(f) {
  lapply(seq_len(f$engine$num.trees), function(k) {
    ti <- ranger::treeInfo(f$engine, k)
    ti[, c("nodeID", "leftChild", "rightChild", "splitvarName", "splitval",
           "terminal")]
  })
}

# One unit in the last place of x.
rfo_ulp <- function(x) 2^floor(log2(abs(x))) * .Machine$double.eps

test_that("fit_rf_model() grows the same forest on a response at an offset", {
  d0 <- rfo_layer()
  f0 <- rfo_fit(d0)
  s  <- stats::sd(d0$z)
  expect_identical(f0$info$response_centre, 0)
  expect_identical(f0$info$response_scale, 1)
  for (off in c(1e8, 1e12)) {
    d <- d0
    d$z <- d0$z + off
    f <- rfo_fit(d)
    # The rounding of the stored response itself is the floor.
    tol <- max(1e-10 * s, 4 * rfo_ulp(off))
    expect_true(abs(f$info$response_centre - off) < 10 * s)
    expect_identical(f$info$response_scale, 1)
    expect_identical(rfo_trees(f), rfo_trees(f0))
    expect_lt(max(abs((fitted(f) - off) - fitted(f0))), tol)
    expect_lt(max(abs((predict(f, newdata = d) - off) -
                        predict(f0, newdata = d0))), tol)
    expect_lt(max(abs(residuals(f) - residuals(f0))), tol)
    expect_lt(abs(f$info$oob_rmse - f0$info$oob_rmse), tol)
    expect_equal(f$info$oob_r_squared, f0$info$oob_r_squared, tolerance = 1e-6)
    # Squared units: the rounding of the response moves the OOB errors the
    # importance is a difference of.
    expect_lt(max(abs(f$info$importance - f0$info$importance)), 4 * s * tol)
  }
})

test_that("fit_rf_model() grows the same forest on a response whose squares overflow or underflow", {
  d0 <- rfo_layer()
  f0 <- rfo_fit(d0)
  for (k in c(-530, 530)) {
    by <- 2^k
    d <- d0
    d$z <- d0$z * by
    f <- rfo_fit(d)
    expect_identical(f$info$response_centre, 0)
    expect_false(f$info$response_scale == 1)
    expect_identical(rfo_trees(f), rfo_trees(f0))
    # A power of two scales exactly, so the forest is the base's to the bit.
    expect_identical(fitted(f) / by, fitted(f0))
    expect_identical(predict(f, newdata = d) / by, predict(f0, newdata = d0))
    expect_identical(residuals(f) / by, residuals(f0))
    expect_identical(f$info$oob_rmse / by, f0$info$oob_rmse)
    expect_identical(f$info$oob_r_squared, f0$info$oob_r_squared)
  }
})

test_that("an importance beyond the range of doubles is said, not left silent", {
  d0 <- rfo_layer()
  d <- d0
  d$z <- d0$z * 2^530
  lines <- capture_spatialkit_log(f <- rfo_fit(d), level = logger::WARN)
  expect_true(all(is.infinite(f$info$importance)))
  expect_true(log_has(lines, "importance of .* lies above the largest double"))
  # The OOB RMSE, in the response's units, is held.
  expect_true(is.finite(f$info$oob_rmse))
  # At 2^-530 the importance is near 1e-320, which keeps few of its bits.
  d$z <- d0$z * 2^-530
  lines <- capture_spatialkit_log(f <- rfo_fit(d), level = logger::WARN)
  expect_true(log_has(lines, "below the smallest normal double"))
})

test_that("cv_rf() scores the same forests on a response at an offset", {
  d0 <- rfo_layer()
  folds <- rep_len(1:4, nrow(d0))
  cv <- function(d) suppressWarnings(
    cv_rf(d, "z", c("a", "b"), folds = folds, num_trees = 50L,
          num_threads = 1L, seed = 1L))
  c0 <- cv(d0)
  key <- function(r) r$predictions[order(r$predictions$..row_id), ]
  for (off in c(1e8, 1e12)) {
    d <- d0
    d$z <- d0$z + off
    c1 <- cv(d)
    tol <- max(1e-10 * stats::sd(d0$z), 4 * rfo_ulp(off))
    expect_lt(max(abs((key(c1)$yhat - off) - key(c0)$yhat)), tol)
    expect_lt(abs(c1$overall$RMSE - c0$overall$RMSE), tol)
    expect_lt(max(abs(c1$fold_metrics$RMSE - c0$fold_metrics$RMSE)), tol)
  }
  d <- d0
  d$z <- d0$z * 2^530
  lines <- capture_spatialkit_log(c1 <- cv(d), level = logger::WARN)
  expect_identical(key(c1)$yhat / 2^530, key(c0)$yhat)
  # The fold forests' importance is never read, so nothing is said of it.
  expect_false(log_has(lines, "importance"))
})

test_that("a response of ordinary size is handed to ranger as it is", {
  d0 <- rfo_layer()
  for (shift in c(0, 100)) {
    d <- d0
    d$z <- d0$z + shift
    f <- rfo_fit(d)
    expect_identical(f$info$response_centre, 0)
    expect_identical(f$info$response_scale, 1)
    X <- sf::st_drop_geometry(f$data_sf)[, c("a", "b")]
    r <- ranger::ranger(x = X, y = f$data_sf$z, num.trees = 50L,
                        importance = "permutation", seed = 1L,
                        replace = TRUE, sample.fraction = 1, num.threads = 1L)
    expect_identical(f$engine$forest, r$forest)
    expect_identical(f$engine$predictions, r$predictions)
    expect_identical(fitted(f), r$predictions)
    expect_identical(f$info$importance, r$variable.importance)
    expect_identical(f$info$oob_rmse, sqrt(r$prediction.error))
  }
})

test_that("splitrule = 'beta' keeps the response in (0, 1) as given", {
  d0 <- rfo_layer()
  d <- d0
  # Centre 0.5, spread about 1e-5: far from zero against its spread.
  d$z <- 0.5 + 1e-5 * d0$z / max(abs(d0$z))
  f <- rfo_fit(d, splitrule = "beta")
  expect_identical(f$info$response_centre, 0)
  expect_identical(f$info$response_scale, 1)
  expect_true(all(fitted(f) > 0 & fitted(f) < 1, na.rm = TRUE))
})

test_that("an importance ranger reports without units is not scaled", {
  d0 <- rfo_layer()
  # Permutation importance divided by its standard error, and impurity
  # importance under maxstat (a rank statistic): neither moves with the
  # units of the response, so a power of two gives the base's importance
  # (to the rounding of the standard error, which ranger forms on the
  # response it is given, a power of two times the base's).
  cases <- list(
    list(importance = "permutation", scale.permutation.importance = TRUE),
    list(importance = "impurity", splitrule = "maxstat"))
  for (args in cases) {
    f0 <- do.call(rfo_fit, c(list(d0), args))
    for (k in c(-400, 400)) {
      d <- d0
      d$z <- d0$z * 2^k
      lines <- capture_spatialkit_log(f <- do.call(rfo_fit, c(list(d), args)),
                                      level = logger::WARN)
      expect_false(f$info$response_scale == 1)
      expect_equal(f$info$importance, f0$info$importance, tolerance = 1e-12)
      expect_false(log_has(lines, "importance"))
    }
  }
  # Impurity importance under the default splitrule is in squared units.
  f0 <- rfo_fit(d0, importance = "impurity")
  d <- d0
  d$z <- d0$z * 2^-400
  f <- rfo_fit(d, importance = "impurity")
  expect_identical(f$info$importance / 2^-400 / 2^-400, f0$info$importance)
})

test_that("a response taken about a centre and scaled is put back in that order", {
  d0 <- rfo_layer()
  # On a grid of 2^-13, so 2^-480 + z * 2^-500 holds z exactly.
  d0$z <- round(d0$z * 2^13) / 2^13
  f0 <- rfo_fit(d0)
  d <- d0
  d$z <- 2^-480 + d0$z * 2^-500
  expect_identical((d$z - 2^-480) / 2^-500, d0$z)
  f <- rfo_fit(d)
  expect_false(f$info$response_centre == 0)
  expect_false(f$info$response_scale == 1)
  expect_identical(rfo_trees(f), rfo_trees(f0))
  # The rounding of the centre, in units of z.
  tol <- 16 * rfo_ulp(2^-480) / 2^-500
  expect_lt(max(abs((fitted(f) - 2^-480) / 2^-500 - fitted(f0))), tol)
  expect_lt(max(abs((predict(f, newdata = d) - 2^-480) / 2^-500 -
                      predict(f0, newdata = d0))), tol)
  expect_lt(max(abs(residuals(f) / 2^-500 - residuals(f0))), tol)
  expect_lt(abs(f$info$oob_rmse / 2^-500 - f0$info$oob_rmse), tol)
})

test_that("print() shows an importance above the largest double as Inf", {
  d0 <- rfo_layer()
  d <- d0
  d$z <- d0$z * 2^530
  f <- rfo_fit(d)
  txt <- paste(utils::capture.output(print(f)), collapse = "\n")
  expect_match(txt, "Importance (permutation): ", fixed = TRUE)
  expect_match(txt, "a=Inf", fixed = TRUE)
  expect_match(txt, "b=Inf", fixed = TRUE)
  expect_match(txt, "(Inf: above the largest double in squared units of the response)",
               fixed = TRUE)
  expect_no_match(txt, "undefined", fixed = TRUE)
})
