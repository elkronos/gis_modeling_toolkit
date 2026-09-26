# tests/testthat/test-review2-eval.R
# ---------------------------------------------------------------------------
# Regressions from the second review of the metrics, compare_models(),
# compare_models_cv(), residual_morans_i() and forward selection.  Each test
# names the finding it closes.
# ---------------------------------------------------------------------------

.r2e_pts <- function(n = 90, seed = 1, extent = 1000) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = runif(n, 0, extent), y = runif(n, 0, extent), w = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 5 + 2 * d$w + rnorm(n, 0, 0.5)
  d
}

.r2e_lm <- function(tr) lm_spatial_fit(tr, "z", "w")

.r2e_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}


# ---------------------------------------------------------------------------
# evaluation MED-LOW-3: one R2 baseline for newdata and cross-validation
# ---------------------------------------------------------------------------

test_that("model_metrics(newdata =) scores R2 against the training mean, as CV does", {
  d <- .r2e_pts(200, seed = 11)
  x <- sf::st_coordinates(d)[, 1]
  d$z <- 0.01 * x + d$w + rnorm(200, 0, 0.3)      # a trend w cannot explain
  tr <- which(x < 600); te <- which(x >= 600)
  fit <- lm_spatial_fit(d[tr, ], "z", "w")
  mm  <- model_metrics(fit, newdata = d[te, ])
  yh  <- predict(fit, newdata = d[te, ])
  expect_equal(mm$R2, 1 - sum((d$z[te] - yh)^2) / sum((d$z[te] - mean(d$z[tr]))^2))
  # The same predictions through cv_spatial() on the same single split.
  cv <- cv_spatial(d, "z", "w", fit_fn = .r2e_lm,
                   folds = list(list(train = tr, test = te)))
  expect_equal(cv$predictions$yhat, unname(yh))
  expect_equal(mm$R2, cv$overall$R2)
  # In sample the two means coincide: the ordinary R2.
  ins <- model_metrics(fit)
  expect_equal(ins$R2, summary(fit$engine)$r.squared)
})


# ---------------------------------------------------------------------------
# evaluation LOW-4 / aoa-utils L6: scale-relative tolerances
# ---------------------------------------------------------------------------

test_that("R2 and MAPE do not depend on the units of the response", {
  set.seed(12)
  y <- rnorm(100, 5, 1); yhat <- y + rnorm(100, 0, 0.3)
  ref <- spatialkit:::.compute_reg_metrics(y, yhat)
  for (s in c(1e-9, 1e-15, 1e6)) {
    m <- spatialkit:::.compute_reg_metrics(y * s, yhat * s)
    expect_equal(m$R2, ref$R2, info = format(s))
    expect_equal(m$MAPE, ref$MAPE, info = format(s))
    expect_identical(m$n_MAPE, 100L, info = format(s))
  }
  # A large offset with an ordinary spread keeps its R2 ...
  big <- spatialkit:::.compute_reg_metrics(1e8 + y / 10, 1e8 + yhat / 10)
  expect_equal(big$R2, ref$R2, tolerance = 1e-6)
  # ... and a constant response still has none.
  expect_true(is.na(spatialkit:::.compute_reg_metrics(rep(5e-9, 4),
                                                      5e-9 * c(1, 1.1, 0.9, 1))$R2))
})

test_that("select_features_forward(metric = 'R2') works on a response in small units", {
  d <- .r2e_pts(120, seed = 13)
  d$a <- rnorm(120); d$b <- rnorm(120)
  d$z <- (3 * d$a + rnorm(120, 0, 0.5)) * 1e-9
  fs <- suppressWarnings(select_features_forward(
    d, "z", c("a", "b"), fit_fn = function(tr, v) lm_spatial_fit(tr, "z", v),
    k = 3, metric = "R2", quiet = TRUE))
  expect_true("a" %in% fs$selected)
  expect_true(is.finite(fs$score))
})


# ---------------------------------------------------------------------------
# evaluation LOW-5: evaluate_insample() and compare_models() label the basis
# ---------------------------------------------------------------------------

test_that("evaluate_insample() says whether each row is in-sample, out-of-bag or newdata", {
  d <- .r2e_pts(60, seed = 14)
  f_in  <- lm_spatial_fit(d, "z", "w")
  f_oob <- f_in; f_oob$info$fitted_are_oob <- TRUE
  out <- evaluate_insample(list(lm = f_in, forest = f_oob))
  expect_identical(out$metric_basis, c("in-sample", "out-of-bag"))
  expect_identical(evaluate_insample(list(lm = f_in), newdata = d[1:20, ])$metric_basis,
                   "newdata")
  lines <- capture_spatialkit_log(
    cmp <- suppressWarnings(compare_models(list(lm = f_in, forest = f_oob))))
  expect_identical(cmp$metric_basis[match(c("lm", "forest"), cmp$model)],
                   c("in-sample", "out-of-bag"))
  expect_true(log_has(lines, "the metrics mix bases"))
})


# ---------------------------------------------------------------------------
# gaps G1.2: information criteria compared only on the same rows
# ---------------------------------------------------------------------------

test_that("compare_models() blanks LOOIC for fits on different rows, with a warning", {
  d <- .r2e_pts(70, seed = 15)
  stub <- function(data, looic) {
    fit <- lm_spatial_fit(data, "z", "w")
    fit$info$looic <- looic
    class(fit) <- c("lmsurf_fit", "bayesian_fit", "spatial_fit")
    fit
  }
  a <- stub(d, 54.8); b <- stub(d[1:50, ], 32.2)
  res <- .r2e_warnings(compare_models(list(A = a, B = b)))
  expect_true(any(grepl("LOOIC is a sum over the rows.*A: n = 70, B: n = 50",
                        res$warnings)))
  expect_true(all(is.na(res$value$LOOIC)))
  # The same rows, in another order: comparable, kept, no such warning.
  b2 <- stub(d[70:1, ], 55.1)
  res2 <- .r2e_warnings(compare_models(list(A = a, B = b2)))
  expect_false(any(grepl("LOOIC is a sum", res2$warnings)))
  expect_equal(res2$value$LOOIC[match(c("A", "B"), res2$value$model)], c(54.8, 55.1))
})


# ---------------------------------------------------------------------------
# gaps G1.10: negative residual autocorrelation is not "missed structure"
# ---------------------------------------------------------------------------

test_that("compare_models() reads significantly negative Moran's I as over-fitting", {
  set.seed(5)
  gx <- rep(1:10, each = 10); gy <- rep(1:10, 10)
  pts <- sf::st_as_sf(data.frame(x = gx * 100 + runif(100, -5, 5),
                                 y = gy * 100 + runif(100, -5, 5), w = rnorm(100)),
                      coords = c("x", "y"), crs = 32632)
  pts$z <- rnorm(100)
  fit <- lm_spatial_fit(pts, "z", "w")
  # Residuals alternating in sign column by column: anti-correlated.
  registerS3method("residuals", "r2_altresid", function(object, ...) {
    xy <- sf::st_coordinates(object$data_sf)
    (-1)^round(xy[, 1] / 100) + stats::rnorm(nrow(xy), 0, 0.1)
  })
  class(fit) <- c("r2_altresid", class(fit))
  lines <- capture_spatialkit_log(cmp <- compare_models(list(alt = fit)))
  expect_lt(cmp$resid_morans_z, 0)
  expect_lt(cmp$resid_morans_p, 0.05)
  expect_true(log_has(lines, "significant negative spatial autocorrelation.*over-fitting"))
  expect_false(log_has(lines, "may not fully capture"))
})


# ---------------------------------------------------------------------------
# gaps G6.8: residual_morans_i() says why it has no residuals
# ---------------------------------------------------------------------------

test_that("residual_morans_i() reports why the residuals could not be used", {
  d <- .r2e_pts(100, seed = 16)
  fit <- new_spatial_fit(subclass = "r2_nomethod", engine = NULL,
                         formula = z ~ w, response_var = "z",
                         predictor_vars = "w", data_sf = d)
  expect_warning(out <- residual_morans_i(fit), "has no residuals\\(\\) method")
  expect_null(out)

  registerS3method("residuals", "r2_throws",
                   function(object, ...) stop("engine lost its QR"))
  thrower <- fit; class(thrower) <- c("r2_throws", class(fit))
  expect_warning(residual_morans_i(thrower),
                 "residuals\\(\\) failed on this fit: engine lost its QR")

  registerS3method("residuals", "r2_short",
                   function(object, ...) rep(0.1, nrow(object$data_sf) - 1L))
  short <- fit; class(short) <- c("r2_short", class(fit))
  expect_warning(residual_morans_i(short),
                 "returned 99 value\\(s\\) for the 100 row\\(s\\)")
})


# ---------------------------------------------------------------------------
# first-pass FP6: a shared fold set that cannot be built is an error
# ---------------------------------------------------------------------------

test_that("compare_models_cv() errors when the requested blocks cannot be built", {
  skip_if_not_installed("ranger")
  d <- .r2e_pts(90, seed = 17)
  # One block covers everything.  It used to fall back to each backend's
  # default five-fold blocks, with no R condition.
  expect_error(suppressMessages(
    compare_models_cv(d, "z", "w", models = "RF", k = 3, block_size = 1e6,
                      rf_args = list(num_trees = 20), quiet = TRUE)),
    "compare_models_cv\\(\\): could not build the shared fold set")
})


# ---------------------------------------------------------------------------
# first-pass FP10: shared blocks placed by the caller's pointize
# ---------------------------------------------------------------------------

test_that("compare_models_cv() assigns polygon rows to blocks by `pointize`", {
  skip_if_not_installed("ranger")
  set.seed(3)
  n <- 80
  ell <- function(x0, y0, s) sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s / 4), c(x0 + s / 4, y0 + s / 4),
    c(x0 + s / 4, y0 + s), c(x0, y0 + s), c(x0, y0))))
  xs <- runif(n, 0, 900); ys <- runif(n, 0, 900); ss <- runif(n, 40, 120)
  d <- sf::st_sf(a = rnorm(n), geometry = sf::st_sfc(
    lapply(seq_len(n), function(i) ell(xs[i], ys[i], ss[i])), crs = 32632))
  d$z <- 2 * d$a + rnorm(n, 0, 0.3)
  fold_of <- function(ff) {
    out <- integer(n)
    for (j in seq_along(ff)) out[ff[[j]]$test] <- j
    out
  }
  for (pz in c("centroid", "auto")) {
    res <- suppressWarnings(suppressMessages(
      compare_models_cv(d, "z", "a", models = "RF", k = 3, pointize = pz,
                        rf_args = list(num_trees = 20), quiet = TRUE)))
    ref <- suppressWarnings(suppressMessages(
      make_folds(coerce_to_points(d, pz), k = 3, method = "block_kfold",
                 seed = 123)))
    expect_identical(fold_of(res$rf_cv$folds),
                     ref$assignment$fold[order(ref$assignment$row_id)],
                     info = pz)
    # The provenance probe is still taken on the polygons.
    expect_equal(res$rf_cv$n_folds_succeeded, 3L, info = pz)
  }
})


# ---------------------------------------------------------------------------
# gaps G4.6: select_features_forward()'s split indexes the layer as passed
# ---------------------------------------------------------------------------

test_that("select_features_forward()'s split positions index the caller's layer", {
  d <- .r2e_pts(200, seed = 18)
  d$a <- rnorm(200); d$b <- rnorm(200)
  d$z <- 2 * d$a + rnorm(200, 0, 0.5)
  d$b[c(3, 40, 77, 120, 150, 181, 199)] <- NA
  fs <- suppressWarnings(select_features_forward(
    d, "z", c("a", "b"), fit_fn = function(tr, v) lm_spatial_fit(tr, "z", v),
    k = 3, quiet = TRUE, select_on = "split"))
  both <- c(fs$split$selection, fs$split$estimation)
  complete <- which(!is.na(d$b))
  # Every complete row in exactly one half, and no dropped row in either.
  expect_setequal(both, complete)
  expect_false(anyDuplicated(both) > 0L)
  expect_false(anyNA(d$b[fs$split$estimation]))
})
