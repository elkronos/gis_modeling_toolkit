# select_features_forward()'s training-fold-mean baseline at a large offset
# (BUG-208).
#
# Where the learner cannot fit an empty predictor set, the sweep's baseline
# is the training-fold mean.  At an offset of 1e15 a double holds the
# response and its training mean to an eighth, so scored on the raw values
# the baseline read that rounding: RMSE 2.4737792 at 1e15 against 2.4701727
# at 0 on the layer below, while every candidate set scored the same at
# both.  The values are rounded as at 1e15 at every offset, so the stored
# spread is the same and every score must agree.  Scores are compared to a
# relative 1e-8: the centred sums are exact at both offsets up to the last
# bits of the platform's sum and sqrt.

fm_layer <- function(o, n = 300, seed = 3) {
  set.seed(seed)
  x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4); w <- rnorm(n); v <- rnorm(n)
  z <- 2 * w + sin(x / 2000) + cos(y / 1500) + rnorm(n, sd = 0.5)
  q <- function(a) (1e15 + a) - 1e15
  sf::st_as_sf(data.frame(x = x, y = y, z = o + q(z), p = q(w), v = q(v)),
               coords = c("x", "y"), crs = 32615)
}

# A GWR learner that refuses the empty set, so the sweep takes the fold mean.
# fit_gwr_model() needs sp, so every block that calls it skips without sp.
fm_gwr <- function(dd, v) {
  if (!length(v)) stop("no predictors")
  fit_gwr_model(dd, "z", v, bandwidth = 60)
}

fm_sweep <- function(o, fit_fn, metric = "RMSE", rule = "tol") {
  msgs <- character(0)
  s <- withCallingHandlers(
    suppressWarnings(select_features_forward(
      fm_layer(o), "z", c("p", "v"), fit_fn = fit_fn, k = 3,
      metric = metric, rule = rule, quiet = FALSE)),
    message = function(e) {
      msgs <<- c(msgs, conditionMessage(e))
      invokeRestart("muffleMessage")
    })
  fm <- grep("baseline is the training-fold mean", msgs, value = TRUE)
  list(s = s,
       figure = as.numeric(sub(".*\\(\\w+ = ([0-9.]+)\\).*", "\\1", fm[1L])))
}

fm_rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))

test_that("the fold-mean baseline of GWR folds scores the same at 0 and 1e15", {
  skip_on_cran()
  skip_if_not_installed("sp")
  base <- function(o) {
    pr <- suppressWarnings(suppressMessages(cv_spatial(
      fm_layer(o), "z", "p", fit_fn = function(d) fm_gwr(d, "p"),
      k = 3, seed = 1)))$predictions
    .cv_overall_metrics(.fs_fold_mean_predictions(list(pr)))
  }
  m0 <- base(0)
  expect_equal(m0$RMSE, 2.4701727, tolerance = 1e-7)
  expect_equal(m0$MAE, 1.9004936, tolerance = 1e-7)
  # Before the fix: RMSE 2.4737792 and MAE 1.9033333 at 1e15, and RMSE
  # 2.4701674 at 1e12.  fit_gwr_model() takes this response about a centre
  # from 1e4 on, and the baseline reads the training means about it.
  for (o in c(1e8, 1e12, 1e15)) {
    m <- base(o)
    expect_lt(fm_rel(m$RMSE, m0$RMSE), 1e-8, label = paste("RMSE at", o))
    expect_lt(fm_rel(m$MAE, m0$MAE), 1e-8, label = paste("MAE at", o))
    expect_identical(m$R2, 0)
  }
})

test_that("select_features_forward()'s fold-mean baseline does not move with the offset (GWR)", {
  skip_on_cran()
  skip_if_not_installed("sp")
  for (metric in c("RMSE", "MAE")) for (rule in c("tol", "se")) {
    a <- fm_sweep(0, fm_gwr, metric, rule)
    b <- fm_sweep(1e15, fm_gwr, metric, rule)
    lab <- paste(metric, rule)
    expect_identical(a$s$baseline$type, "fold mean", label = lab)
    expect_identical(b$s$baseline$type, "fold mean", label = lab)
    expect_lt(fm_rel(b$s$baseline$score, a$s$baseline$score), 1e-8, label = lab)
    h0 <- a$s$history; h15 <- b$s$history
    expect_lt(fm_rel(h15$score[h15$step == 0L], h0$score[h0$step == 0L]), 1e-8,
              label = lab)
    f0 <- a$s$fold_scores; f15 <- b$s$fold_scores
    s0 <- f0$score[f0$step == 0L]; s15 <- f15$score[f15$step == 0L]
    expect_identical(f15$fold[f15$step == 0L], f0$fold[f0$step == 0L], label = lab)
    expect_lt(fm_rel(s15, s0), 1e-8, label = lab)
    # The message prints four decimals.
    expect_true(is.finite(a$figure) && is.finite(b$figure), label = lab)
    expect_identical(b$figure, a$figure, label = lab)
    expect_identical(b$s$selected, a$s$selected, label = lab)
  }
})

test_that("a forest's fold-mean baseline does not move with the offset", {
  skip_on_cran()
  skip_if_not_installed("ranger")
  rf <- function(dd, v) fit_rf_model(dd, "z", v, num_trees = 50L, seed = 1L)
  a <- fm_sweep(0, rf)
  b <- fm_sweep(1e15, rf)
  expect_identical(a$s$baseline$type, "fold mean")
  expect_identical(b$s$baseline$type, "fold mean")
  # Before the fix the baseline at 1e15 read the training means rounded to
  # an eighth.  The forest's own candidate scores are its learner's and are
  # not compared.
  expect_lt(fm_rel(b$s$baseline$score, a$s$baseline$score), 1e-8)
  f0 <- a$s$fold_scores; f15 <- b$s$fold_scores
  expect_lt(fm_rel(f15$score[f15$step == 0L], f0$score[f0$step == 0L]), 1e-8)
  expect_identical(b$figure, a$figure)
})

fm_plain <- c("..row_id", "fold", "y", "y_train_mean", "yhat")

test_that("the baseline frame is the plain one where no GWR fold has a centre", {
  skip_on_cran()
  skip_if_not_installed("sp")
  pr <- suppressWarnings(suppressMessages(cv_spatial(
    fm_layer(0), "z", "p", fit_fn = function(d) fm_gwr(d, "p"),
    k = 3, seed = 1)))$predictions
  bp <- .fs_fold_mean_predictions(list(pr))
  expect_identical(names(bp), fm_plain)
  expect_identical(bp$yhat, bp$y_train_mean)
})

test_that("the baseline frame is the plain one where no lm fold has a centre", {
  skip_on_cran()
  # A fit with no centre, asked for the centre as the sweep asks: at 0 and
  # 1e8 the raw training response keeps its spread, so the folds record no
  # centre and the frame is the plain one.
  for (o in c(0, 1e8)) {
    .cv_request$fold_mean_centre <- TRUE
    pr <- suppressWarnings(suppressMessages(cv_spatial(
      fm_layer(o), "z", "p", fit_fn = function(d) lm_spatial_fit(d, "z", "p"),
      k = 3, seed = 1)))$predictions
    expect_false(any(.cv_mean_centre_cols %in% names(pr)))
    bp <- .fs_fold_mean_predictions(list(pr))
    expect_identical(names(bp), fm_plain)
    expect_identical(bp$yhat, bp$y_train_mean)
  }
})

test_that("cv_spatial() adds the centre columns only when the sweep asks, once", {
  skip_on_cran()
  lm_fit <- function(d) lm_spatial_fit(d, "z", "p")
  dd <- fm_layer(1e15)
  run <- function() suppressWarnings(suppressMessages(cv_spatial(
    dd, "z", "p", fit_fn = lm_fit, k = 3, seed = 1)))
  plain <- run()
  expect_false(any(.cv_mean_centre_cols %in% names(plain$predictions)))
  .cv_request$fold_mean_centre <- TRUE
  asked <- run()
  expect_null(.cv_request$fold_mean_centre)
  expect_true(all(.cv_mean_centre_cols %in% names(asked$predictions)))
  # The training mean about the centre holds what the raw mean rounds away.
  ap <- asked$predictions
  expect_true(all(ap[[.cv_mean_centre_cols[1L]]] != 0))
  # The request is used once: the next call is the plain one again.
  expect_identical(run(), plain)
  # Everything cv_spatial() returns but the two columns is as before.
  ap[.cv_mean_centre_cols] <- NULL
  expect_identical(ap, plain$predictions)
  expect_identical(asked$overall, plain$overall)
  expect_identical(asked$fold_metrics, plain$fold_metrics)
})

test_that("select_features_forward() leaves no request behind", {
  skip_on_cran()
  skip_if_not_installed("sp")
  fm_sweep(1e15, fm_gwr)
  expect_null(.cv_request$fold_mean_centre)
})

test_that("a sweep whose folds record a centre in some folds only runs", {
  skip_on_cran()
  # One row lies 2^17 above the rest.  A training set that holds it keeps
  # its spread at 1e15 (sd / max about 8e-12, above the 1e-12 of
  # .regression_centre()), and the one training set without it does not
  # (about 2e-15), so one fold records the centre columns and the others
  # carry NA in them once cv_spatial() binds the folds.  The margins are
  # wide, so the split is the same on every platform.
  mixed <- function(o) {
    d <- fm_layer(o)
    d$z[1L] <- o + 2^17
    d
  }
  lm_refuse <- function(dd, v) {
    if (!length(v)) stop("no predictors")
    lm_spatial_fit(dd, "z", v)
  }
  .cv_request$fold_mean_centre <- TRUE
  pr <- suppressWarnings(suppressMessages(cv_spatial(
    mixed(1e15), "z", "p", fit_fn = function(d) lm_refuse(d, "p"),
    k = 5, seed = 1)))$predictions
  centred <- tapply(is.finite(pr[[.cv_mean_centre_cols[1L]]]), pr$fold, all)
  expect_identical(sum(centred), 1L)
  expect_true(all(is.na(pr[[.cv_mean_centre_cols[1L]]][
    pr$fold %in% names(centred)[!centred]])))
  bp <- .fs_fold_mean_predictions(list(pr))
  expect_s3_class(bp, "data.frame")
  expect_identical(nrow(bp), nrow(pr))
  on <- bp$fold %in% names(centred)[centred]
  expect_true(all(bp$y_centre[on] != 0))
  expect_true(all(bp$y_centre[!on] == 0))
  expect_true(all(is.finite(bp$yhat_centred)))
  expect_true(all(is.finite(unlist(.cv_overall_metrics(bp)[c("RMSE", "MAE")]))))

  # The sweep takes the same frame outside its try(), so a frame it cannot
  # build stops the sweep.  Its baseline is scored: the uncentred folds
  # read their training means, which a double holds at 1e15 to an eighth,
  # so the score is within 0.125 of the score at 0.
  sw <- function(o) suppressWarnings(suppressMessages(select_features_forward(
    mixed(o), "z", c("p", "v"), fit_fn = lm_refuse, k = 5, seed = 1,
    quiet = TRUE)))
  a <- sw(0); b <- sw(1e15)
  expect_identical(b$baseline$type, "fold mean")
  expect_true(is.finite(b$baseline$score))
  expect_lte(abs(b$baseline$score - a$baseline$score), 0.125)
  expect_null(.cv_request$fold_mean_centre)
})
